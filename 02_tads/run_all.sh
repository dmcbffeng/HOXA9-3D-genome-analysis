#!/usr/bin/env bash
# Convenience driver. Usage:
#   ./run_all.sh                       # run every sample in samples.tsv
#   ./run_all.sh 4943 MOLM13           # run a subset
#
# Outputs:
#   temp_matrix/<sample>_<chr>.matrix       (TopDom-format)
#   results/<sample>_<chr>.{bed,binSignal,domain}
#   merged/<sample>.bed                     (concatenated TAD calls)
#   logs/factor_<sample>.txt                (auto-detected common factor)

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PY="${PY:-python3}"
# R: prefer the lmod-managed R if available, otherwise the fithichip conda env.
if command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="Rscript"
elif [[ -x /home/fanfeng/anaconda3/envs/fithichip/bin/Rscript ]]; then
    RSCRIPT="/home/fanfeng/anaconda3/envs/fithichip/bin/Rscript"
else
    echo "Could not find Rscript. module load R/4.4.3 first or set RSCRIPT env var." >&2
    exit 1
fi

SAMPLE_ARGS=()
if [[ $# -gt 0 ]]; then
    SAMPLE_ARGS=(--samples "$@")
fi

echo "[1/4] Generating temp matrices (with factor auto-detection) ..."
"$PY" 01_process_matrix.py "${SAMPLE_ARGS[@]}" 2>&1 | tee -a logs/01_process_matrix.log

echo "[2/4] Running TopDom ..."
"$RSCRIPT" 02_run_TopDom.R "${SAMPLE_ARGS[@]}" 2>&1 | tee -a logs/02_run_TopDom.log

echo "[3/4] Checking for empty matrices / missing beds ..."
"$PY" 03_check_empty.py 2>&1 | tee -a logs/03_check_empty.log

echo "[4/4] Merging per-chromosome beds into per-sample beds ..."
"$PY" 04_merge.py 2>&1 | tee -a logs/04_merge.log

echo "Done. Per-sample TAD bed files are in: $HERE/merged/"
