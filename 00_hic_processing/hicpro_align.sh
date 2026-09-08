#!/bin/bash
#SBATCH --job-name=HiCPro_O3_24h
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-v100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120g
#SBATCH --time=96:00:00
#SBATCH --mail-user=fanfeng@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=%x_%j.log

set -eo pipefail

# ---------------------------------------------------------------------------
# HiC-Pro alignment for HiC/O3_24h  (sample DW_HiC_24h)
#   Genome : hg19  (bowtie2 index at MOLM13/bowtie2index/hg19)
#   Enzyme : MboI/DpnII (GATC), ligation junction GATCGATC  -- determined
#            empirically (~38% of reads carry GATCGATC; HindIII absent)
#   Temp/working area : /scratch  (large intermediates stay here)
#   Final results     : copied back to this folder on /nfs/turbo
# ---------------------------------------------------------------------------

RESULTS=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_24h/03.HiCPro_alignment
WORK=/scratch/drjieliu_root/drjieliu/fanfeng/temp/O3_24h_HiCPro
OUT=$WORK/output
CONFIG=$RESULTS/config-hicpro.txt
CHROMSIZES=/home/fanfeng/HiCPro/envs3/HiC-Pro_3.1.0/annotation/chrom_hg19.sizes
JUICER_JAR=/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar

echo "Job started : $(date)"
echo "Host        : $(hostname)"

source /home/fanfeng/anaconda3/etc/profile.d/conda.sh
conda activate fithichip2
set -u   # enable only after conda activate (its activate.d trips on unbound vars)

# HiC-Pro requires the output dir to NOT pre-exist
rm -rf "$OUT"

echo "===== Running HiC-Pro ====="
HiC-Pro -i "$WORK/rawdata" -o "$OUT" -c "$CONFIG"
echo "===== HiC-Pro finished : $(date) ====="

# --- Generate a .hic (Juicebox) from all valid pairs (best-effort) ----------
VP="$OUT/hic_results/data/DW_HiC_24h/DW_HiC_24h.allValidPairs"
if [[ -f "$VP" ]]; then
    echo "===== hicpro2juicebox : building .hic ====="
    /home/fanfeng/HiCPro/envs3/HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh \
        -i "$VP" -g "$CHROMSIZES" -j "$JUICER_JAR" -t "$WORK/juicer_tmp" \
        || echo "WARNING: hicpro2juicebox failed (non-fatal)"
fi

# --- Copy final results back to turbo ---------------------------------------
echo "===== Copying results to $RESULTS ====="
mkdir -p "$RESULTS/results"
# QC reports, stats and plots (small, the key deliverable)
cp -r "$OUT/hic_results/pic"    "$RESULTS/results/" 2>/dev/null || true
cp -r "$OUT/hic_results/stats"  "$RESULTS/results/" 2>/dev/null || true
# Valid pairs (main alignment output) and contact matrices
cp -r "$OUT/hic_results/data"   "$RESULTS/results/" 2>/dev/null || true
cp -r "$OUT/hic_results/matrix" "$RESULTS/results/" 2>/dev/null || true
# Logs
cp -r "$OUT/logs" "$RESULTS/results/logs" 2>/dev/null || true
# .hic if produced
cp "$OUT/hic_results/data/DW_HiC_24h/"*.hic "$RESULTS/results/" 2>/dev/null || true

echo "Job finished : $(date)"
echo "Intermediate BAMs / tmp remain under: $WORK"
