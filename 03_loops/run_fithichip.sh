#!/bin/bash
#SBATCH --job-name=FitHiChIP
#SBATCH --mem=128g
#SBATCH --time=96:00:00
#SBATCH --output=slurm-%x-%j.out
#
# Call HOXA9-anchored chromatin loops with FitHiChIP.
#
#   sbatch -J FH_O3_WT run_fithichip.sh configs/configfile_O3_WT_HOXA9_10000_0.05
#   ./run_fithichip.sh configs/configfile_AML4943_HOXA9_10000_0.05   # or run directly
#
# FitHiChIP itself is the unmodified upstream release
# (https://github.com/ay-lab/FitHiChIP); set FITHICHIP_HOME to its checkout.
# Every parameter that differs from the FitHiChIP defaults lives in the config
# files under configs/ -- see 03_loops/README.md for the table.

set -euo pipefail

CONFIG="${1:?usage: run_fithichip.sh <configfile>}"
FITHICHIP_HOME="${FITHICHIP_HOME:-$HOME/FitHiChIP}"

# conda env with FitHiChIP's R/python dependencies (see ../README.md)
if [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate "${FITHICHIP_ENV:-fithichip}"
fi

sh "$FITHICHIP_HOME/FitHiChIP_HiCPro.sh" -C "$CONFIG"
