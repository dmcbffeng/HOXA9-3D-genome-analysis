#!/bin/bash
#SBATCH --job-name=insul
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-v100
#SBATCH --cpus-per-task=2
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH --output=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/ins_scores/slurm-%x-%j.out

cd /nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/ins_scores
python insulation.py "$1" "$2"
