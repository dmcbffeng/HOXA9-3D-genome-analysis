#!/bin/bash
#SBATCH --job-name=ds4dn
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-v100
#SBATCH --cpus-per-task=20
#SBATCH --mem=110g
#SBATCH --time=48:00:00
#SBATCH --mail-user=fanfeng@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_4DN_downsampled_O3WT/logs/%x_%j.log
# Uniformly downsample a 4DN valid-pairs file by FRAC (matching cis depth to O3-WT), rebuild .hic.
# Mirrors merge_4dn.sh exactly for the .hic build so the map is directly comparable.
# Usage: sbatch -J <SAMPLE>_ds downsample_4dn.sh <SAMPLE> <FRAC> <IN_pairs.gz> <OUTDIR>
set -eo pipefail
SAMPLE="$1"; FRAC="$2"; INPAIRS="$3"; OUT="$4"; CPU=20; SEED=1
CHROM=/home/fanfeng/HiCPro/envs3/HiC-Pro_3.1.0/annotation/chrom_hg19.sizes
JUICER_JAR=/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar
W=/scratch/drjieliu_root/drjieliu/fanfeng/temp/${SAMPLE}_ds
source /home/fanfeng/anaconda3/etc/profile.d/conda.sh
rm -rf "$W"; mkdir -p "$W/tmp"; cd "$W"; mkdir -p "$OUT"
echo "=== $SAMPLE DOWNSAMPLE frac=$FRAC seed=$SEED start $(date) on $(hostname) ==="
echo "input: $INPAIRS"
conda activate hic_snv
set -u
echo "=== pairtools sample $(date) ==="
pairtools sample --seed $SEED --nproc-in $CPU --nproc-out $CPU -o valid.pairs.gz "$FRAC" "$INPAIRS"
pairtools stats valid.pairs.gz -o valid.stats.txt
pairix -f valid.pairs.gz 2>/dev/null || true
echo "=== juicer .hic $(date) ==="
zcat valid.pairs.gz | grep -v '^#' | \
  awk -v cf="$CHROM" 'BEGIN{OFS="\t"; while((getline l < cf)>0){split(l,a,"\t"); ok[a[1]]=1}}
       ($2 in ok)&&($4 in ok){s1=($6=="+")?0:1; s2=($7=="+")?0:1; print s1,$2,$3,0,s2,$4,$5,1}' > juicer_short.txt
LC_ALL=C sort -k2,2 -k6,6 -S 60G --parallel=$CPU -T "$W/tmp" juicer_short.txt > juicer_short.sorted.txt
rm -f juicer_short.txt
/usr/bin/java -Xmx100g -jar "$JUICER_JAR" pre \
  -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 \
  juicer_short.sorted.txt "${SAMPLE}.hic" "$CHROM"
rm -f juicer_short.sorted.txt
echo "=== copy to turbo $(date) ==="
cp valid.pairs.gz "${SAMPLE}.hic" valid.stats.txt "$OUT/" 2>/dev/null || true
[ -f valid.pairs.gz.px2 ] && cp valid.pairs.gz.px2 "$OUT/" 2>/dev/null || true
echo "=== $SAMPLE stats ==="
grep -E "^total\b|^cis\b|^trans\b|cis_20kb" valid.stats.txt || true
echo "=== $SAMPLE DOWNSAMPLE done $(date) ==="
