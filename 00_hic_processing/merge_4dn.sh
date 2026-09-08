#!/bin/bash
#SBATCH --job-name=merge4dn
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-v100
#SBATCH --cpus-per-task=20
#SBATCH --mem=110g
#SBATCH --time=48:00:00
#SBATCH --mail-user=fanfeng@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_4DN_merged/logs/%x_%j.log
# Merge deduplicated 4DN valid-pairs across libraries of the SAME sample, re-dedup, rebuild .hic.
# Usage: sbatch -J <SAMPLE>_merge merge_4dn.sh <SAMPLE> <OUTDIR> <pairs1.gz> <pairs2.gz> [pairs3 ...]
set -eo pipefail
SAMPLE="$1"; OUT="$2"; shift 2; INPUTS=("$@"); CPU=20
CHROM=/home/fanfeng/HiCPro/envs3/HiC-Pro_3.1.0/annotation/chrom_hg19.sizes
JUICER_JAR=/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar
W=/scratch/drjieliu_root/drjieliu/fanfeng/temp/${SAMPLE}_merge
source /home/fanfeng/anaconda3/etc/profile.d/conda.sh
rm -rf "$W"; mkdir -p "$W/tmp"; cd "$W"; mkdir -p "$OUT"
echo "=== $SAMPLE MERGE start $(date) on $(hostname) ==="
printf 'input: %s\n' "${INPUTS[@]}"
conda activate hic_snv
set -u
echo "=== pairtools merge $(date) ==="
pairtools merge --nproc $CPU --tmpdir "$W/tmp" -o merged.pairs.gz "${INPUTS[@]}"
echo "=== pairtools dedup (across libraries) $(date) ==="
pairtools dedup --mark-dups --output-stats dedup.stats.txt -o dedup.pairs.gz merged.pairs.gz
rm -f merged.pairs.gz
echo "=== select valid UU/UR/RU $(date) ==="
pairtools select "(pair_type=='UU') or (pair_type=='UR') or (pair_type=='RU')" dedup.pairs.gz -o valid.pairs.gz
rm -f dedup.pairs.gz
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
cp valid.pairs.gz valid.pairs.gz.px2 "${SAMPLE}.hic" dedup.stats.txt valid.stats.txt "$OUT/" 2>/dev/null || true
echo "=== $SAMPLE MERGE SUMMARY ==="
grep -E "^total\b|^total_dups|^total_nodups|^cis\b|^trans\b" dedup.stats.txt
echo "=== $SAMPLE MERGE done $(date) ==="
