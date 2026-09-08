#!/bin/bash
#SBATCH --job-name=O3_4DN_full
#SBATCH --account=drjieliu
#SBATCH --partition=drjieliu-v100
#SBATCH --cpus-per-task=20
#SBATCH --mem=110g
#SBATCH --time=96:00:00
#SBATCH --mail-user=fanfeng@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/drjieliu_root/drjieliu/fanfeng/temp/O3_24h_4DN/full_4dn_%j.log

# ---------------------------------------------------------------------------
# 4DN-method Hi-C processing for DW_HiC_24h (hg19, DpnII/MboI)
#   bwa mem -SP5M  ->  pairtools parse2 (5unique walks)  ->  sort -> dedup
#   -> select valid -> cooler .mcool (balanced) -> juicer .hic (best-effort)
#   Mirrors the 4DN dcic hi_c-processing-pipeline, run locally.
# ---------------------------------------------------------------------------
set -eo pipefail
CPU=20
SRC=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_24h/01.RawData/DW_HiC_24h
R1="$SRC/DW_HiC_24h_CKDL260012172-1A_23JK7LLT3_L5_1.fq.gz"
R2="$SRC/DW_HiC_24h_CKDL260012172-1A_23JK7LLT3_L5_2.fq.gz"
IDX=/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/MOLM13/bwaindex/hg19
# hg19 chrom.sizes, 2 columns (chr1-22, X, Y, M). A copy ships alongside this
# script as hg19.chrom.sizes -- point CHROM at it, or export CHROM to override.
CHROM="${CHROM:-/path/to/hg19.chrom.sizes}"
JUICER_JAR=/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar
W=/scratch/drjieliu_root/drjieliu/fanfeng/temp/O3_24h_4DN/full
RES=/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis/HiC/O3_24h/04.4DN_alignment
source /home/fanfeng/anaconda3/etc/profile.d/conda.sh

rm -rf "$W"; mkdir -p "$W/tmp"; cd "$W"
echo "Job start $(date) on $(hostname)"

conda activate hic_snv
set -u
echo "=== bwa mem -SP5M | pairtools parse2 $(date) ==="
bwa mem -SP5M -t $CPU "$IDX" "$R1" "$R2" 2> bwa.log | \
  pairtools parse \
    --min-mapq 10 --max-inter-align-gap 30 --walks-policy 5unique \
    --drop-seq --drop-sam \
    --chroms-path "$CHROM" --assembly hg19 \
    --nproc-in $CPU --nproc-out $CPU \
    -o parsed.pairs.gz

echo "=== sort $(date) ==="
pairtools sort --nproc $CPU --tmpdir "$W/tmp" -o sorted.pairs.gz parsed.pairs.gz
rm -f parsed.pairs.gz

echo "=== dedup $(date) ==="
pairtools dedup --mark-dups --output-stats dedup.stats.txt -o dedup.pairs.gz sorted.pairs.gz
rm -f sorted.pairs.gz

echo "=== select valid (UU+UR+RU) $(date) ==="
pairtools select "(pair_type=='UU') or (pair_type=='UR') or (pair_type=='RU')" dedup.pairs.gz -o valid.pairs.gz
bgzip -@ $CPU -c /dev/null >/dev/null 2>&1 || true
pairtools stats valid.pairs.gz -o valid.stats.txt
# index for random access
pairix -f valid.pairs.gz 2>/dev/null || true

echo "=== cooler .mcool (balanced, multi-res) $(date) ==="
conda activate fithichip2   # cooler 0.9.1 + pandas 1.3.5 (the 'cooler' env's pandas 2.0 breaks cload)
BASE=1000
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly hg19 "$CHROM":$BASE valid.$BASE.cool
cooler zoomify -p $CPU --balance --resolutions 1000,2000,5000,10000,20000,40000,100000,500000,1000000 \
  -o DW_HiC_24h.mcool valid.$BASE.cool

echo "=== juicer .hic (best-effort) $(date) ==="
conda activate hic_snv
# pairs -> juicer short format (str chr pos frag), sorted by chr1 chr2
zcat valid.pairs.gz | grep -v '^#' | \
  awk 'BEGIN{OFS="\t"}{s1=($6=="+")?0:1; s2=($7=="+")?0:1; print s1,$2,$3,0,s2,$4,$5,1}' > juicer_short.txt
sort -k2,2 -k6,6 --parallel=$CPU -S 20G juicer_short.txt > juicer_short.sorted.txt
java -Xmx90g -jar "$JUICER_JAR" pre juicer_short.sorted.txt DW_HiC_24h.hic "$CHROM" \
  && rm -f juicer_short.txt juicer_short.sorted.txt \
  || echo "WARNING: juicer pre failed (non-fatal); .mcool is the primary deliverable"

echo "=== copy results to turbo $(date) ==="
mkdir -p "$RES/results"
cp valid.pairs.gz valid.pairs.gz.px2 "$RES/results/" 2>/dev/null || true
cp DW_HiC_24h.mcool "$RES/results/" 2>/dev/null || true
cp DW_HiC_24h.hic "$RES/results/" 2>/dev/null || true
cp dedup.stats.txt valid.stats.txt bwa.log "$RES/results/" 2>/dev/null || true

echo "############## 4DN FULL SUMMARY ##############"
echo "----- dedup stats (key) -----"
grep -E "^total\b|^total_mapped|^total_unmapped|^total_single_sided|^total_dups|^total_nodups|^cis\b|^trans\b" dedup.stats.txt 2>/dev/null
VALID=$(zcat valid.pairs.gz | grep -cv '^#')
echo "valid_pairs (UU+UR+RU, deduped): $VALID"
echo "Job done $(date)"
