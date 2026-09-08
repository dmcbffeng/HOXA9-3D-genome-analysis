# Hi-C processing

Methods §"Hi-C library preparation, Hi-C data processing and embedding in 2D/3D space".
Upstream of every other directory. QC metrics from these runs are table S8.

| Script | What it does |
|---|---|
| `hicpro_align.sh` + `config-hicpro.txt` | HiC-Pro 3.1.0, hg19, MboI/DpnII (`GATC`, junction `GATCGATC`), 20 CPUs. Produces `allValidPairs`, then `hicpro2juicebox.sh` → `.hic`. This is the pipeline used for the OCI-AML3, OCI-AML3-HOXA9-KO, HSPC and AML4943 maps. |
| `fourdn_align.sh` | 4DN-style alternative: `bwa mem -SP5M` → `pairtools parse` (`--min-mapq 10 --walks-policy 5unique`) → sort → dedup → select `UU/UR/RU` → `cooler zoomify` (balanced `.mcool`) → `juicer pre` (`.hic`). Used for the 24 h dTAG-13 library — see the note below. |
| `merge_4dn.sh` | Merge deduplicated valid-pairs across libraries of the *same* biological sample, re-dedup across libraries, rebuild one `.hic`. Produced the merged HOXA9-KO (348M valid pairs) and 24 h dTAG-13 (595M) maps. |
| `downsample_4dn.sh` | `pairtools sample` at a fixed fraction and seed, then the same `.hic` build as `merge_4dn.sh`, so the downsampled map is directly comparable. Used to match dTAG-13 depth to OCI-AML3 WT as a control for the loop and APA comparisons. |
| `hic2cool.py` | `.hic` -> per-chromosome `.cool` at 50 kb via HiCExplorer `hicConvertFormat`. Input step for `hicrep_scc.py`. |
| `hicrep_scc.py` | Pairwise HiCRep stratum-adjusted correlation between contact maps (dejunlin/hicrep), per chromosome, `h = 1`, `dBPMax = 2 Mb`. Writes one `chr*.csv` similarity matrix per chromosome. As run, this covered the OCI-AML3 / HSPC maps together with the published primary-AML and HSPC cohort maps used for the AML-specific TAD comparison in fig. S1B. Distances from these matrices were embedded with scHiCTools (`liu-bioinfo-lab/scHiCTools`). |

Both alignment scripts stage large intermediates on node-local scratch and copy only the
final valid-pairs / `.mcool` / `.hic` back. Edit the path constants at the top before use.

## Why two alignment pipelines

The OCI-AML3, OCI-AML3-HOXA9-KO, CD34+ HSPC and AML4943 maps were produced with HiC-Pro.
The 24 h dTAG-13 library, and the merged HOXA9-KO map, were produced with the 4DN pipeline
instead: HiC-Pro's end-to-end bowtie2 step reports ~0.11% valid pairs on that library, which
is an artifact of end-to-end mapping on it and not a bad library — the same FASTQs give ~76%
valid pairs through `bwa mem -SP5M` + pairtools. Both pipelines write hg19 `.hic` files that
are used interchangeably downstream, and every downstream analysis here is O/E-normalised,
so the two are directly comparable.
