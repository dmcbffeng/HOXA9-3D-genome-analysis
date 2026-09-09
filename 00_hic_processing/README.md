# Hi-C processing (4DN pipeline)

Methods §"Hi-C library preparation, Hi-C data processing and embedding in 2D/3D space".
Upstream of every other directory. QC metrics from these runs are table S8.

Hi-C reads are aligned with the 4DN DCIC processing pipeline
(https://github.com/4dn-dcic/docker-4dn-hic), run locally rather than through CWL:

    bwa mem -SP5M  ->  pairtools parse  ->  sort  ->  dedup  ->  select valid
                   ->  cooler (.mcool, balanced)  and  juicer pre (.hic)

`bwa mem -SP5M` is the 4DN setting for Hi-C: `-S`/`-P` skip mate rescue and pairing so each
side is placed independently, and `-5M` picks the 5'-most alignment as primary so chimeric
reads spanning a ligation junction are assigned to the correct fragment end. Chimera
handling therefore happens at the alignment step, and no restriction-fragment file or
ligation-junction sequence is required.

| Script | What it does |
|---|---|
| `fourdn_align.sh` | FASTQ → `.hic`/`.mcool` for one library. `bwa mem -SP5M` → `pairtools parse` (`--min-mapq 10 --max-inter-align-gap 30 --walks-policy 5unique`) → `pairtools sort` → `pairtools dedup` → `pairtools select` keeping `UU`/`UR`/`RU` pairs → `cooler cload` at 1 kb + `cooler zoomify --balance` (multi-resolution balanced `.mcool`) → `juicer pre` (`.hic`). |
| `merge_4dn.sh` | Merge deduplicated valid-pairs across libraries of the *same* biological sample, re-dedup across libraries so cross-library PCR duplicates are removed, rebuild one `.hic`. Produced the merged HOXA9-KO (348M valid pairs) and 24 h dTAG-13 (595M valid pairs) maps. |
| `downsample_4dn.sh` | `pairtools sample` at a fixed fraction and seed, then the same `.hic` build as `merge_4dn.sh`, so the downsampled map is directly comparable. Used to match dTAG-13 depth to OCI-AML3 WT as a control for the loop and APA comparisons. |
| `hic2cool.py` | `.hic` → per-chromosome `.cool` at 50 kb via HiCExplorer `hicConvertFormat`. Input step for `hicrep_scc.py`. |
| `hicrep_scc.py` | Pairwise HiCRep stratum-adjusted correlation between contact maps (dejunlin/hicrep), per chromosome, `h = 1`, `dBPMax = 2 Mb`. Writes one `chr*.csv` similarity matrix per chromosome. Distances from these matrices were embedded with scHiCTools (`liu-bioinfo-lab/scHiCTools`). |
| `hg19.chrom.sizes` | hg19 primary chromosomes (chr1–22, X, Y, M), used by `pairtools parse`, `cooler` and `juicer pre`. Set `CHROM` in the three shell scripts to its path. |

All three shell scripts stage large intermediates on node-local scratch and copy only the
final valid-pairs / `.mcool` / `.hic` back. Edit the path constants at the top (`R1`, `R2`,
`IDX`, `CHROM`, `JUICER_JAR`, `W`, `RES`) before use.

## Scope

These scripts cover the contact maps aligned for this study. Contact maps generated earlier
and reused here were aligned separately and are not reproduced by this pipeline.

The FitHiChIP entry point in `../03_loops/` is named `FitHiChIP_HiCPro.sh` and its config
files carry upstream comments referring to HiC-Pro inputs. That is FitHiChIP's own naming
for its generic driver; the runs here supply contact data through the `Bed=` and `HIC=`
options, not through HiC-Pro output.
