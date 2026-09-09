# TAD calling (TopDom)

Fig. 1A–C, figs. S1C/S1G/S4C, table S2. Methods §"Compartment, TAD calling and analysis of
TAD boundary".

TopDom v0.0.2 at **50 kb**, `window.size = 10` (500 kb), `statFilter = TRUE`, chr1–22 + X.
Only rows annotated `domain` are TADs; `boundary` and `gap` rows are kept in the output for
the boundary-comparison step.

| Script | What it does |
|---|---|
| `samples.tsv` | The sample table — `sample_name<TAB>hic_path`, one row per contact map. This is the only file to edit to add a sample. |
| `01_process_matrix.py` | `juicer_tools dump oe KR` per chromosome → TopDom-format dense matrix in `temp_matrix/<sample>_<chr>.matrix`. |
| `02_run_TopDom.R` | TopDom v0.0.2 (source inlined, unmodified, credited in the header) over every temp matrix → `results/<sample>_<chr>.{bed,binSignal,domain}`. |
| `03_check_empty.py` | QC — reports empty/missing matrices and chromosomes that produced no BED. |
| `04_merge.py` | Concatenates per-chromosome BEDs → `merged/<sample>.bed`. |
| `run_all.sh` | Runs all four steps. `./run_all.sh` for every sample, `./run_all.sh OCI-AML3 HSPC` for a subset. |

**Why O/E:** two of the deposited maps (`4943.hic`, `O3.hic`) have raw observed counts that
are exact multiples of 60. Using `oe` makes that factor cancel within each distance
stratum, so no rescaling is needed. If you switch this pipeline to `observed NONE` or
`observed KR`, divide those two maps by 60 first.

**Not included here:** classifying each TAD as unchanged / fused / split / shifted
(BEDTools `windowBed` on the boundary calls this pipeline produces).
