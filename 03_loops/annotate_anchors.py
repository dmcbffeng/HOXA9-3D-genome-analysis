#!/usr/bin/env python
"""
Annotate FitHiChIP interaction anchors (merged calls only) for the runs in
FitHiChIP_202608_newO3.

Two OCI-AML3 HOXA9 CUT&RUN FitHiChIP runs, each on its own new 4DN-merged Hi-C
contact map (used directly via HIC=), in two conditions:
  * O3_HOXA9KO_HOXA9     : HOXA9-KO Hi-C (O3_HOXA9KO.hic)
  * O3_24h_dTAG13_HOXA9  : 24 h dTAG-13 Hi-C (O3_24h_dTAG13.hic)

Both use the same HOXA9 peak set (peaks/HOXA9_peaks.bed, the 2,331-peak
O3_HOXA9_HA_overlap_peaks_IgGcorr.bed).

Each interaction has two anchors: (chr1,s1,e1) and (chr2,s2,e2), each exactly
one 10 kb bin. For every anchor we append Yes/No labels (plus gene names where
relevant) for four features:

  1. HOXA9_peak   : anchor bin overlaps a HOXA9 peak.
  2. TAD_boundary : anchor bin overlaps a TopDom 'boundary' region of the
                    matching condition's own merged TAD calls
                    (O3_HOXA9KO_TAD_TopDom_50kb.bed for the KO run,
                    O3_24h_dTAG13_TAD_TopDom_50kb.bed for the 24h run).
  3. promoter     : anchor bin overlaps the TSS +/- 2 kb window of >=1 gene;
                    the promoter gene name(s) are recorded.
  4. gene_body    : anchor bin overlaps the body of >=1 gene *and is not already
                    a promoter of that same gene*; the gene name(s) are recorded.

Gene coordinates / strand / TSS come from the shared table written by
resolve_gene_strand_tss.py (all_genes/hg19_genes_from_David_with_strand.tsv).

Overlap is evaluated on the 10 kb bin grid: every feature interval is indexed
into each 10 kb bin it spans, and an anchor (one bin) looks up its bin index.

Output: gene_annotation/<run>_merged.annotated.csv (comma-delimited; gene
lists within a field are ';'-separated so no quoting is needed).
"""
import os
import csv

RES = 10000
PROMOTER_FLANK = 2000  # +/- around TSS
CHROMS = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]
BASE = '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis'

GENE_STRAND_TABLE = f'{BASE}/all_genes/hg19_genes_from_David_with_strand.tsv'
# Per-condition merged TAD boundary calls (TopDom 50 kb), alongside each .hic
TAD_BED_HOXA9KO = f'{BASE}/HiC/O3_4DN_merged/O3_HOXA9KO/O3_HOXA9KO_TAD_TopDom_50kb.bed'
TAD_BED_24H = f'{BASE}/HiC/O3_4DN_merged/O3_24h_dTAG13/O3_24h_dTAG13_TAD_TopDom_50kb.bed'
RUN_DIR = f'{BASE}/FitHiChIP_202608_newO3'
PEAK_DIR = f'{RUN_DIR}/peaks'
OUT_DIR = f'{RUN_DIR}/gene_annotation'


# ---------- bin index ---------------------------------------------------------

def new_index():
    return {ch: {} for ch in CHROMS}


def add_interval(idx, ch, start, end, label):
    if ch not in idx:
        return
    for b in range(start // RES, end // RES + 1):
        idx[ch].setdefault(b, set()).add(label)


def lookup(idx, ch, start):
    return idx.get(ch, {}).get(start // RES, set())


# ---------- feature loaders ---------------------------------------------------

def load_gene_indices():
    """Build promoter and gene-body bin indices from the strand-annotated gene
    table written by resolve_gene_strand_tss.py. Gene body uses [start,end];
    promoter uses TSS +/- PROMOTER_FLANK (only for genes with a known strand).
    Returns (prom_idx, body_idx, stats)."""
    if not os.path.exists(GENE_STRAND_TABLE):
        raise SystemExit(
            f'Missing {GENE_STRAND_TABLE}\n'
            f'Run resolve_gene_strand_tss.py first.')
    prom_idx, body_idx = new_index(), new_index()
    n_total = n_resolved = 0
    with open(GENE_STRAND_TABLE) as f:
        next(f)  # header
        for line in f:
            name, ch, start, end, strand, tss = line.rstrip('\n').split('\t')
            if ch not in prom_idx:
                continue
            n_total += 1
            start, end = int(start), int(end)
            add_interval(body_idx, ch, start, end, name)
            if strand in ('+', '-'):
                n_resolved += 1
                tss = int(tss)
                add_interval(prom_idx, ch, max(0, tss - PROMOTER_FLANK),
                             tss + PROMOTER_FLANK, name)
    stats = dict(total=n_total, resolved=n_resolved,
                 unresolved=n_total - n_resolved)
    return prom_idx, body_idx, stats


def load_boundary_index(tad_bed):
    """TopDom 'boundary' rows from a TAD bed."""
    idx = new_index()
    for line in open(tad_bed):
        c = line.rstrip('\n').split('\t')
        if len(c) < 4 or c[3] != 'boundary':
            continue
        add_interval(idx, c[0], int(c[1]), int(c[2]), 'boundary')
    return idx


def load_peak_index(peak_file):
    idx = new_index()
    for line in open(peak_file):
        c = line.rstrip('\n').split('\t')
        add_interval(idx, c[0], int(c[1]), int(c[2]), c[3])
    return idx


# ---------- annotation --------------------------------------------------------

def anchor_labels(ch, start, peak_idx, bnd_idx, prom_idx, body_idx):
    peak = 'Yes' if lookup(peak_idx, ch, start) else 'No'
    tad = 'Yes' if lookup(bnd_idx, ch, start) else 'No'
    prom_genes = lookup(prom_idx, ch, start)
    # a gene that is a promoter here is not also counted as gene body here
    body_genes = lookup(body_idx, ch, start) - prom_genes
    return [
        peak,
        tad,
        'Yes' if prom_genes else 'No',
        ';'.join(sorted(prom_genes)),
        'Yes' if body_genes else 'No',
        ';'.join(sorted(body_genes)),
    ]


def annotate(in_file, out_file, tf_name, peak_idx, bnd_idx, prom_idx, body_idx):
    def cols(n):
        return [f'Anchor{n}_{tf_name}_peak', f'Anchor{n}_TAD_boundary',
                f'Anchor{n}_promoter', f'Anchor{n}_promoter_genes',
                f'Anchor{n}_gene_body', f'Anchor{n}_gene_body_genes']

    with open(in_file) as fin, open(out_file, 'w', newline='') as fout:
        writer = csv.writer(fout)
        header = fin.readline().rstrip('\n').split('\t')
        writer.writerow(header + cols(1) + cols(2))
        n = 0
        for line in fin:
            lst = line.rstrip('\n').split('\t')
            ch1, s1 = lst[0], int(lst[1])
            ch2, s2 = lst[3], int(lst[4])
            row = (lst
                   + anchor_labels(ch1, s1, peak_idx, bnd_idx, prom_idx, body_idx)
                   + anchor_labels(ch2, s2, peak_idx, bnd_idx, prom_idx, body_idx))
            writer.writerow(row)
            n += 1
    return n


# ---------- driver ------------------------------------------------------------

MERGED_REL = ('FitHiChIP_Peak2ALL_b10000_L20000_U3000000/P2PBckgr_0/'
              'Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/'
              'FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed')

# run/folder stem -> (TF column name, TAD bed). Peak file is peaks/HOXA9_peaks.bed,
# merged input is <stem>_10000_0.05/<MERGED_REL>, output <stem>_merged.annotated.csv.
RUNS = [
    ('O3_HOXA9KO_HOXA9',    'HOXA9', TAD_BED_HOXA9KO),
    ('O3_24h_dTAG13_HOXA9', 'HOXA9', TAD_BED_24H),
]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print('Building gene promoter / body indices ...')
    prom_idx, body_idx, stats = load_gene_indices()
    print(f'  genes: {stats["total"]}, strand resolved: {stats["resolved"]}, '
          f'unresolved (body-only): {stats["unresolved"]}')

    peak_idx = load_peak_index(f'{PEAK_DIR}/HOXA9_peaks.bed')

    for stem, tf, tad_bed in RUNS:
        in_file = f'{RUN_DIR}/{stem}_10000_0.05/{MERGED_REL}'
        out_file = os.path.join(OUT_DIR, f'{stem}_merged.annotated.csv')
        if not os.path.exists(in_file):
            print(f'  {stem:20s} ({tf}): SKIP - merged file not found yet')
            continue
        bnd_idx = load_boundary_index(tad_bed)
        n = annotate(in_file, out_file, tf, peak_idx, bnd_idx, prom_idx, body_idx)
        print(f'  {stem:20s} ({tf}): {n:6d} interactions '
              f'[TAD={os.path.basename(tad_bed)}] -> {out_file}')


if __name__ == '__main__':
    main()
