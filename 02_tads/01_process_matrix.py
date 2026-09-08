"""
Dump per-chromosome contact matrices from .hic files for TopDom TAD calling.

For each sample listed in samples.tsv we dump `oe KR` (observed/expected,
KR-normalized) values with juicer_tools and write an (N x (N+3))
TopDom-format matrix to temp_matrix/<sample>_<chr>.matrix.

We do NOT divide by any "common factor" here: a few of the input .hic files
have raw observed counts that are all multiples of 30 or 60 (see claude.md),
but because we extract `oe`, the per-stratum expected divides out the same
constant factor and the OE matrix is unaffected. If you ever switch to
`observed KR` or `observed NONE`, revisit this.

Usage:
    python 01_process_matrix.py               # all samples, all autosomes + chrX
    python 01_process_matrix.py --samples 4943 MOLM13
    python 01_process_matrix.py --resolution 25000
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd


HERE = os.path.dirname(os.path.abspath(__file__))
JUICER = '/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar'
DEFAULT_CHROMS = [f'chr{c}' for c in list(range(1, 23)) + ['X']]


# ---------- sample list -------------------------------------------------------

def load_samples(tsv_path):
    samples = []
    with open(tsv_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = ln.split('\t')
            if len(parts) < 2:
                print(f'[WARN] skipping malformed line: {ln!r}', file=sys.stderr)
                continue
            samples.append((parts[0], parts[1]))
    return samples


# ---------- matrix dump -------------------------------------------------------

def hic2matrix(hic_file, chrom, sample_name, resolution, out_dir):
    max_size = 250_000_000 // resolution
    contacts = np.zeros((max_size, max_size), dtype=np.float64)
    actual_max = 0

    tmp = os.path.join(out_dir, f'_tmp_{sample_name}_{chrom}.txt')
    cmd = ['java', '-jar', JUICER, 'dump', 'oe', 'KR',
           hic_file, chrom, chrom, 'BP', str(resolution), tmp]
    print(f'[dump]  {sample_name} {chrom}')
    res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(f'[dump]  FAILED for {sample_name} {chrom}', file=sys.stderr)
        if res.stderr:
            print(res.stderr.decode(errors='ignore')[:400], file=sys.stderr)
        if os.path.exists(tmp):
            os.remove(tmp)
        return

    with open(tmp) as fh:
        for i, line in enumerate(fh):
            if i and i % 1_000_000 == 0:
                print(f'  loading line {i:,}')
            lst = line.strip().split()
            if len(lst) < 3 or lst[2].lower() == 'nan':
                continue
            p1 = int(lst[0]) // resolution
            p2 = int(lst[1]) // resolution
            v = float(lst[2])
            contacts[p1, p2] += v
            if p1 != p2:
                contacts[p2, p1] += v
            if p1 + 1 > actual_max:
                actual_max = p1 + 1
            if p2 + 1 > actual_max:
                actual_max = p2 + 1
    os.remove(tmp)

    if actual_max == 0:
        print(f'[dump]  WARN empty matrix for {sample_name} {chrom}', file=sys.stderr)
        # write an empty file so downstream check picks it up
        open(os.path.join(out_dir, f'{sample_name}_{chrom}.matrix'), 'w').close()
        return

    contacts = contacts[:actual_max, :actual_max]
    first = pd.DataFrame({
        'chr':        [chrom] * actual_max,
        'from.coord': [i * resolution for i in range(actual_max)],
        'to.coord':   [(i + 1) * resolution for i in range(actual_max)],
    })
    body = pd.DataFrame(contacts, columns=[str(i) for i in range(actual_max)])
    df = pd.concat([first, body], axis=1)
    out_path = os.path.join(out_dir, f'{sample_name}_{chrom}.matrix')
    df.to_csv(out_path, index=False, header=False, sep='\t')


# ---------- driver ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples-tsv', default=os.path.join(HERE, 'samples.tsv'))
    ap.add_argument('--out-dir', default=os.path.join(HERE, 'temp_matrix'))
    ap.add_argument('--log-dir', default=os.path.join(HERE, 'logs'))
    ap.add_argument('--resolution', type=int, default=50000)
    ap.add_argument('--samples', nargs='*', default=None,
                    help='If given, only run these sample names.')
    ap.add_argument('--chroms', nargs='*', default=DEFAULT_CHROMS)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)

    all_samples = load_samples(args.samples_tsv)
    if args.samples:
        wanted = set(args.samples)
        all_samples = [s for s in all_samples if s[0] in wanted]
        missing = wanted - {s[0] for s in all_samples}
        if missing:
            print(f'[WARN] requested samples not in {args.samples_tsv}: {missing}',
                  file=sys.stderr)

    if not all_samples:
        print('No samples to process.', file=sys.stderr)
        sys.exit(1)

    print(f'Will process {len(all_samples)} sample(s) at {args.resolution} bp:')
    for name, path in all_samples:
        print(f'  {name}\t{path}')

    for sample_name, hic_file in all_samples:
        if not os.path.exists(hic_file):
            print(f'[ERROR] {hic_file} not found — skipping {sample_name}',
                  file=sys.stderr)
            continue
        for chrom in args.chroms:
            out_file = os.path.join(args.out_dir, f'{sample_name}_{chrom}.matrix')
            if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
                print(f'[skip]  {out_file} already exists')
                continue
            hic2matrix(hic_file, chrom, sample_name, args.resolution,
                       args.out_dir)


if __name__ == '__main__':
    main()
