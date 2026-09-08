"""Concatenate per-chromosome TopDom .bed files into one bed per sample."""
import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CHROMS = [f'chr{c}' for c in list(range(1, 23)) + ['X']]


def load_sample_names(tsv_path):
    names = []
    with open(tsv_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            names.append(ln.split('\t')[0])
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples-tsv', default=os.path.join(HERE, 'samples.tsv'))
    ap.add_argument('--results-dir', default=os.path.join(HERE, 'results'))
    ap.add_argument('--merged-dir',  default=os.path.join(HERE, 'merged'))
    ap.add_argument('--chroms', nargs='*', default=DEFAULT_CHROMS)
    args = ap.parse_args()

    os.makedirs(args.merged_dir, exist_ok=True)
    samples = load_sample_names(args.samples_tsv)

    for s in samples:
        out = os.path.join(args.merged_dir, f'{s}.bed')
        n = 0
        with open(out, 'w') as fout:
            for c in args.chroms:
                p = os.path.join(args.results_dir, f'{s}_{c}.bed')
                if os.path.exists(p):
                    with open(p) as fin:
                        fout.write(fin.read())
                    n += 1
        print(f'{s}: merged {n}/{len(args.chroms)} chroms -> {out}')


if __name__ == '__main__':
    main()
