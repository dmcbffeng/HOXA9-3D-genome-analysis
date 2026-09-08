"""Report empty temp_matrix files and missing TopDom .bed outputs."""
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
    ap.add_argument('--matrix-dir',  default=os.path.join(HERE, 'temp_matrix'))
    ap.add_argument('--results-dir', default=os.path.join(HERE, 'results'))
    ap.add_argument('--chroms', nargs='*', default=DEFAULT_CHROMS)
    args = ap.parse_args()

    samples = load_sample_names(args.samples_tsv)
    n_empty_mat = 0
    n_missing_bed = 0
    print('--- empty / missing temp_matrix files ---')
    for s in samples:
        for c in args.chroms:
            f = os.path.join(args.matrix_dir, f'{s}_{c}.matrix')
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                print(f'  {s}\t{c}\t({"missing" if not os.path.exists(f) else "empty"})')
                n_empty_mat += 1

    print('--- missing TopDom .bed outputs ---')
    for s in samples:
        for c in args.chroms:
            bed = os.path.join(args.results_dir, f'{s}_{c}.bed')
            if not os.path.exists(bed):
                print(f'  {s}\t{c}')
                n_missing_bed += 1

    print(f'\nempty/missing matrices : {n_empty_mat}')
    print(f'missing TopDom .bed    : {n_missing_bed}')


if __name__ == '__main__':
    main()
