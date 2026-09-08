import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def hic2matrix(hic_files, ch, ch_size, resolution=50000):
    all_contacts = np.zeros((ch_size, ch_size))

    overflow = False
    for hic_file in hic_files:
        print('======')
        print(f'Extract Hi-C file for {hic_file.split("/")[-1]}:{ch}')
        juicer = '/nfs/turbo/umms-drjieliu/juicer_tools_1.11.04_jcuda.0.8.jar'
        cmd = f'java -jar {juicer} dump oe NONE {hic_file} {ch} {ch} BP {resolution} temp.txt'
        os.system(cmd)
        print(f'Re-format the data into a matrix for {ch}')
        cnt_lines = 0
        for line in open('temp.txt'):
            if cnt_lines % 1000000 == 0:
                print(' Loading:', cnt_lines)
            cnt_lines += 1
            lst = line.strip().split()
            if lst[2].lower() == 'nan':
                continue
            p1, p2, v = int(lst[0]) // resolution, int(lst[1]) // resolution, float(lst[2])
            if p1 >= ch_size or p2 >= ch_size:
                overflow = True
                continue
            all_contacts[p1, p2] += v
            if p1 != p2:
                all_contacts[p2, p1] += v
        os.remove('temp.txt')
    if overflow:
        print('Warning: chr size might be wrong. Some values are out of bounds.')
    if np.sum(all_contacts) < 1e-5:
        print('Warning: no contacts found.')
    return all_contacts


def call_eigenvecs(
        hic_file, chromosomes, bg_prefix, hic_prefix,
        resolution=500_000
):
    # Load chromosome sizes
    chrom_sizes = {}
    for line in open('hg19.chrom.sizes'):
        lst = line.strip().split()
        ch, length = lst[0], int(lst[1])
        ll = length // resolution
        if length % resolution != 0:
            ll += 1
        chrom_sizes[ch] = ll

    table = []
    for ch in chromosomes:
        # Obtain OE mat
        mat = hic2matrix([hic_file], ch, chrom_sizes[ch], resolution=resolution)
        # Find non-zero columns
        orig_size = len(mat)
        col_sum = np.sum(mat, axis=1)
        non_zeros = [i for i, val in enumerate(col_sum) if val > 1e-5]
        mat = mat[non_zeros, :]
        mat = mat[:, non_zeros]

        # Get correlation matrix
        corr_mat = np.corrcoef(mat)
        # Get PC1 and PC2
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(corr_mat)
        _pc1, _pc2 = pcs[:, 0], pcs[:, 1]

        # insert the zeros back
        pc1, pc2 = np.zeros((orig_size,)), np.zeros((orig_size,))
        for idx, val in zip(non_zeros, _pc1):
            pc1[idx] = val
        for idx, val in zip(non_zeros, _pc2):
            pc2[idx] = val

        # Check H3K27ac
        bg_vec = np.load(f'00_bg_vectors/{bg_prefix}_{ch}.npy')
        assert len(bg_vec) == len(pc1)
        pc1_pos = np.where(pc1 > 1e-5, 1, 0)
        pc1_neg = np.where(pc1 < -1e-5, 1, 0)
        pc1_pos_marker = np.inner(pc1_pos, bg_vec)
        pc1_neg_marker = np.inner(pc1_neg, bg_vec)
        if pc1_pos_marker < pc1_neg_marker:
            pc1 = -pc1
            pc1_pos_marker, pc1_neg_marker = pc1_neg_marker, pc1_pos_marker
        pc2_pos = np.where(pc2 > 1e-5, 1, 0)
        pc2_neg = np.where(pc2 < -1e-5, 1, 0)
        pc2_pos_marker = np.inner(pc2_pos, bg_vec)
        pc2_neg_marker = np.inner(pc2_neg, bg_vec)
        if pc2_pos_marker < pc2_neg_marker:
            pc2 = -pc2
            pc2_pos_marker, pc2_neg_marker = pc2_neg_marker, pc2_pos_marker
        table.append(
            [ch, pc1_pos_marker, pc1_neg_marker, pc2_pos_marker, pc2_neg_marker]
        )
        pcs = np.array([pc1, pc2])
        np.save(f'01_AB_vec/{hic_prefix}_{ch}.npy', pcs)
    df = pd.DataFrame(table, columns=['chr', 'pc1_A', 'pc1_B', 'pc2_A', 'pc2_B'])
    df.to_csv(f'01_AB_vec/{hic_prefix}.csv', index=False)


if __name__ == '__main__':
    resolution = 500_000
    chromosomes = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]
    path = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/'

    call_eigenvecs(
        hic_file=f'{path}/GSM2861708_HIC1258_aligned_inter_30.hic',
        chromosomes=chromosomes,
        bg_prefix='O3_500kb',
        hic_prefix='HSPC',
        resolution=resolution
    )
    call_eigenvecs(
        hic_file=f'{path}/O3.hic',
        chromosomes=chromosomes,
        bg_prefix='O3_500kb',
        hic_prefix='O3',
        resolution=resolution
    )
    call_eigenvecs(
        hic_file=f'{path}/O3_HOXA9KO.allValidPairs.hic',
        chromosomes=chromosomes,
        bg_prefix='O3_500kb',
        hic_prefix='O3KO',
        resolution=resolution
    )
    call_eigenvecs(
        hic_file=f'{path}/4943.hic',
        chromosomes=chromosomes,
        bg_prefix='O3_500kb',
        hic_prefix='4943',
        resolution=resolution
    )

