import numpy as np


chroms = [f'chr{i}' for i in list(range(1, 23)) + ['X']]


def load_chrom_regions(resolution=10000):
    all_regions = []
    ch2sizes = {}
    for line in open('hg19.chrom.sizes'):
        lst = line.strip().split()
        if lst:
            ch, length = lst[0], int(lst[1]) // resolution
            if ch in chroms:
                ch2sizes[ch] = length
                chr_regions = [(ch, i) for i in range(length)]
                all_regions.extend(chr_regions)
    return all_regions, ch2sizes


def load_HOXA_sites(file_path='O3_HOXA9_HA_overlap_peaks_IgGcorr.bed', resolution=10000):
    HOXA_sites = {ch: {} for ch in chroms}
    for line in open(file_path):
        lst = line.strip().split()
        ch, p1, p2, name = lst[0], int(lst[1]), int(lst[2]), lst[3]
        if ch not in chroms:
            continue
        st, ed = p1 // resolution, p2 // resolution
        for i in range(st, ed + 1):
            if i not in HOXA_sites[ch]:
                HOXA_sites[ch][i] = []
            HOXA_sites[ch][i].append((ch, p1, p2, name))
    return HOXA_sites


def sample_region(input_file, output_prefix, resolution=10000):
    np.random.seed(0)
    all_regions, ch2sizes = load_chrom_regions(resolution=resolution)
    all_regions_idx = np.arange(len(all_regions))
    HOXA_sites = load_HOXA_sites(resolution=resolution)

    # Load all distances and HOXA sites
    hoxa_sites = []
    distances = []
    for line in open(input_file):
        lst = line.strip().split()
        if lst[1] == 's1':
            continue
        ch, p1, p2 = lst[0], int(lst[1]) // resolution, int(lst[4]) // resolution
        if ch not in chroms:
            continue
        distances.append(abs(p2 - p1))
        anchor1, anchor2 = p1 in HOXA_sites[ch], p2 in HOXA_sites[ch]
        assert anchor1 or anchor2
        if anchor1 and anchor2:  # random sample one of them
            if np.random.random() > 0.5:
                anchor1 = False
            else:
                anchor2 = False
        if anchor1:
            hoxa_sites.append((ch, p1))
        else:
            hoxa_sites.append((ch, p2))

    np.random.shuffle(distances)

    # Sample method 1: total random
    random_contact = []
    for i in range(len(distances)):
        while True:
            idx = np.random.choice(all_regions_idx)
            ch, p1 = all_regions[idx]
            p2 = p1 + distances[i]
            if p2 < ch2sizes[ch]:
                break
        random_contact.append((ch, p1, p2))
    random_contact.sort()
    with open(f'{output_prefix}_random1.txt', 'w') as f:
        for (ch, p1, p2) in random_contact:
            f.write(f'{ch}\t{p1 * resolution}\t{p2 * resolution}\n')

    # Sample method 2: keep HOXA9 anchor
    random_contact = []
    for (ch, p0), dis in zip(hoxa_sites, distances):
        if p0 < dis:
            random_contact.append((ch, p0, p0 + dis))
        elif p0 >= ch2sizes[ch] - dis:
            random_contact.append((ch, p0 - dis, p0))
        elif np.random.random() > 0.5:
            random_contact.append((ch, p0, p0 + dis))
        else:
            random_contact.append((ch, p0 - dis, p0))
    random_contact.sort()
    with open(f'{output_prefix}_random2.txt', 'w') as f:
        for (ch, p1, p2) in random_contact:
            f.write(f'{ch}\t{p1 * resolution}\t{p2 * resolution}\n')


if __name__ == '__main__':
    sample_region(
        input_file='../O3_HOXA9_10000_0.05/FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed',
        output_prefix='O3',
        resolution=10000
    )

    sample_region(
        input_file='../O3KO_HOXA9_10000_0.05/FitHiChIP.interactions_FitHiC_Q0.05_MergeNearContacts.bed',
        output_prefix='O3KO',
        resolution=10000
    )


