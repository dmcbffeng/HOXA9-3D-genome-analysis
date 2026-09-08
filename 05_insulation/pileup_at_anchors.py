"""
O/E contact-map pile-ups centred on loop anchors -- the heatmaps under Fig. 3B / 4B.

For each named class of sites (one BED or peak-name file per class), dump the
observed/expected KR-normalized matrix at 10 kb (`juicer_tools dump oe KR`), take the
41 x 41 window (+/- 20 bins = +/- 200 kb) centred on the midpoint of every site, and
average over all sites genome-wide. Writes one `APA_<class>.npy` per class.

Fig. 3B classes: HOXA9 anchors of E-P, E-E and E-B loops, plus RAD21 and CTCF sites.
Fig. 4B classes: HOXA9 vs HOXA9/MEIS1 co-bound E-B anchors, plus MEIS1/PBX2/PBX3-only
controls -- same script, different `peak_name2file` entries.

Input files may be 4-column BED (`chr start end name`) or a single column of HOMER-style
peak names carrying the coordinates. Edit `hic_file` and `peak_name2file` at the bottom.
"""
import os
import subprocess
import numpy as np
import pandas as pd


# chromosomes to use
chroms = [f'chr{i}' for i in list(range(1, 23)) + ['X']]


def load_chrom_sizes(file='hg19.chrom.sizes'):
    chrom_sizes = {}
    for line in open(file):
        lst = line.strip().split()
        if lst[0] in chroms:
            chrom_sizes[lst[0]] = int(lst[1])
    return chrom_sizes


def parse_boundary_location(peak_file, resolution=10000):
    df_peaks = pd.read_csv(peak_file, sep='\s+', header=None)
    # If it is a bed file with 4 columns
    if len(df_peaks.columns) >= 4:
        first_4_cols = df_peaks.columns[:4]
        col_names = ['chr', 'start', 'end', 'name']
        mapping = dict(zip(first_4_cols, col_names))
        df_peaks = df_peaks.rename(columns=mapping)
        df_peaks['loc'] = (df_peaks['start'] + df_peaks['end']) / 2  # use the mid point of the peak
    # If it is just peak name, with coordinates in it
    elif len(df_peaks.columns) == 1:
        df_peaks.columns = ['name']
        df_peaks[['UnUsed1', 'chr', 'loc', 'UnUsed2']] = df_peaks['name'].str.split('-', expand=True)
        df_peaks['loc'] = df_peaks['loc'].astype(int)
    peak_dict = {}
    for ch in chroms:
        subdf = df_peaks[df_peaks['chr'] == ch]
        locations = [int(elm // resolution) for elm in subdf['loc']]
        peak_dict[ch] = locations
    return peak_dict


def APA_peaks(hic_file, peak_name2file, resolution=10000, window=20):
    # Load peaks from all files
    name2peaks = {
        name: parse_boundary_location(
            peak_name2file[name], resolution=resolution
        ) for name in peak_name2file
    }
    # load chrom sizes
    chrom_sizes = load_chrom_sizes()

    apa_number_regions_per_chrom = {
        name: {ch: 0 for ch in chroms} for name in name2peaks
    }
    apa_sum_contact_map_per_chrom = {
        name: {ch: np.zeros((2 * window + 1, 2 * window + 1)) for ch in chroms} for name in name2peaks
    }

    for ch in chroms:
        print('APA on', hic_file, ch)

        # Load numpy matrix from hic file
        juicer = './juicer_tools_1.11.04_jcuda.0.8.jar'
        # OE and KR normalization
        cmd = f'java -jar {juicer} dump oe KR {hic_file} {ch} {ch} BP {resolution} temp.txt'
        print(cmd)
        cmd = [
            "java", "-jar", juicer,
            "dump", "oe", "KR", hic_file,
            ch, ch, "BP", "10000", "temp.txt"
        ]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"Error running {ch}:")
            print(result.stderr)
        else:
            print(f"{ch} done.")

        if not os.path.exists(f'temp.txt'):
            print(f'KR norm failed for {ch}, skipped this chromosome')
            continue

        max_size = chrom_sizes[ch] // resolution + 1
        mat = np.zeros((max_size, max_size))
        for line in open(f'temp.txt'):
            lst = line.strip().split()
            if lst[2].lower() == 'nan':
                continue
            p1, p2, v = int(lst[0]) // resolution, int(lst[1]) // resolution, float(lst[2])
            mat[p1, p2] += v
            if p1 != p2:
                mat[p2, p1] += v
        os.remove(f'temp.txt')

        if np.sum(mat) == 0:
            print(f'KR norm failed for {ch}, skipped this chromosome')
            continue

        # Run apa for each group of peaks
        for name in name2peaks:
            peaks_this_chr = name2peaks[name][ch]
            for loc in peaks_this_chr:
                if loc < window or loc >= max_size - window:  # too close to start and end of the chromosome
                    continue
                submat = mat[loc - window:loc + window + 1, loc - window:loc + window + 1]
                apa_number_regions_per_chrom[name][ch] += 1
                apa_sum_contact_map_per_chrom[name][ch] += submat

        # Take average and save the results
        for name in name2peaks:
            sum_mat, number_mats = np.zeros((2 * window + 1, 2 * window + 1)), 0
            for ch in chroms:
                sum_mat += apa_sum_contact_map_per_chrom[name][ch]
                number_mats += apa_number_regions_per_chrom[name][ch]
            average_mat = sum_mat / number_mats
            # Save the matrices
            np.save(f'APA_{name}.npy', average_mat)



if __name__ == '__main__':
    hic_file = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/O3.hic'
    peak_name2file = {
        'RAD21': 'OCI-AML3_RAD21_peaks.txt',
        'CTCF': 'OCI-AML3_CTCF_peaks.txt',
        'HOXA9_EP': 'OCI-AML3_HOXA9_EP_peaks.txt',
        'HOXA9_EE': 'OCI-AML3_HOXA9_EE_peaks.txt',
        'HOXA9_EB': 'OCI-AML3_HOXA9_EB_peaks.txt',
    }

    APA_peaks(
        hic_file, peak_name2file,
        resolution=10000, window=20
    )

