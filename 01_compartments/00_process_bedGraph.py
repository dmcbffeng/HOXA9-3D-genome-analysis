import numpy as np


# Convert BedGraph files into vectors at given resolutions
resolution = 500_000
chromosomes = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]

input_file = 'O3_H3K27ac.bedGraph'
output_file = 'O3_500kb'
# input_file = 'HSPC_H3K27ac_hg19.bedGraph'
# output_file = 'HSPC_500kb'

# Load chromosome sizes
chrom_sizes = {}
for line in open('hg19.chrom.sizes'):
    lst = line.strip().split()
    ch, length = lst[0], int(lst[1])
    ll = length // resolution
    if length % resolution != 0:
        ll += 1
    chrom_sizes[ch] = ll

# Convert each chromosome into one vector and store as
bg_vectors = {
    ch: np.zeros((chrom_sizes[ch],)) for ch in chromosomes
}
for line in open(input_file):
    lst = line.strip().split()
    ch, start, end, val = lst[0], int(lst[1]), int(lst[2]), float(lst[3])
    if ch not in bg_vectors:
        continue
    st, ed = start // resolution, end // resolution
    for bin_id in range(st, ed + 1):
        real_start = max(bin_id * resolution, start)
        real_end = min((bin_id + 1) * resolution, end)
        normalized_val = val * (real_end - real_start) / resolution
        bg_vectors[ch][bin_id] += normalized_val
for ch in chromosomes:
    np.save(
        f'00_bg_vectors/{output_file}_{ch}.npy', bg_vectors[ch]
    )
