"""
`.hic` -> per-chromosome `.cool` at a fixed resolution, via HiCExplorer's
`hicConvertFormat`. This is the input step for `hicrep_scc.py`, which reads
`.cool` files through `hicrep.utils.readMcool`.

Edit `hic_path`, `temp_path`, `chr_size_file` and `resolution` before running.
"""
import os


chroms = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
hic_path = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic'
temp_path = '/scratch/drjieliu_root/drjieliu/fanfeng/Xiaotian_data/AML'
# output_path = '/nfs/turbo/umms-drjieliu/usr/temp_Fan/AML_Analysis'
chr_size_file = '/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/Yue_AML/pairs/hg38.chrom.sizes'
resolution = 50000


for ch in chroms:
    print(ch)
    cool_files = []
    # Step 1: convert .hic into .cool
    for path, subdir, files in os.walk(hic_path):
        for file in files:
            if not file.endswith('.hic'):
                continue
            print('Converting', file)
            f_input = os.path.join(path, file)
            prefix = file.split('.')[0]
            f_output = os.path.join(temp_path, prefix + f'_{ch}.cool')
            cool_files.append(f_output)

            cmd = f'hicConvertFormat --matrices {f_input} --outFileName {f_output} --inputFormat hic --outputFormat cool'
            cmd += f' --chromosome {ch} --resolutions {resolution} --chromosomeSizes {chr_size_file}'

            os.system(cmd)

