import numpy as np
import pandas as pd
from hicrep.utils import readMcool
from hicrep import hicrepSCC


# names = [
#     '4943', '5577', '5638', '5832_CD34neg', '5832_CD34pos',
#     '6281_AZA', '6281_DMSO', '6527', 'CD34_HSPC_GLCdel', 'GSM4604271_576_Hi-C',
#     'GSM4604272_103_Hi-C', 'GSM4604273_673_Hi-C', 'GSM4604274_629_Hi-C', 'GSM4604275_270_Hi-C', 'GSM4604276_868_Hi-C',
#     'GSM4604277_773_Hi-C', 'GSM4604278_1953_Hi-C', 'GSM4604279_566_Hi-C', 'GSM4604280_018_Hi-C', 'GSM4604281_260_Hi-C',
#     'GSM4604282_546_Hi-C', 'GSM4604283_1021_Hi-C', 'GSM4604284_1916_Hi-C', 'GSM4604285_472_Hi-C', 'GSM4604286_936_Hi-C',
#     'GSM4604287_1360_Hi-C', 'GSM4604288_168_Hi-C', 'GSM4604289_027_Hi-C', 'GSM4604290_990_Hi-C', 'GSM5688515_1071',
#     'GSM5688516_1324', 'GSM5688517_424', 'GSM5688518_496', 'GSM5688519_504', 'GSM5688520_798',
#     'GSM5688521_838', 'GSM5688522_933', 'GSM5688523_HSPC-050', 'GSM5688524_HSPC-213', 'GSM5688525_HSPC-341',
#     'GSM5688526_HSPC-822', 'HSPC_downsample_allValidPairs', 'Kasumi', 'MLL1KO', 'MV411',
#     'O3']
names = [
    'GSM4604271_576_Hi-C',
    'GSM4604272_103_Hi-C', 'GSM4604273_673_Hi-C', 'GSM4604274_629_Hi-C', 'GSM4604275_270_Hi-C', 'GSM4604276_868_Hi-C',
    'GSM4604277_773_Hi-C', 'GSM4604278_1953_Hi-C', 'GSM4604279_566_Hi-C', 'GSM4604280_018_Hi-C', 'GSM4604281_260_Hi-C',
    'GSM4604282_546_Hi-C', 'GSM4604283_1021_Hi-C', 'GSM4604284_1916_Hi-C', 'GSM4604285_472_Hi-C', 'GSM4604286_936_Hi-C',
    'GSM4604287_1360_Hi-C', 'GSM4604288_168_Hi-C', 'GSM4604289_027_Hi-C', 'GSM4604290_990_Hi-C', 'GSM5688515_1071',
    'GSM5688516_1324', 'GSM5688517_424', 'GSM5688518_496', 'GSM5688519_504', 'GSM5688520_798',
    'GSM5688521_838', 'GSM5688522_933', 'GSM5688523_HSPC-050', 'GSM5688524_HSPC-213', 'GSM5688525_HSPC-341',
    'GSM5688526_HSPC-822']
chroms = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
temp_path = '/scratch/drjieliu_root/drjieliu/fanfeng/Xiaotian_data/AML'
resolution = 50000


for ch in chroms:
    print(ch)

    all_cools = []
    for name in names:
        cool_file = f'{temp_path}/{name}_{ch}_50000.cool'
        print('Loading', cool_file.split('/')[-1])
        cool, binSize = readMcool(cool_file, -1)
        all_cools.append(cool)

    mat = np.identity(len(names))
    for i in range(len(all_cools) - 1):
        for j in range(i + 1, len(all_cools)):
            try:
                scc = hicrepSCC(all_cools[i], all_cools[j], h=1, dBPMax=40 * resolution, bDownSample=False,
                                chrNames=[ch[3:]])
            except ValueError:
                scc = hicrepSCC(all_cools[i], all_cools[j], h=1, dBPMax=40 * resolution, bDownSample=False,
                                chrNames=[ch])
            print(scc)
            mat[i, j] = scc[0]
            mat[j, i] = scc[0]

    df = pd.DataFrame(mat, columns=names)
    df.to_csv(f'{ch}.csv', index=False)



