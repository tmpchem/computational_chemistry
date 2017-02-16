
"""Functions and tables for AMBER94 molecular mechanics parameters."""

import sys

# relative atomic masses of elements (in atomic mass units [g/mol]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {  'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067,
'F' : 18.9984, 'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
'Be': 9.01218, 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
'Si': 28.0855, 'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X' : 0.00000}

# covalent (or ionic) radii by atomic element [Angstroms] from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
cov_radii = {'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
 'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
 'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
 'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
 'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
 'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
 'Se': 1.17, 'Kr': 1.03, 'X' : 0.00}

# AMBER94 molecular mechanics van der waals parameters for atom types:
# 1 -> (float) ro/2 [Angstrom], van der waals radius (divided by 2)
# 2 -> (float) eps [kcal/mol], van der waals attraction magnitude
vdw_params = {          'C' : [1.9080, 0.0860], 'CA': [1.9080, 0.0860],
'CM': [1.9080, 0.0860], 'CC': [1.9080, 0.0860], 'CV': [1.9080, 0.0860],
'CW': [1.9080, 0.0860], 'CR': [1.9080, 0.0860], 'CB': [1.9080, 0.0860],
'C*': [1.9080, 0.0860], 'CN': [1.9080, 0.0860], 'CK': [1.9080, 0.0860],
'CQ': [1.9080, 0.0860], 'CT': [1.9080, 0.1094], 'CS': [3.3950, 0.0000806],
'F' : [1.7500, 0.0610], 'H' : [0.6000, 0.0157], 'H1': [1.3870, 0.0157],
'H2': [1.2870, 0.0157], 'H3': [1.1870, 0.0157], 'H4': [1.4090, 0.0150],
'H5': [1.3590, 0.0150], 'HA': [1.4590, 0.0150], 'HC': [1.4870, 0.0157],
'HO': [0.0001, 0.0000], 'HP': [1.1000, 0.0157], 'HS': [0.6000, 0.0157],
'HW': [0.0001, 0.0000], 'Li': [1.1370, 0.0183], 'IP': [1.8680, 0.00277],
'Li': [1.1370, 0.0183], 'N' : [1.8240, 0.1700], 'K' : [2.6580, 0.000328],
'O' : [1.6612, 0.2100], 'O2': [1.6612, 0.2100], 'OH': [1.7210, 0.2104],
'OS': [1.6837, 0.1700], 'OW': [1.7683, 0.1520], 'P' : [2.1000, 0.2000],
'S' : [2.0000, 0.2500], 'SH': [2.0000, 0.2500], 'N' : [1.8240, 0.1700],
'NA': [1.8240, 0.1700], 'NB': [1.8240, 0.1700], 'NC': [1.8240, 0.1700],
'N*': [1.8240, 0.1700], 'N2': [1.8240, 0.1700], 'N3': [1.8750, 0.1700],
'He': [1.5800, 0.0112], 'Ar': [1.8436, 0.4466], 'HH': [1.5000, 0.0000]}

# AMBER94 molecular mechanics bond parameters for atom type pairs:
# 1 -> (float) k_b [kcal/(mol*A^2)], bond spring constant
# 2 -> (float) r_eq [Angstrom], equilibrium bond length
bond_params = {         'C CA': [469.0, 1.409], 'C CB': [447.0, 1.419],
'C CM': [410.0, 1.444], 'C CT': [317.0, 1.522], 'C N ': [490.0, 1.335],
'C N*': [424.0, 1.383], 'C NA': [418.0, 1.388], 'C NC': [457.0, 1.358],
'C O ': [570.0, 1.229], 'C O2': [656.0, 1.250], 'C OH': [450.0, 1.364],
'C*CB': [388.0, 1.459], 'C*CT': [317.0, 1.495], 'C*CW': [546.0, 1.352],
'C*HC': [367.0, 1.080], 'CACA': [469.0, 1.400], 'CACB': [469.0, 1.404],
'CACM': [427.0, 1.433], 'CACN': [469.0, 1.400], 'CACT': [317.0, 1.510],
'CAH4': [367.0, 1.080], 'CAHA': [367.0, 1.080], 'CAN2': [481.0, 1.340],
'CANA': [427.0, 1.381], 'CANC': [483.0, 1.339], 'CBCB': [520.0, 1.370],
'CBCN': [447.0, 1.419], 'CBN*': [436.0, 1.374], 'CBNB': [414.0, 1.391],
'CBNC': [461.0, 1.391], 'CCCT': [317.0, 1.504], 'CCCV': [512.0, 1.375],
'CCCW': [518.0, 1.371], 'CCNA': [422.0, 1.385], 'CCNB': [410.0, 1.394],
'CKH5': [367.0, 1.080], 'CKN*': [440.0, 1.371], 'CKNB': [529.0, 1.304],
'CMCM': [549.0, 1.350], 'CMCT': [317.0, 1.510], 'CMH4': [367.0, 1.080],
'CMH5': [367.0, 1.080], 'CMHA': [367.0, 1.080], 'CMN*': [448.0, 1.365],
'CNNA': [428.0, 1.380], 'CQH5': [367.0, 1.080], 'CQNC': [502.0, 1.324],
'CRH5': [367.0, 1.080], 'CRNA': [477.0, 1.343], 'CRNB': [488.0, 1.335],
'CTCT': [310.0, 1.526], 'CTF ': [367.0, 1.380], 'CTH1': [340.0, 1.090],
'CTH2': [340.0, 1.090], 'CTH3': [340.0, 1.090], 'CTHC': [340.0, 1.090],
'CTHP': [340.0, 1.090], 'CTN ': [337.0, 1.449], 'CTN*': [337.0, 1.475],
'CTN2': [337.0, 1.350], 'CTN3': [367.0, 1.471], 'CTOH': [320.0, 1.410],
'CTOS': [320.0, 1.410], 'CTS ': [227.0, 1.810], 'CTSH': [237.0, 1.810],
'CVH4': [367.0, 1.080], 'CVNB': [410.0, 1.394], 'CWH4': [367.0, 1.080],
'CWNA': [410.0, 1.394], 'H N ': [434.0, 1.010], 'H N*': [434.0, 1.010],
'H N2': [434.0, 1.010], 'H N3': [434.0, 1.010], 'H NA': [434.0, 1.010],
'HOOH': [553.0, 0.960], 'HOOS': [553.0, 0.960], 'HSSH': [274.0, 1.336],
'O2P ': [525.0, 1.480], 'OHP ': [230.0, 1.610], 'OSP ': [230.0, 1.610],
'OWHW': [553.0,0.9572], 'S S ': [166.0, 2.038], 'HHHH': [100.0, 0.740]}

# AMBER94 molecular mechanics angle parameters for atom type triplets:
# 1 -> (float) k_a [kcal/(mol*rad^2)], angle spring constant
# 2 -> (float) a_eq [degrees], equilibrium bond angle
angle_params = {           'C CACA': [ 63.0, 120.00],
'C CAHA': [ 35.0, 120.00], 'C CBCB': [ 63.0, 119.20],
'C CBNB': [ 70.0, 130.00], 'C CMCM': [ 63.0, 120.70],
'C CMCT': [ 70.0, 119.70], 'C CMH4': [ 35.0, 119.70],
'C CMHA': [ 35.0, 119.70], 'C CTCT': [ 63.0, 111.10],
'C CTH1': [ 50.0, 109.50], 'C CTHC': [ 50.0, 109.50],
'C CTHP': [ 50.0, 109.50], 'C CTN ': [ 63.0, 110.10],
'C CTN3': [ 80.0, 111.20], 'C N CT': [ 50.0, 121.90],
'C N H ': [ 30.0, 120.00], 'C N*CM': [ 70.0, 121.60],
'C N*CT': [ 70.0, 117.60], 'C N*H ': [ 30.0, 119.20],
'C NAC ': [ 70.0, 126.40], 'C NACA': [ 70.0, 125.20],
'C NAH ': [ 30.0, 116.80], 'C NCCA': [ 70.0, 120.50],
'C OHHO': [ 35.0, 113.00], 'C*CBCA': [ 63.0, 134.90],
'C*CBCN': [ 63.0, 108.80], 'C*CTCT': [ 63.0, 115.60],
'C*CTHC': [ 50.0, 109.50], 'C*CWH4': [ 35.0, 120.00],
'C*CWNA': [ 70.0, 108.70], 'CAC CA': [ 63.0, 120.00],
'CAC OH': [ 70.0, 120.00], 'CACACA': [ 63.0, 120.00],
'CACACB': [ 63.0, 120.00], 'CACACN': [ 63.0, 120.00],
'CACACT': [ 70.0, 120.00], 'CACAH4': [ 35.0, 120.00],
'CACAHA': [ 35.0, 120.00], 'CACBCB': [ 63.0, 117.30],
'CACBCN': [ 63.0, 116.20], 'CACBNB': [ 70.0, 132.40],
'CACMCM': [ 63.0, 117.00], 'CACMH4': [ 35.0, 123.30],
'CACMHA': [ 35.0, 123.30], 'CACNCB': [ 63.0, 122.70],
'CACNNA': [ 70.0, 132.80], 'CACTCT': [ 63.0, 114.00],
'CACTHC': [ 50.0, 109.50], 'CAN2CT': [ 50.0, 123.20],
'CAN2H ': [ 35.0, 120.00], 'CANAH ': [ 30.0, 118.00],
'CANCCB': [ 70.0, 112.20], 'CANCCQ': [ 70.0, 118.60],
'CBC NA': [ 70.0, 111.30], 'CBC O ': [ 80.0, 128.80],
'CBC*CT': [ 70.0, 128.60], 'CBC*CW': [ 63.0, 106.40],
'CBCAH4': [ 35.0, 120.00], 'CBCAHA': [ 35.0, 120.00],
'CBCAN2': [ 70.0, 123.50], 'CBCANC': [ 70.0, 117.30],
'CBCBN*': [ 70.0, 106.20], 'CBCBNB': [ 70.0, 110.40],
'CBCBNC': [ 70.0, 127.70], 'CBCNNA': [ 70.0, 104.40],
'CBN*CK': [ 70.0, 105.40], 'CBN*CT': [ 70.0, 125.80],
'CBN*H ': [ 30.0, 125.80], 'CBNBCK': [ 70.0, 103.80],
'CBNCCQ': [ 70.0, 111.00], 'CCCTCT': [ 63.0, 113.10],
'CCCTHC': [ 50.0, 109.50], 'CCCVH4': [ 35.0, 120.00],
'CCCVNB': [ 70.0, 120.00], 'CCCWH4': [ 35.0, 120.00],
'CCCWNA': [ 70.0, 120.00], 'CCNACR': [ 70.0, 120.00],
'CCNAH ': [ 30.0, 120.00], 'CCNBCR': [ 70.0, 117.00],
'CKN*CT': [ 70.0, 128.80], 'CKN*H ': [ 30.0, 128.80],
'CMC NA': [ 70.0, 114.10], 'CMC O ': [ 80.0, 125.30],
'CMCAN2': [ 70.0, 120.10], 'CMCANC': [ 70.0, 121.50],
'CMCMCT': [ 70.0, 119.70], 'CMCMH4': [ 35.0, 119.70],
'CMCMHA': [ 35.0, 119.70], 'CMCMN*': [ 70.0, 121.20],
'CMCTHC': [ 50.0, 109.50], 'CMN*CT': [ 70.0, 121.20],
'CMN*H ': [ 30.0, 121.20], 'CNCAHA': [ 35.0, 120.00],
'CNNACW': [ 70.0, 111.60], 'CNNAH ': [ 30.0, 123.10],
'CRNACW': [ 70.0, 120.00], 'CRNAH ': [ 30.0, 120.00],
'CRNBCV': [ 70.0, 117.00], 'CTC N ': [ 70.0, 116.60],
'CTC O ': [ 80.0, 120.40], 'CTC O2': [ 70.0, 117.00],
'CTC*CW': [ 70.0, 125.00], 'CTCCCV': [ 70.0, 120.00],
'CTCCCW': [ 70.0, 120.00], 'CTCCNA': [ 70.0, 120.00],
'CTCCNB': [ 70.0, 120.00], 'CTCTCT': [ 40.0, 109.50],
'CTCTH1': [ 50.0, 109.50], 'CTCTH2': [ 50.0, 109.50],
'CTCTHC': [ 50.0, 109.50], 'CTCTHP': [ 50.0, 109.50],
'CTCTN ': [ 80.0, 109.70], 'CTCTN*': [ 50.0, 109.50],
'CTCTN2': [ 80.0, 111.20], 'CTCTN3': [ 80.0, 111.20],
'CTCTOH': [ 50.0, 109.50], 'CTCTOS': [ 50.0, 109.50],
'CTCTS ': [ 50.0, 114.70], 'CTCTSH': [ 50.0, 108.60],
'CTN CT': [ 50.0, 118.00], 'CTN H ': [ 30.0, 118.04],
'CTN2H ': [ 35.0, 118.40], 'CTN3H ': [ 50.0, 109.50],
'CTOHHO': [ 55.0, 108.50], 'CTOSCT': [ 60.0, 109.50],
'CTOSP ': [100.0, 120.50], 'CTS CT': [ 62.0,  98.90],
'CTS S ': [ 68.0, 103.70], 'CTSHSH': [ 43.0,  96.00],
'CVCCNA': [ 70.0, 120.00], 'CWCCNA': [ 70.0, 120.00],
'CWCCNB': [ 70.0, 120.00], 'CWNAH ': [ 30.0, 120.00],
'F CTF ': [ 77.0, 109.10], 'F CTH1': [ 35.0, 109.50],
'H N H ': [ 35.0, 120.00], 'H N2H ': [ 35.0, 120.00],
'H N3H ': [ 35.0, 109.50], 'H1CTH1': [ 35.0, 109.50],
'H1CTN ': [ 50.0, 109.50], 'H1CTN*': [ 50.0, 109.50],
'H1CTN2': [ 50.0, 109.50], 'H1CTOH': [ 50.0, 109.50],
'H1CTOS': [ 50.0, 109.50], 'H1CTS ': [ 50.0, 109.50],
'H1CTSH': [ 50.0, 109.50], 'H2CTH2': [ 35.0, 109.50],
'H2CTN*': [ 50.0, 109.50], 'H2CTOS': [ 50.0, 109.50],
'H4CMN*': [ 35.0, 119.10], 'H4CVNB': [ 35.0, 120.00],
'H4CWNA': [ 35.0, 120.00], 'H5CKN*': [ 35.0, 123.05],
'H5CKNB': [ 35.0, 123.05], 'H5CQNC': [ 35.0, 115.45],
'H5CRNA': [ 35.0, 120.00], 'H5CRNB': [ 35.0, 120.00],
'HCCTHC': [ 35.0, 109.50], 'HOOHP ': [ 45.0, 108.50],
'HPCTHP': [ 35.0, 109.50], 'HPCTN3': [ 50.0, 109.50],
'HSSHHS': [ 35.0,  92.07], 'HWOWHW': [100.0, 104.52],
'N C O ': [ 80.0, 122.90], 'N*C NA': [ 70.0, 115.40],
'N*C NC': [ 70.0, 118.60], 'N*C O ': [ 80.0, 120.90],
'N*CBNC': [ 70.0, 126.20], 'N*CKNB': [ 70.0, 113.90],
'N*CTOS': [ 50.0, 109.50], 'N2CAN2': [ 70.0, 120.00],
'N2CANA': [ 70.0, 116.00], 'N2CANC': [ 70.0, 119.30],
'NAC O ': [ 80.0, 120.60], 'NACANC': [ 70.0, 123.30],
'NACRNA': [ 70.0, 120.00], 'NACRNB': [ 70.0, 120.00],
'NCC O ': [ 80.0, 122.50], 'NCCQNC': [ 70.0, 129.10],
'O C O ': [ 80.0, 126.00], 'O2C O2': [ 80.0, 126.00],
'O2P O2': [140.0, 119.90], 'O2P OH': [ 45.0, 108.23],
'O2P OS': [100.0, 108.23], 'OHP OS': [ 45.0, 102.60],
'OSP OS': [ 45.0, 102.60], 'P OSP ': [100.0, 120.50]} 

# AMBER94 molecular mechanics torsion parameters for atom type quartets
# where only the 2 central atom types are known:
# 1 -> (float) vn/2 [kcal/mol], rotation barrier height
# 2 -> (float) gamma [degrees], barrier minimum offset angle
# 3 -> (int) n [unitless], frequency of barrier
# 4 -> (int) paths [unitless], number of unique torsion paths
torsion23_params = {          'C CA': [14.50, 180.0, 2, 4],
'C CB': [12.00, 180.0, 2, 4], 'C CM': [ 8.70, 180.0, 2, 4],
'C CT': [ 0.00,   0.0, 2, 4], 'C N ': [10.00, 180.0, 2, 4],
'C N*': [ 5.80, 180.0, 2, 4], 'C NA': [ 5.40, 180.0, 2, 4],
'C NC': [ 8.00, 180.0, 2, 2], 'C OH': [ 1.80, 180.0, 2, 2],
'C*CB': [ 6.70, 180.0, 2, 4], 'C*CT': [ 0.00,   0.0, 2, 6],
'C*CW': [26.10, 180.0, 2, 4], 'CACA': [14.50, 180.0, 2, 4],
'CACB': [14.00, 180.0, 2, 4], 'CACM': [10.20, 180.0, 2, 4],
'CACN': [14.50, 180.0, 2, 4], 'CACT': [ 0.00,   0.0, 2, 6],
'CAN2': [ 9.60, 180.0, 2, 4], 'CANA': [ 6.00, 180.0, 2, 2],
'CANC': [ 9.60, 180.0, 2, 2], 'CBCB': [21.80, 180.0, 2, 4],
'CBCN': [12.00, 180.0, 2, 4], 'CBN*': [ 6.60, 180.0, 2, 4],
'CBNB': [ 5.10, 180.0, 2, 2], 'CBNC': [ 8.30, 180.0, 2, 2],
'CCCT': [ 0.00,   0.0, 2, 6], 'CCCV': [20.60, 180.0, 2, 4],
'CCCW': [21.50, 180.0, 2, 4], 'CCNA': [ 5.60, 180.0, 2, 4],
'CCNB': [ 4.80, 180.0, 2, 2], 'CKN*': [ 6.80, 180.0, 2, 4],
'CKNB': [20.00, 180.0, 2, 2], 'CMCM': [26.60, 180.0, 2, 4],
'CMCT': [ 0.00,   0.0, 3, 6], 'CMN*': [ 7.40, 180.0, 2, 4],
'CNNA': [ 6.10, 180.0, 2, 4], 'CQNC': [13.60, 180.0, 2, 2],
'CRNA': [ 9.30, 180.0, 2, 4], 'CRNB': [10.00, 180.0, 2, 2],
'CTCT': [ 1.40,   0.0, 3, 9], 'CTN ': [ 0.00,   0.0, 2, 6],
'CTN*': [ 0.00,   0.0, 2, 6], 'CTN2': [ 0.00,   0.0, 3, 6],
'CTN3': [ 1.40,   0.0, 3, 9], 'CTOH': [ 0.50,   0.0, 3, 3],
'CTOS': [ 1.15,   0.0, 3, 3], 'CTS ': [ 1.00,   0.0, 3, 3],
'CTSH': [ 0.75,   0.0, 3, 3], 'CVNB': [ 4.80, 180.0, 2, 2],
'CWNA': [ 6.00, 180.0, 2, 4], 'OHP ': [ 0.75,   0.0, 3, 3],
'OSP ': [ 4.80, 180.0, 2, 2]}

# AMBER94 molecular mechanics torsion parameters for atom type quartets
# where all 4 atom types are known (see above for parameter descriptions).
torsion1234_params = {
'C N CTC ': [[ 0.00,   0.0, -4, 1], [ 0.00, 180.0, -3, 1],
             [ 0.20, 180.0, -2, 1], [ 0.00, 180.0,  1, 1]],
'CTCTC N ': [[ 0.10,   0.0, -4, 1], [ 0.00,   0.0, -3, 1],
             [ 0.07,   0.0, -2, 1], [ 0.00, 180.0,  1, 1]],
'CTCTN C ': [[ 0.50, 180.0, -4, 1], [ 0.15, 180.0, -3, 1],
             [ 0.00, 180.0, -2, 1], [ 0.53,   0.0,  1, 1]],
'CTCTOSCT': [[0.383,   0.0, -3, 1], [ 0.10, 180.0,  2, 1]],
'CTS S CT': [[ 0.60,   0.0,  3, 1], [ 3.50,   0.0, -2, 1]],
'H N C O ': [[ 2.50, 180.0, -2, 1], [ 2.00,   0.0,  1, 1]],
'N CTC N ': [[ 0.40, 180.0, -4, 1], [ 0.00,   0.0, -3, 1],
             [ 1.35, 180.0, -2, 1], [ 0.75, 180.0,  1, 1]],
'OHCTCTOH': [[0.144,   0.0, -3, 1], [ 1.00,   0.0,  2, 1]],
'OHP OSCT': [[ 0.25,   0.0, -3, 1], [ 1.20,   0.0,  2, 1]],
'OSCTCTOH': [[0.144,   0.0, -3, 1], [ 1.00,   0.0,  2, 1]],
'OSCTCTOS': [[0.144,   0.0, -3, 1], [ 1.00,   0.0,  3, 1]],
'OSCTN*CK': [[ 0.50, 180.0, -2, 1], [ 2.50,   0.0,  1, 1]],
'OSCTN*CM': [[ 0.50, 180.0, -2, 1]],
'OSP OSCT': [[ 0.25,   0.0, -3, 1], [ 1.20,   0.0, 2.0, 1]],
'S CTN*CM': [[ 2.50,   0.0,  1, 1]]}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where only final 2 atom types are known:
# 1 -> (float) vn/2 [kcal/mol], rotation barrier height
oop34_params = {'C O ': 10.5, 'CAH4':  1.1, 'CAH5':  1.1, 'CAHA':  1.1,
  'CKH5':  1.1, 'CMH4':  1.1, 'CMHA':  1.1, 'CQH5':  1.1, 'CRH5':  1.1,
  'CVH4':  1.1, 'CWH4':  1.1, 'N H ':  1.0, 'N2H ':  1.0, 'NAH ':  1.0}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where only final 3 atom types are known (see above for
# parameter descriptions).
oop234_params = {'CTN CT':  1.0, 'N2CAN2': 10.5, 'O2C O2': 10.5}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where al 4 atom types are known (see above for parameter
# descriptions).
oop1234_params = {'CACAC OH':  1.1, 'CACACACT':  1.1, 'CBNCCAN2':  1.1,
'CKCBN*CT':  1.0, 'CMC CMCT':  1.1, 'CMC N*CT':  1.0, 'CTCMCMC' :  1.1,
'CWCBC*CT':  1.1, 'NCCMCAN2':  1.1, 'NACVCCCT':  1.1, 'NACWCCCT':  1.1,
'NANCCAN2':  1.1, 'NBCWCCCT':  1.1}

def get_vdw_param(at_type):
    """Find van der waals parameter for specified AMBER mm atom type.
    
    Args:
        at_type (str): AMBER mm atom type.
    
    Returns:
        vdw (float, float): van der waals parameters for atom type:
            ro/2 and eps.
    """
    try:
        vdw = vdw_params[at_type]
    except KeyError:
        print('Error: atom type (%s) not found!' % (at_type))
        sys.exit()
    return vdw

def get_at_mass(element):
    """Find the mass of an atom of a given element (periodic table avg).
    
    Args:
        element (str): atomic element string.
    
    Returns:
        at_mass (float): average periodic table mass of element.
    """
    try:
        at_mass = at_masses[element]
    except KeyError:
        print('Error: atomic mass for element (%s) not found!' % (element))
        sys.exit()
    return at_mass

def get_cov_rad(element):
    """Find the covalent radius of an atom of a given element.
    
    Args:
        element (str): atomic element string.
    
    Returns:
        cov_rad (float): covalent radius [Angstrom] of atom.
    """
    try:
        cov_rad = cov_radii[element]
    except KeyError:
        print('Error: covalent radius for element (%s) not found!' % (element))
        sys.exit()
    return cov_rad

def get_bond_param(at1_type, at2_type):
    """Find bond parameters for 2 AMBER94 mm atom types.
    
    Args:
        at1_type (str): atom 1 AMBER94 mm atom type.
        at2_type (str): atom 2 AMBER94 mm atom type.
    
    Returns:
        bond (float, float): AMBER94 mm bond parameters for atom types:
            k_b [kcal/(mol*A^2)] and r_eq [Angstrom].
    """
    try:
        bond_types = '%-2s%-2s' % (at1_type, at2_type)
        bond = bond_params[bond_types]
    except KeyError:
        try:
            bond_types = '%-2s%-2s' % (at2_type, at1_type)
            bond = bond_params[bond_types]
        except KeyError:
            print('Error: bond type (%s, %s) not recognized!' % (at1_type,
                at2_type))
            sys.exit()
    return bond

def get_angle_param(at1_type, at2_type, at3_type):
    """Find angle parameters for 3 AMBER94 mm atom types.
    
    Args:
        at1_type (str): atom 1 AMBER94 mm atom type.
        at2_type (str): atom 2 AMBER94 mm atom type.
        at3_type (str): atom 3 AMBER94 mm atom type.

    Returns:
        angle (float, float): AMBER94 mm angle parameters for atom types:
            k_a [kcal/(mol*rad^2)] and a_eq [degrees].
    """
    try:
        angle_types = '%-2s%-2s%-2s' % (at1_type, at2_type, at3_type)
        angle = angle_params[angle_types]
    except KeyError:
        try:
            angle_types = '%-2s%-2s%-2s' % (at3_type, at2_type, at1_type)
            angle = angle_params[angle_types]
        except KeyError:
            print('Error: angle type (%s, %s, %s) not recognized!' % (
                at1_type, at2_type, at3_type))
            sys.exit()
    return angle

def get_torsion_param(at1_type, at2_type, at3_type, at4_type):
    """Find torsion parameters for 4 AMBER94 mm atom types:
    
    Args:
        at1_type (str): atom 1 AMBER94 mm atom type.
        at2_type (str): atom 2 AMBER94 mm atom type.
        at3_type (str): atom 3 AMBER94 mm atom type.
        at4_type (str): atom 4 AMBER94 mm atom type.

    Returns:
        torsion (float, float, int, int): AMBER94 mm torsion parameters
            for atom types: vn/2 [kcal/mol], gamma [degrees], n [],
            and paths [].
    """
    torsion = []
  
    # two-atom torsion potentials
    try: 
        torsion23_types = '%-2s%-2s' % (at2_type, at3_type)
        torsion23 = torsion23_params[torsion23_types]
        torsion.append(torsion23)
    except KeyError:
        try:
            torsion23_types = '%-2s%-2s' % (at3_type, at2_type)
            torsion23 = torsion23_params[torsion23_types]
            torsion.append(torsion23)
        except KeyError:
            print('Error: torsion (X-%s-%s-X) not recognized!' % (at2_type,
                at3_type))
            sys.exit()
  
    # four-atom torsion potentials
    try:
        torsion1234_types = '%-2s%-2s%-2s%-2s' % (at1_type, at2_type,
            at3_type, at4_type)
        torsion1234 = torsion1234_params[torsion1234_types]
    except KeyError:
        try:
            torsion1234_types = '%-2s%-2s%-2s%-2s' % (at4_type, at3_type,
                at2_type, at1_type)
            torsion1234 = torsion1234_params[torsion1234_types]
        except KeyError:
            torsion1234 = []
    for i in range(len(torsion1234)):
        torsion.append(torsion1234[i])
  
    return torsion[0]

def get_outofplane_param(at1_type, at2_type, at3_type, at4_type):
    """Find outofplane parameters for 4 AMBER94 mm atom types:
    
    Args:
        at1_type (str): atom 1 AMBER94 mm atom type.
        at2_type (str): atom 2 AMBER94 mm atom type.
        at3_type (str): atom 3 AMBER94 mm atom type.
        at4_type (str): atom 4 AMBER94 mm atom type.

    Returns:
        torsion (float): AMBER94 mm outofplane parameters for atom types:
            vn/2 [kcal/mol].
    """
    oop = 0.0
  
    # two-atom out-of-plane potentials
    try:
        oop34_types = '%-2s%-2s' % (at3_type, at4_type)
        oop34 = oop34_params[oop34_types]
    except KeyError:
        try:
            oop34_types = '%-2s%-2s' % (at2_type, at1_type)
            oop34 = oop34_params[oop34_types]
        except KeyError:
            oop34 = []
    if (len(oop34) > 0):
        oop = oop34
  
    # three-atom out-of-plane potentials
    try:
        oop234_types = '%-2s%-2s%-2s' % (at2_type, at3_type, at4_type)
        oop234 = oop234_params[oop234_types]
    except KeyError:
        try:
            oop234_types = '%-2s%-2s%-2s' % (at3_type, at2_type, at1_type)
            oop234 = oop234_params[oop234_types]
        except KeyError:
            oop234 = []
    if (len(oop234) > 0):
        oop = oop234
  
    # four-atom out-of-plane potentials
    try:
        oop1234_types = '%-2s%-2s%-2s%-2s' % (at1_type, at2_type, at3_type,
            at4_type)
        oop1234 = oop1234_params[oop1234_types]
    except KeyError:
        try:
            oop1234_types = '%-2s%-2s%-2s%-2s' % (at4_type, at3_type,
                at2_type, at1_type)
            oop1234 = oop1234_params[oop1234_types]
        except KeyError:
            oop1234 = []
    if (len(oop1234) > 0):
        oop = oop1234
  
    return oop

# end of module

