"""Functions and dictionaries for AMBER94 molecular mechanics parameters.

Includes parameters for atomic mass, atomic covalent radii, atomic van der
waals, mm bonds, mm angles, mm torsions, and mm outofplanes.
"""

# relative atomic masses of elements (in atomic mass units [g/mol]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
_ATOMIC_MASSES = {
    'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067, 'F' : 18.9984,
    'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040, 'I' : 126.904,
    'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100, 'Be': 9.01218,
    'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815, 'Si': 28.0855,
    'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670, 'V' : 50.9415,
    'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332, 'Ni': 58.6934,
    'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400, 'As': 74.9216,
    'Se': 78.9600, 'Kr': 83.7980, 'X' : 0.00000}

# covalent (or ionic) radii by atomic element [Angstroms] from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
_COVALENT_RADII = {
    'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71, 'P' : 1.10,
    'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30, 'Ne': 0.84,
    'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02, 'Mg': 0.72,
    'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75, 'Ti': 0.86,
    'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64, 'Ni': 0.55,
    'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22, 'Se': 1.17,
    'Kr': 1.03, 'X' : 0.00}

# AMBER94 molecular mechanics van der waals parameters for atom types:
#   1 -> (float) ro/2 [Angstrom], van der waals radius (divided by 2)
#   2 -> (float) eps [kcal/mol], van der waals attraction magnitude
_VAN_DER_WAALS_PARAMETERS = {
    'C' : (1.9080, 0.0860), 'CA': (1.9080, 0.0860), 'CM': (1.9080, 0.0860),
    'CC': (1.9080, 0.0860), 'CV': (1.9080, 0.0860), 'CW': (1.9080, 0.0860),
    'CR': (1.9080, 0.0860), 'CB': (1.9080, 0.0860), 'C*': (1.9080, 0.0860),
    'CN': (1.9080, 0.0860), 'CK': (1.9080, 0.0860), 'CQ': (1.9080, 0.0860),
    'CT': (1.9080, 0.1094), 'CS': (3.3950, 0.0000806), 'F' : (1.7500, 0.0610),
    'H' : (0.6000, 0.0157), 'H1': (1.3870, 0.0157), 'H2': (1.2870, 0.0157),
    'H3': (1.1870, 0.0157), 'H4': (1.4090, 0.0150), 'H5': (1.3590, 0.0150), 
    'HA': (1.4590, 0.0150), 'HC': (1.4870, 0.0157), 'HO': (0.0001, 0.0000), 
    'HP': (1.1000, 0.0157), 'HS': (0.6000, 0.0157), 'HW': (0.0001, 0.0000),
    'Li': (1.1370, 0.0183), 'IP': (1.8680, 0.00277), 'Li': (1.1370, 0.0183),
    'N' : (1.8240, 0.1700), 'K' : (2.6580, 0.000328), 'O' : (1.6612, 0.2100),
    'O2': (1.6612, 0.2100), 'OH': (1.7210, 0.2104), 'OS': (1.6837, 0.1700),
    'OW': (1.7683, 0.1520), 'P' : (2.1000, 0.2000), 'S' : (2.0000, 0.2500),
    'SH': (2.0000, 0.2500), 'N' : (1.8240, 0.1700), 'NA': (1.8240, 0.1700),
    'NB': (1.8240, 0.1700), 'NC': (1.8240, 0.1700), 'N*': (1.8240, 0.1700),
    'N2': (1.8240, 0.1700), 'N3': (1.8750, 0.1700), 'He': (1.5800, 0.0112),
    'Ar': (1.8436, 0.4466), 'HH': (1.5000, 0.0000), 'X' : (0.0000, 0.0000)}

# AMBER94 molecular mechanics bond parameters for atom type pairs:
#   1 -> (float) k_b [kcal/(mol*A^2)], bond spring constant
#   2 -> (float) r_eq [Angstrom], equilibrium bond length
_BOND_LENGTH_PARAMETERS = {
    ('C' , 'CA'): (469.0, 1.409), ('C' , 'CB'): (447.0, 1.419),
    ('C' , 'CM'): (410.0, 1.444), ('C' , 'CT'): (317.0, 1.522),
    ('C' ,  'N'): (490.0, 1.335), ('C' , 'N*'): (424.0, 1.383),
    ('C' , 'NA'): (418.0, 1.388), ('C' , 'NC'): (457.0, 1.358),
    ('C' , 'O' ): (570.0, 1.229), ('C' , 'O2'): (656.0, 1.250),
    ('C' , 'OH'): (450.0, 1.364), ('C*', 'CB'): (388.0, 1.459),
    ('C*', 'CT'): (317.0, 1.495), ('C*', 'CW'): (546.0, 1.352),
    ('C*', 'HC'): (367.0, 1.080), ('CA', 'CA'): (469.0, 1.400),
    ('CA', 'CB'): (469.0, 1.404), ('CA', 'CM'): (427.0, 1.433),
    ('CA', 'CN'): (469.0, 1.400), ('CA', 'CT'): (317.0, 1.510),
    ('CA', 'H4'): (367.0, 1.080), ('CA', 'HA'): (367.0, 1.080),
    ('CA', 'N2'): (481.0, 1.340), ('CA', 'NA'): (427.0, 1.381),
    ('CA', 'NC'): (483.0, 1.339), ('CB', 'CB'): (520.0, 1.370),
    ('CB', 'CN'): (447.0, 1.419), ('CB', 'N*'): (436.0, 1.374),
    ('CB', 'NB'): (414.0, 1.391), ('CB', 'NC'): (461.0, 1.391),
    ('CC', 'CT'): (317.0, 1.504), ('CC', 'CV'): (512.0, 1.375),
    ('CC', 'CW'): (518.0, 1.371), ('CC', 'NA'): (422.0, 1.385),
    ('CC', 'NB'): (410.0, 1.394), ('CK', 'H5'): (367.0, 1.080),
    ('CK', 'N*'): (440.0, 1.371), ('CK', 'NB'): (529.0, 1.304),
    ('CM', 'CM'): (549.0, 1.350), ('CM', 'CT'): (317.0, 1.510),
    ('CM', 'H4'): (367.0, 1.080), ('CM', 'H5'): (367.0, 1.080),
    ('CM', 'HA'): (367.0, 1.080), ('CM', 'N*'): (448.0, 1.365),
    ('CN', 'NA'): (428.0, 1.380), ('CQ', 'H5'): (367.0, 1.080),
    ('CQ', 'NC'): (502.0, 1.324), ('CR', 'H5'): (367.0, 1.080),
    ('CR', 'NA'): (477.0, 1.343), ('CR', 'NB'): (488.0, 1.335),
    ('CT', 'CT'): (310.0, 1.526), ('CT', 'F' ): (367.0, 1.380),
    ('CT', 'H1'): (340.0, 1.090), ('CT', 'H2'): (340.0, 1.090),
    ('CT', 'H3'): (340.0, 1.090), ('CT', 'HC'): (340.0, 1.090),
    ('CT', 'HP'): (340.0, 1.090), ('CT', 'N' ): (337.0, 1.449),
    ('CT', 'N*'): (337.0, 1.475), ('CT', 'N2'): (337.0, 1.350),
    ('CT', 'N3'): (367.0, 1.471), ('CT', 'OH'): (320.0, 1.410),
    ('CT', 'OS'): (320.0, 1.410), ('CT', 'S' ): (227.0, 1.810),
    ('CT', 'SH'): (237.0, 1.810), ('CV', 'H4'): (367.0, 1.080),
    ('CV', 'NB'): (410.0, 1.394), ('CW', 'H4'): (367.0, 1.080),
    ('CW', 'NA'): (410.0, 1.394), ('H' , 'N' ): (434.0, 1.010),
    ('H' , 'N*'): (434.0, 1.010), ('H' , 'N2'): (434.0, 1.010),
    ('H' , 'N3'): (434.0, 1.010), ('H' , 'NA'): (434.0, 1.010),
    ('HO', 'OH'): (553.0, 0.960), ('HO', 'OS'): (553.0, 0.960),
    ('HS', 'SH'): (274.0, 1.336), ('O2', 'P' ): (525.0, 1.480),
    ('OH', 'P' ): (230.0, 1.610), ('OS', 'P' ): (230.0, 1.610),
    ('OW', 'HW'): (553.0,0.9572), ('S' , 'S' ): (166.0, 2.038),
    ('HH', 'HH'): (100.0, 0.740), ('X' , 'X' ): (  0.0, 0.000)}

# AMBER94 molecular mechanics angle parameters for atom type triplets:
#   1 -> (float) k_a [kcal/(mol*rad^2)], angle spring constant
#   2 -> (float) a_eq [degrees], equilibrium bond angle
_BOND_ANGLE_PARAMETERS= {
    ('C' , 'CA', 'CA'): ( 63.0, 120.00), ('C' , 'CA', 'HA'): ( 35.0, 120.00),
    ('C' , 'CB', 'CB'): ( 63.0, 119.20), ('C' , 'CB', 'NB'): ( 70.0, 130.00),
    ('C' , 'CM', 'CM'): ( 63.0, 120.70), ('C' , 'CM', 'CT'): ( 70.0, 119.70),
    ('C' , 'CM', 'H4'): ( 35.0, 119.70), ('C' , 'CM', 'HA'): ( 35.0, 119.70),
    ('C' , 'CT', 'CT'): ( 63.0, 111.10), ('C' , 'CT', 'H1'): ( 50.0, 109.50),
    ('C' , 'CT', 'HC'): ( 50.0, 109.50), ('C' , 'CT', 'HP'): ( 50.0, 109.50),
    ('C' , 'CT', 'N' ): ( 63.0, 110.10), ('C' , 'CT', 'N3'): ( 80.0, 111.20),
    ('C' , 'N' , 'CT'): ( 50.0, 121.90), ('C' , 'N' , 'H' ): ( 30.0, 120.00),
    ('C' , 'N*', 'CM'): ( 70.0, 121.60), ('C' , 'N*', 'CT'): ( 70.0, 117.60),
    ('C' , 'N*', 'H' ): ( 30.0, 119.20), ('C' , 'NA', 'C' ): ( 70.0, 126.40),
    ('C' , 'NA', 'CA'): ( 70.0, 125.20), ('C' , 'NA', 'H' ): ( 30.0, 116.80),
    ('C' , 'NC', 'CA'): ( 70.0, 120.50), ('C' , 'OH', 'HO'): ( 35.0, 113.00),
    ('C*', 'CB', 'CA'): ( 63.0, 134.90), ('C*', 'CB', 'CN'): ( 63.0, 108.80),
    ('C*', 'CT', 'CT'): ( 63.0, 115.60), ('C*', 'CT', 'HC'): ( 50.0, 109.50),
    ('C*', 'CW', 'H4'): ( 35.0, 120.00), ('C*', 'CW', 'NA'): ( 70.0, 108.70),
    ('CA', 'C' , 'CA'): ( 63.0, 120.00), ('CA', 'C' , 'OH'): ( 70.0, 120.00),
    ('CA', 'CA', 'CA'): ( 63.0, 120.00), ('CA', 'CA', 'CB'): ( 63.0, 120.00),
    ('CA', 'CA', 'CN'): ( 63.0, 120.00), ('CA', 'CA', 'CT'): ( 70.0, 120.00),
    ('CA', 'CA', 'H4'): ( 35.0, 120.00), ('CA', 'CA', 'HA'): ( 35.0, 120.00),
    ('CA', 'CB', 'CB'): ( 63.0, 117.30), ('CA', 'CB', 'CN'): ( 63.0, 116.20),
    ('CA', 'CB', 'NB'): ( 70.0, 132.40), ('CA', 'CM', 'CM'): ( 63.0, 117.00),
    ('CA', 'CM', 'H4'): ( 35.0, 123.30), ('CA', 'CM', 'HA'): ( 35.0, 123.30),
    ('CA', 'CN', 'CB'): ( 63.0, 122.70), ('CA', 'CN', 'NA'): ( 70.0, 132.80),
    ('CA', 'CT', 'CT'): ( 63.0, 114.00), ('CA', 'CT', 'HC'): ( 50.0, 109.50),
    ('CA', 'N2', 'CT'): ( 50.0, 123.20), ('CA', 'N2', 'H' ): ( 35.0, 120.00),
    ('CA', 'NA', 'H' ): ( 30.0, 118.00), ('CA', 'NC', 'CB'): ( 70.0, 112.20),
    ('CA', 'NC', 'CQ'): ( 70.0, 118.60), ('CB', 'C' , 'NA'): ( 70.0, 111.30),
    ('CB', 'C' , 'O' ): ( 80.0, 128.80), ('CB', 'C*', 'CT'): ( 70.0, 128.60),
    ('CB', 'C*', 'CW'): ( 63.0, 106.40), ('CB', 'CA', 'H4'): ( 35.0, 120.00),
    ('CB', 'CA', 'HA'): ( 35.0, 120.00), ('CB', 'CA', 'N2'): ( 70.0, 123.50),
    ('CB', 'CA', 'NC'): ( 70.0, 117.30), ('CB', 'CB', 'N*'): ( 70.0, 106.20),
    ('CB', 'CB', 'NB'): ( 70.0, 110.40), ('CB', 'CB', 'NC'): ( 70.0, 127.70),
    ('CB', 'CN', 'NA'): ( 70.0, 104.40), ('CB', 'N*', 'CK'): ( 70.0, 105.40),
    ('CB', 'N*', 'CT'): ( 70.0, 125.80), ('CB', 'N*', 'H' ): ( 30.0, 125.80),
    ('CB', 'NB', 'CK'): ( 70.0, 103.80), ('CB', 'NC', 'CQ'): ( 70.0, 111.00),
    ('CC', 'CT', 'CT'): ( 63.0, 113.10), ('CC', 'CT', 'HC'): ( 50.0, 109.50),
    ('CC', 'CV', 'H4'): ( 35.0, 120.00), ('CC', 'CV', 'NB'): ( 70.0, 120.00),
    ('CC', 'CW', 'H4'): ( 35.0, 120.00), ('CC', 'CW', 'NA'): ( 70.0, 120.00),
    ('CC', 'NA', 'CR'): ( 70.0, 120.00), ('CC', 'NA', 'H' ): ( 30.0, 120.00),
    ('CC', 'NB', 'CR'): ( 70.0, 117.00), ('CK', 'N*', 'CT'): ( 70.0, 128.80),
    ('CK', 'N*', 'H' ): ( 30.0, 128.80), ('CM', 'C' , 'NA'): ( 70.0, 114.10),
    ('CM', 'C' , 'O' ): ( 80.0, 125.30), ('CM', 'CA', 'N2'): ( 70.0, 120.10),
    ('CM', 'CA', 'NC'): ( 70.0, 121.50), ('CM', 'CM', 'CT'): ( 70.0, 119.70),
    ('CM', 'CM', 'H4'): ( 35.0, 119.70), ('CM', 'CM', 'HA'): ( 35.0, 119.70),
    ('CM', 'CM', 'N*'): ( 70.0, 121.20), ('CM', 'CT', 'HC'): ( 50.0, 109.50),
    ('CM', 'N*', 'CT'): ( 70.0, 121.20), ('CM', 'N*', 'H' ): ( 30.0, 121.20),
    ('CN', 'CA', 'HA'): ( 35.0, 120.00), ('CN', 'NA', 'CW'): ( 70.0, 111.60),
    ('CN', 'NA', 'H' ): ( 30.0, 123.10), ('CR', 'NA', 'CW'): ( 70.0, 120.00),
    ('CR', 'NA', 'H' ): ( 30.0, 120.00), ('CR', 'NB', 'CV'): ( 70.0, 117.00),
    ('CT', 'C' , 'N' ): ( 70.0, 116.60), ('CT', 'C' , 'O' ): ( 80.0, 120.40),
    ('CT', 'C' , 'O2'): ( 70.0, 117.00), ('CT', 'C*', 'CW'): ( 70.0, 125.00),
    ('CT', 'CC', 'CV'): ( 70.0, 120.00), ('CT', 'CC', 'CW'): ( 70.0, 120.00),
    ('CT', 'CC', 'NA'): ( 70.0, 120.00), ('CT', 'CC', 'NB'): ( 70.0, 120.00),
    ('CT', 'CT', 'CT'): ( 40.0, 109.50), ('CT', 'CT', 'H1'): ( 50.0, 109.50),
    ('CT', 'CT', 'H2'): ( 50.0, 109.50), ('CT', 'CT', 'HC'): ( 50.0, 109.50),
    ('CT', 'CT', 'HP'): ( 50.0, 109.50), ('CT', 'CT', 'N' ): ( 80.0, 109.70),
    ('CT', 'CT', 'N*'): ( 50.0, 109.50), ('CT', 'CT', 'N2'): ( 80.0, 111.20),
    ('CT', 'CT', 'N3'): ( 80.0, 111.20), ('CT', 'CT', 'OH'): ( 50.0, 109.50),
    ('CT', 'CT', 'OS'): ( 50.0, 109.50), ('CT', 'CT', 'S' ): ( 50.0, 114.70),
    ('CT', 'CT', 'SH'): ( 50.0, 108.60), ('CT', 'N' , 'CT'): ( 50.0, 118.00),
    ('CT', 'N' , 'H' ): ( 30.0, 118.04), ('CT', 'N2', 'H' ): ( 35.0, 118.40),
    ('CT', 'N3', 'H' ): ( 50.0, 109.50), ('CT', 'OH', 'HO'): ( 55.0, 108.50),
    ('CT', 'OS', 'CT'): ( 60.0, 109.50), ('CT', 'OS', 'P' ): (100.0, 120.50),
    ('CT', 'S' , 'CT'): ( 62.0,  98.90), ('CT', 'S' , 'S' ): ( 68.0, 103.70),
    ('CT', 'SH', 'SH'): ( 43.0,  96.00), ('CV', 'CC', 'NA'): ( 70.0, 120.00),
    ('CW', 'CC', 'NA'): ( 70.0, 120.00), ('CW', 'CC', 'NB'): ( 70.0, 120.00),
    ('CW', 'NA', 'H' ): ( 30.0, 120.00), ('F' , 'CT', 'F' ): ( 77.0, 109.10),
    ('F' , 'CT', 'H1'): ( 35.0, 109.50), ('H' , 'N' , 'H' ): ( 35.0, 120.00),
    ('H' , 'N2', 'H' ): ( 35.0, 120.00), ('H' , 'N3', 'H' ): ( 35.0, 109.50),
    ('H1', 'CT', 'H1'): ( 35.0, 109.50), ('H1', 'CT', 'N' ): ( 50.0, 109.50),
    ('H1', 'CT', 'N*'): ( 50.0, 109.50), ('H1', 'CT', 'N2'): ( 50.0, 109.50),
    ('H1', 'CT', 'OH'): ( 50.0, 109.50), ('H1', 'CT', 'OS'): ( 50.0, 109.50),
    ('H1', 'CT', 'S' ): ( 50.0, 109.50), ('H1', 'CT', 'SH'): ( 50.0, 109.50),
    ('H2', 'CT', 'H2'): ( 35.0, 109.50), ('H2', 'CT', 'N*'): ( 50.0, 109.50),
    ('H2', 'CT', 'OS'): ( 50.0, 109.50), ('H4', 'CM', 'N*'): ( 35.0, 119.10),
    ('H4', 'CV', 'NB'): ( 35.0, 120.00), ('H4', 'CW', 'NA'): ( 35.0, 120.00),
    ('H5', 'CK', 'N*'): ( 35.0, 123.05), ('H5', 'CK', 'NB'): ( 35.0, 123.05),
    ('H5', 'CQ', 'NC'): ( 35.0, 115.45), ('H5', 'CR', 'NA'): ( 35.0, 120.00),
    ('H5', 'CR', 'NB'): ( 35.0, 120.00), ('HC', 'CT', 'HC'): ( 35.0, 109.50),
    ('HO', 'OH', 'P' ): ( 45.0, 108.50), ('HP', 'CT', 'HP'): ( 35.0, 109.50),
    ('HP', 'CT', 'N3'): ( 50.0, 109.50), ('HS', 'SH', 'HS'): ( 35.0,  92.07),
    ('HW', 'OW', 'HW'): (100.0, 104.52), ('N' , 'C' , 'O' ): ( 80.0, 122.90),
    ('N*', 'C' , 'NA'): ( 70.0, 115.40), ('N*', 'C' , 'NC'): ( 70.0, 118.60),
    ('N*', 'C' , 'O' ): ( 80.0, 120.90), ('N*', 'CB', 'NC'): ( 70.0, 126.20),
    ('N*', 'CK', 'NB'): ( 70.0, 113.90), ('N*', 'CT', 'OS'): ( 50.0, 109.50),
    ('N2', 'CA', 'N2'): ( 70.0, 120.00), ('N2', 'CA', 'NA'): ( 70.0, 116.00),
    ('N2', 'CA', 'NC'): ( 70.0, 119.30), ('NA', 'C' , 'O' ): ( 80.0, 120.60),
    ('NA', 'CA', 'NC'): ( 70.0, 123.30), ('NA', 'CR', 'NA'): ( 70.0, 120.00),
    ('NA', 'CR', 'NB'): ( 70.0, 120.00), ('NC', 'C' , 'O' ): ( 80.0, 122.50),
    ('NC', 'CQ', 'NC'): ( 70.0, 129.10), ('O' , 'C' , 'O' ): ( 80.0, 126.00),
    ('O2', 'C' , 'O2'): ( 80.0, 126.00), ('O2', 'P' , 'O2'): (140.0, 119.90),
    ('O2', 'P' , 'OH'): ( 45.0, 108.23), ('O2', 'P' , 'OS'): (100.0, 108.23),
    ('OH', 'P' , 'OS'): ( 45.0, 102.60), ('OS', 'P' , 'OS'): ( 45.0, 102.60),
    ('P' , 'OS', 'P' ): (100.0, 120.50), ('X' , 'X' , 'X' ): (  0.0,   0.00)} 

# AMBER94 molecular mechanics torsion parameters for atom type quartets
# where only the 2 central atom types are known:
#   1 -> (float) vn/2 [kcal/mol], rotation barrier height
#   2 -> (float) gamma [degrees], barrier minimum offset angle
#   3 -> (int) n [unitless], frequency of barrier
#   4 -> (int) paths [unitless], number of unique torsion paths
_TORSION_23_PARAMETERS = {
    ('C' , 'CA'): (14.50, 180.0, 2, 4), ('C' , 'CB'): (12.00, 180.0, 2, 4),
    ('C' , 'CM'): ( 8.70, 180.0, 2, 4), ('C' , 'CT'): ( 0.00,   0.0, 2, 4),
    ('C' , 'N' ): (10.00, 180.0, 2, 4), ('C' , 'N*'): ( 5.80, 180.0, 2, 4),
    ('C' , 'NA'): ( 5.40, 180.0, 2, 4), ('C' , 'NC'): ( 8.00, 180.0, 2, 2),
    ('C' , 'OH'): ( 1.80, 180.0, 2, 2), ('C*', 'CB'): ( 6.70, 180.0, 2, 4),
    ('C*', 'CT'): ( 0.00,   0.0, 2, 6), ('C*', 'CW'): (26.10, 180.0, 2, 4),
    ('CA', 'CA'): (14.50, 180.0, 2, 4), ('CA', 'CB'): (14.00, 180.0, 2, 4),
    ('CA', 'CM'): (10.20, 180.0, 2, 4), ('CA', 'CN'): (14.50, 180.0, 2, 4),
    ('CA', 'CT'): ( 0.00,   0.0, 2, 6), ('CA', 'N2'): ( 9.60, 180.0, 2, 4),
    ('CA', 'NA'): ( 6.00, 180.0, 2, 2), ('CA', 'NC'): ( 9.60, 180.0, 2, 2),
    ('CB', 'CB'): (21.80, 180.0, 2, 4), ('CB', 'CN'): (12.00, 180.0, 2, 4),
    ('CB', 'N*'): ( 6.60, 180.0, 2, 4), ('CB', 'NB'): ( 5.10, 180.0, 2, 2),
    ('CB', 'NC'): ( 8.30, 180.0, 2, 2), ('CC', 'CT'): ( 0.00,   0.0, 2, 6),
    ('CC', 'CV'): (20.60, 180.0, 2, 4), ('CC', 'CW'): (21.50, 180.0, 2, 4),
    ('CC', 'NA'): ( 5.60, 180.0, 2, 4), ('CC', 'NB'): ( 4.80, 180.0, 2, 2),
    ('CK', 'N*'): ( 6.80, 180.0, 2, 4), ('CK', 'NB'): (20.00, 180.0, 2, 2),
    ('CM', 'CM'): (26.60, 180.0, 2, 4), ('CM', 'CT'): ( 0.00,   0.0, 3, 6),
    ('CM', 'N*'): ( 7.40, 180.0, 2, 4), ('CN', 'NA'): ( 6.10, 180.0, 2, 4),
    ('CQ', 'NC'): (13.60, 180.0, 2, 2), ('CR', 'NA'): ( 9.30, 180.0, 2, 4),
    ('CR', 'NB'): (10.00, 180.0, 2, 2), ('CT', 'CT'): ( 1.40,   0.0, 3, 9),
    ('CT', 'N' ): ( 0.00,   0.0, 2, 6), ('CT', 'N*'): ( 0.00,   0.0, 2, 6),
    ('CT', 'N2'): ( 0.00,   0.0, 3, 6), ('CT', 'N3'): ( 1.40,   0.0, 3, 9),
    ('CT', 'OH'): ( 0.50,   0.0, 3, 3), ('CT', 'OS'): ( 1.15,   0.0, 3, 3),
    ('CT', 'S' ): ( 1.00,   0.0, 3, 3), ('CT', 'SH'): ( 0.75,   0.0, 3, 3),
    ('CV', 'NB'): ( 4.80, 180.0, 2, 2), ('CW', 'NA'): ( 6.00, 180.0, 2, 4),
    ('OH', 'P' ): ( 0.75,   0.0, 3, 3), ('OS', 'P' ): ( 4.80, 180.0, 2, 2)}

# AMBER94 molecular mechanics torsion parameters for atom type quartets
# where all 4 atom types are known (see above for parameter descriptions).
_TORSION_1234_PARAMETERS = {
    ('C' , 'N' , 'CT', 'C' ): [( 0.00,   0.0, -4, 1), ( 0.00, 180.0, -3, 1),
                               ( 0.20, 180.0, -2, 1), ( 0.00, 180.0,  1, 1)],
    ('CT', 'CT', 'C' , 'N' ): [( 0.10,   0.0, -4, 1), ( 0.00,   0.0, -3, 1),
                               ( 0.07,   0.0, -2, 1), ( 0.00, 180.0,  1, 1)],
    ('CT', 'CT', 'N' , 'C' ): [( 0.50, 180.0, -4, 1), ( 0.15, 180.0, -3, 1),
                               ( 0.00, 180.0, -2, 1), ( 0.53,   0.0,  1, 1)],
    ('CT', 'CT', 'OS', 'CT'): [(0.383,   0.0, -3, 1), ( 0.10, 180.0,  2, 1)],
    ('CT', 'S' , 'S' , 'CT'): [( 0.60,   0.0,  3, 1), ( 3.50,   0.0, -2, 1)],
    ('H' , 'N' , 'C' , 'O' ): [( 2.50, 180.0, -2, 1), ( 2.00,   0.0,  1, 1)],
    ('N' , 'CT', 'C' , 'N' ): [( 0.40, 180.0, -4, 1), ( 0.00,   0.0, -3, 1),
                               ( 1.35, 180.0, -2, 1), ( 0.75, 180.0,  1, 1)],
    ('OH', 'CT', 'CT', 'OH'): [(0.144,   0.0, -3, 1), ( 1.00,   0.0,  2, 1)],
    ('OH', 'P' , 'OS', 'CT'): [( 0.25,   0.0, -3, 1), ( 1.20,   0.0,  2, 1)],
    ('OS', 'CT', 'CT', 'OH'): [(0.144,   0.0, -3, 1), ( 1.00,   0.0,  2, 1)],
    ('OS', 'CT', 'CT', 'OS'): [(0.144,   0.0, -3, 1), ( 1.00,   0.0,  3, 1)],
    ('OS', 'CT', 'N*', 'CK'): [( 0.50, 180.0, -2, 1), ( 2.50,   0.0,  1, 1)],
    ('OS', 'CT', 'N*', 'CM'): [( 0.50, 180.0, -2, 1)],
    ('OS', 'P' , 'OS', 'CT'): [( 0.25,   0.0, -3, 1), ( 1.20,   0.0, 2.0, 1)],
    ('S' , 'CT', 'N*', 'CM'): [( 2.50,   0.0,  1, 1)]}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where only final 2 atom types are known:
# 1 -> (float) vn/2 [kcal/mol], rotation barrier height
_OUTOFPLANE_34_PARAMETERS = {
    ('C' , 'O' ): 10.5, ('CA', 'H4'):  1.1, ('CA', 'H5'):  1.1,
    ('CA', 'HA'):  1.1, ('CK', 'H5'):  1.1, ('CM', 'H4'):  1.1,
    ('CM', 'HA'):  1.1, ('CQ', 'H5'):  1.1, ('CR', 'H5'):  1.1,
    ('CV', 'H4'):  1.1, ('CW', 'H4'):  1.1, ('N' , 'H' ):  1.0,
    ('N2', 'H' ):  1.0, ('NA', 'H' ):  1.0}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where only final 3 atom types are known (see above for
# parameter descriptions).
_OUTOFPLANE_234_PARAMETERS= {
    ('CT', 'N' , 'CT'):  1.0, ('N2', 'CA', 'N2'): 10.5,
    ('O2', 'C' , 'O2'): 10.5}

# AMBER94 molecular mechanics outofplane parameters for atom type
# quartets where al 4 atom types are known (see above for parameter
# descriptions).
_OUTOFPLANE_1234_PARAMETERS = {
    ('CA', 'CA', 'C' , 'OH'): 1.1, ('CA', 'CA', 'CA', 'CT'): 1.1,
    ('CB', 'NC', 'CA', 'N2'): 1.1, ('CK', 'CB', 'N*', 'CT'): 1.0,
    ('CM', 'C' , 'CM', 'CT'): 1.1, ('CM', 'C' , 'N*', 'CT'): 1.0,
    ('CT', 'CM', 'CM', 'C' ): 1.1, ('CW', 'CB', 'C*', 'CT'): 1.1,
    ('NC', 'CM', 'CA', 'N2'): 1.1, ('NA', 'CV', 'CC', 'CT'): 1.1,
    ('NA', 'CW', 'CC', 'CT'): 1.1, ('NA', 'NC', 'CA', 'N2'): 1.1,
    ('NB', 'CW', 'CC', 'CT'): 1.1}

def GetElement(at_type):
  """Infer atomic element from atom type.

  If atom type is a single character, or second character is uppercase,
  return uppercase first letter. Otherwise, return capitalized first two
  characters.

  Args:
    at_type (str): Atom type.

  Returns:
    at_element (str): Atomic element.
  """
  if len(at_type) == 1 or not at_type[1].islower():
    return at_type[0].upper()
  else:
    return at_type[0:2].capitalize()


def GetMass(element):
  """Find the mass of an atom of a given element (periodic table avg).
  
  Args:
    element (str): atomic element string.
  
  Returns:
    at_mass (float): average periodic table mass [amu] of element.

  Raises:
    ValueError: If parameter not found for element.
  """
  if element in _ATOMIC_MASSES:
    return _ATOMIC_MASSES[element]
  else:
    raise ValueError('No atomic mass found for element: %s' % (element))


def GetCovRad(element):
  """Find the covalent radius of an atom of a given element.
  
  Args:
    element (str): atomic element string.
  
  Returns:
    cov_rad (float): covalent radius [Angstrom] of atom.

  Raises:
    ValueError: If parameter not found for element.
  """
  if element in _COVALENT_RADII:
    return _COVALENT_RADII[element]
  else:
    raise ValueError('No covalent radius found for element: %s' % (element))


def GetVdwParam(at_type):
  """Find van der waals parameter for specified AMBER mm atom type.
  
  Args:
    at_type (str): AMBER mm atom type.
  
  Returns:
    van der waals parameters (float, float) for atom type:
      * ro/2 [Angstrom] 0.0 --> ?
      * eps [kcal/mol] 0.0 --> ?

  Raises:
    ValueError: If parameters not found for atom type.
  """
  if at_type in _VAN_DER_WAALS_PARAMETERS:
    return _VAN_DER_WAALS_PARAMETERS[at_type]
  else:
    raise ValueError('No vdw param found for atom type: %s' % (at_type))


def GetBondParam(at1_type, at2_type):
  """Find bond parameters for 2 AMBER94 mm atom types.
  
  Args:
    at1_type (str): atom 1 AMBER94 mm atom type.
    at2_type (str): atom 2 AMBER94 mm atom type.
  
  Returns:
    bond (float, float) AMBER94 mm bond parameters for atom types:
        * k_b [kcal/(mol*A^2)] 0.0 --> ?
        * r_eq [Angstrom] 0.0 --> ?

  Raises:
    ValueError: If parameters not found for atom types.
  """
  b12_types = at1_type, at2_type
  b21_types = at2_type, at1_type
  if b12_types in _BOND_LENGTH_PARAMETERS:
    return _BOND_LENGTH_PARAMETERS[b12_types]
  elif b21_types in _BOND_LENGTH_PARAMETERS:
    return _BOND_LENGTH_PARAMETERS[b21_types]
  else:
    raise ValueError('No bond parameters found for atom type pair: '
                     '%s, %s' % (at1_type, at2_type))


def GetAngleParam(at1_type, at2_type, at3_type):
  """Find angle parameters for 3 AMBER94 mm atom types.
  
  Args:
    at1_type (str): atom 1 AMBER94 mm atom type.
    at2_type (str): atom 2 AMBER94 mm atom type.
    at3_type (str): atom 3 AMBER94 mm atom type.

  Returns:
    angle (float, float) AMBER94 mm angle parameters for atom types:
        * k_a [kcal/(mol*rad^2)] 0.0 --> ?
        * a_eq [degrees] 0.0 --> 180.0

  Raises:
    ValueError: If parameters not found for atom types.
  """
  a123_types = at1_type, at2_type, at3_type
  a321_types = at3_type, at2_type, at1_type
  if a123_types in _BOND_ANGLE_PARAMETERS:
    return _BOND_ANGLE_PARAMETERS[a123_types]
  elif a321_types in _BOND_ANGLE_PARAMETERS:
    return _BOND_ANGLE_PARAMETERS[a321_types]
  else:
    raise ValueError('No angle parameters found for atom type triplet: '
                     '%s, %s, %s' % (at1_type, at2_type, at3_type))


def GetTorsionParam(at1_type, at2_type, at3_type, at4_type):
  """Find torsion parameters for 4 AMBER94 mm atom types:
  
  Args:
    at1_type (str): atom 1 AMBER94 mm atom type.
    at2_type (str): atom 2 AMBER94 mm atom type.
    at3_type (str): atom 3 AMBER94 mm atom type.
    at4_type (str): atom 4 AMBER94 mm atom type.

  Returns:
    torsion (float, float, int, int)* AMBER94 mm torsion parameter array for
        atom types:
      * vn/2 [kcal/mol] 0.0 --> ?
      * gamma [degrees] -180.0 --> 180.0
      * n [] 1 --> 6
      * paths [] 1 --> 9

  Raises:
    ValueError: If parameters not found for atom types.
  """
  torsion = []

  # two-atom torsion potentials if found
  t23_types = at2_type, at3_type
  t32_types = at3_type, at2_type
  if t23_types in _TORSION_23_PARAMETERS:
    t23 = _TORSION_23_PARAMETERS[t23_types]
    torsion.append(t23)
  elif t32_types in _TORSION_23_PARAMETERS:
    t23 = _TORSION_23_PARAMETERS[t32_types]
    torsion.append(t23)
  else:
    raise ValueError('No torsion parameters found for central atom pair: '
                     '%s, %s' % (at2_type, at3_type))

  # four-atom torsion potentials if found
  t1234_types = at1_type, at2_type, at3_type, at4_type
  t4321_types = at4_type, at3_type, at2_type, at1_type
  t1234 = []
  if t1234_types in _TORSION_1234_PARAMETERS:
    t1234 = _TORSION_1234_PARAMETERS[t1234_types]
  elif t4321_types in _TORSION_1234_PARAMETERS:
    t1234 = _TORSION_1234_PARAMETERS[t4321_types]
  for i in range(len(t1234)):
    torsion.append(t1234[i])

  # return array of all found torsion parameters
  return torsion


def GetOutofplaneParam(at1_type, at2_type, at3_type, at4_type):
  """Find outofplane parameters for 4 AMBER94 mm atom types:
    
  Args:
    at1_type (str): atom 1 AMBER94 mm atom type.
    at2_type (str): atom 2 AMBER94 mm atom type.
    at3_type (str): atom 3 AMBER94 mm atom type.
    at4_type (str): atom 4 AMBER94 mm atom type.

  Returns:
    torsion (float): AMBER94 mm outofplane parameters for atom types:
        * vn/2 [kcal/mol] 0.0 --> ?
  """

  # return four-atom out-of-plane potential if found
  oop1234_types = at1_type, at2_type, at3_type, at4_type
  if oop1234_types in _OUTOFPLANE_1234_PARAMETERS:
    return _OUTOFPLANE_1234_PARAMETERS[oop1234_types]

  # return three-atom out-of-plane potential if found
  oop234_types = at2_type, at3_type, at4_type
  if oop234_types in _OUTOFPLANE_234_PARAMETERS:
    return _OUTOFPLANE_234_PARAMETERS[oop234_types]

  # return two-atom out-of-plane potential if found
  oop34_types = at3_type, at4_type
  if oop34_types in _OUTOFPLANE_34_PARAMETERS:
    return _OUTOFPLANE_34_PARAMETERS[oop34_types]

  # return zero if out-of-plane potential not found
  return 0.0
