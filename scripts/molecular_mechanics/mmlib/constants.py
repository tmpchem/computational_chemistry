"""Constant values for use in other modules within mmlib.

Contains physical constants, arbitrary values, and dictionary mappings for
plotting (listed in alphabetical order).
"""

import math

# Conversion of acceleration from [kcal/(A*g) to [A/(ps^2)].
ACCCONV = 418.400000

# Threshold beyond covalent radii sum to determine bond cutoff.
BONDTHRESHOLD = 1.2

# Conversion of electrostatic energy from [ceu] to [kcal/mol].
CEU2KCAL = 332.06375

# Conversion from degrees to radians
DEG2RAD = math.pi / 180.0

# Boltzmann constant [kcal/(mol*K)].
KB = 0.001987204

# Conversion from [kcal*A^3/mol] to [Pa] for pressure.
KCALAMOL2PA = 69476.95

# Conversion of kinetic energy from [amu*A^2/ps^2] to [kcal/mol].
KIN2KCAL = 0.00239005736

# Number of Cartesian dimensions
NUMDIM = 3

# Displacement distance [Angstrom] for numerical gradient.
NUMDISP = 1.0E-6

# Number of step iterations to take in optimization binary line search.
NUMLINESEARCHSTEPS = 7

# Default optimization criteria keyword dictionary.
# [delta_e, grad_rms, grad_max, disp_rms, disp_max]
OPTCRITERIAREFS = {
    'loose':     [1.0E-4,  1.0E-3, 2.0E-3, 1.0E-2, 2.0E-2],
    'default':   [1.0E-6,  1.0E-4, 2.0E-4, 1.0E-3, 2.0E-3],
    'tight':     [1.0E-8,  1.0E-5, 2.0E-5, 1.0E-4, 2.0E-4],
    'verytight': [1.0E-10, 1.0E-6, 2.0E-6, 1.0E-5, 2.0E-5]}

# Factor by which to adjust the initial line search step size between steps.
OPTSTEPADJUSTOR = math.sqrt(2)

# Fraction of image width which is covered by the plot field.
PERCENTIMAGEPLOT = 0.75

# Unit conversion between points and inches
POINTSPERINCH = 72

# Legend labels, line colors, and plotting priority for properties.
# [energy_term, print_priority, line_color, index]
PROPERTYDICTIONARY = {
    'e_total':    ['Total',      12, '#000000', 1],
    'e_kin':      ['Kinetic',    11, '#007D34', 2],
    'e_pot':      ['Potential',   1, '#C10020', 3],
    'e_nonbond':  ['Non-bonded',  7, '#0000FF', 4],
    'e_bonded':   ['Bonded',      2, '#FF6800', 5],
    'e_boundary': ['Boundary',   10, '#551A8B', 6],
    'e_vdw':      ['Vdw',         9, '#00BFFF', 7],
    'e_elst':     ['Elst',        8, '#EEC900', 8],
    'e_bond':     ['Bonds',       3, '#F08080', 9],
    'e_angle':    ['Angles',      4, '#90EE90', 10],
    'e_tors':     ['Torsions',    6, '#FF83FA', 11],
    'e_oop':      ['Outofplanes', 5, '#A9A9A9', 12]}

# Physical property keys for output file data labels.
PROPERTYKEYS = [
    'e_total', 'e_kin', 'e_pot', 'e_nonbond', 'e_bonded', 'e_boundary',
    'e_vdw', 'e_elst', 'e_bond', 'e_angle', 'e_tors', 'e_oop', 'temperature',
    'pressure']

# Conversion from radians to degrees.
RAD2DEG = 180.0 / math.pi

# Gas constant in units of [amu*A^2/(ps^2*K)].
RGAS = 0.83144598

# Dictionary of order-of-magnitude axis tick labels.
TICCHARS = {0: '', 1: 'k', 2: 'M', 3: 'B', 4: 'T', 5: 'P'}
