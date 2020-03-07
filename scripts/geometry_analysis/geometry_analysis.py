import sys, math, numpy

#############################################################################
#                       Welcome to geometry_analysis.py                     #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and outputs      #
# bond lengths, angles, torsions, out-of-planes, center of mass, and        #
# moment of inertia.                                                        #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 12/23/2015                                #
#############################################################################

## CONSTANTS ##

# threshold beyond average of covalent radii to determine bond cutoff
bond_thresh = 1.2

# conversion from radians to degrees and vice versa
rad2deg = 180.0 / math.pi
deg2rad = math.pi / 180.0

# threshold for degenerate principal moments of inertia (1 part per 10**n)
mom_thresh = 3

# absolute threshold for moment of inertia difference
mom_min = 10**(-mom_thresh)

# Planck's constant (J*s)
h = 6.62607 * 10**-34

# avogadro's number (mol^-1)
na = 6.02214 * 10**23

# speed of light (cm/s)
c = 2.99792 * 10**10

# cartesian indices
xyz = {0: 'X', 1: 'Y', 2: 'Z'}

# covalent (or ionic) radii by atomic element (Angstroms) from
# "Inorganic Chemistry" 3rd ed, Housecroft, Appendix 6, pgs 1013-1014
cov_rads = {  'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 'F' : 0.71,
  'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 'Br': 1.14, 'I' : 1.33, 'He': 0.30,
  'Ne': 0.84, 'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 'Na': 1.02,
  'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75,
  'Ti': 0.86, 'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 'Co': 0.64,
  'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 'Ga': 1.22, 'Ge': 1.22, 'As': 1.22,
  'Se': 1.17, 'Kr': 1.03, 'X' : 0.00}

# relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {    'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067,
  'F' : 18.9984, 'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
  'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
  'Be': 9.01218, 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
  'Si': 28.0855, 'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
  'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
  'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
  'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X' : 0.00000}

## IO FUNCTIONS ##

# read file data into a 2-d array
def get_file_string_array(file_name):
    try:
        file = open(file_name, "r")
    except IOError:
        print('Error: file (%s) not found!\n' % (file_name))
        sys.exit()
    lines = file.readlines()
    file.close()
    array = []
    for line in lines:
        array.append(line.split())
    return array

# read in geometry from xyz file
def get_geom(xyz_file_name):
    xyz_array = get_file_string_array(xyz_file_name)
    n_atoms = int(xyz_array[0][0])
    at_types = ['' for i in range(n_atoms)]
    coords = numpy.zeros((n_atoms, 3))
    for i in range(n_atoms):
        at_types[i] = xyz_array[i+2][0]
        for j in range(3):
            coords[i][j] = float(xyz_array[i+2][j+1])
    geom = [at_types, coords]
    return geom

# input syntax and usage warnings
def get_inputs():
    if (not len(sys.argv) == 2):
        print('\nUsage: %s XYZ_FILE\n' % sys.argv[0])
        print('  XYZ_FILE: coordinates of target molecule\n')
        sys.exit()
    else:
        xyz_file_name = sys.argv[1]
        xyz_name = xyz_file_name.split('/')[-1].split('.')[0]
        print('\ngeometry analysis for %s\n' % (xyz_name))
        return xyz_file_name

# print geometry to screen
def print_geom(geom, comment):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    print('%i\n%s\n' % (n_atoms, comment), end='')
    for i in range(n_atoms):
        print('%-2s' % (at_types[i]),end='')
        for j in range(3):
            print( ' %12.6f' % (coords[i][j]), end='')
        print('\n', end='')
    print('\n', end='')

# print bond graph to screen
def print_bond_graph(geom, bond_graph, comment):
    at_types = geom[0]
    n_atoms = len(at_types)
    print('%s\n' % (comment), end='')
    for i in range(n_atoms):
        print(' %4i %-2s -' % (i+1, at_types[i]), end='')
        for j in range(len(bond_graph[i])):
            print(' %i' % (bond_graph[i][j] + 1), end='')
        print('\n', end='')
    print('\n', end='')

# print list of bond lengths to screen
def print_bonds(geom, bonds):
    at_types = geom[0]
    n_bonds = len(bonds)
    print('%i bond length(s) found (Angstrom)' % (n_bonds))
    if (n_bonds > 0):
        print(' atoms            elements         values')
    for q in range(n_bonds):
        n1, n2  = bonds[q][0:2]
        r12 = bonds[q][2]
        nstr = '%i-%i' % (n1+1, n2+1)
        tstr = '(%s-%s) ' % (at_types[n1], at_types[n2])
        print(' %-15s  %-13s    %6.4f\n' % (nstr, tstr, r12), end='')
    print('\n', end='')
    
# print list of bond angles to screen
def print_angles(geom, angles):
    at_types = geom[0]
    n_angles = len(angles)
    print('%i bond angle(s) found (degrees)' % (n_angles))
    if (n_angles > 0):
        print(' atoms            elements         values')
    for q in range(n_angles):
        n1, n2, n3 = angles[q][0:3]
        a123 = angles[q][3]
        nstr = '%i-%i-%i' % (n1+1, n2+1, n3+1)
        tstr = '(%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3])
        print(' %-15s  %-13s   %7.3f\n' % (nstr, tstr, a123), end='')
    print('\n', end='')

# print list of torsion angles to screen
def print_torsions(geom, torsions):
    at_types = geom[0]
    n_torsions = len(torsions)
    print('%i torsion angle(s) found (degrees)' % (n_torsions))
    if (n_torsions > 0):
        print(' atoms            elements         values')
    for q in range(n_torsions):
        n1, n2, n3, n4 = torsions[q][0:4]
        t1234 = torsions[q][4]
        nstr = '%i-%i-%i-%i' % (n1+1, n2+1, n3+1, n4+1)
        tstr = '(%s-%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3], at_types[n4])
        print(' %-15s  %-13s  %8.3f\n' % (nstr, tstr, t1234), end='')
    print('\n', end='')

# print list of out-of-plane angles to screen
def print_outofplanes(geom, outofplanes):
    at_types = geom[0]
    n_outofplanes = len(outofplanes)
    print('%i out-of-plane angle(s) found (degrees)' % (n_outofplanes))
    if (n_outofplanes > 0):
        print(' atoms            elements         values')
    for q in range(n_outofplanes):
        n1, n2, n3, n4 = outofplanes[q][0:4]
        o1234 = outofplanes[q][4]
        nstr = '%i-%i-%i-%i' % (n1+1, n2+1, n3+1, n4+1)
        tstr = '(%s-%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3], at_types[n4])
        print(' %-15s  %-13s  %8.3f\n' % (nstr, tstr, o1234), end='')
    print('\n', end='')
    
# print center of mass coordinates
def print_com(com):
    print('molecular center of mass (Angstrom)\n   ', end='')
    print('       X             Y             Z\n  ', end='')
    for p in range(3):
        print(' %12.6f' % (com[p]), end='')
    print('\n')

# print moment of inertia tensor
def print_moi(moi):
    print('molecular moment of inertia tensor (amu * A^2)\n', end='')
    print('           X             Y             Z\n', end='')
    for p in range(3):
        print('%-2s' % (xyz[p]), end='')
        for q in range(3):
            print(' %12.6f' % (moi.item((p, q))), end='')
        print('\n', end='')
    print('\n', end='')
    
# print principal moments of inertia (eigenvalues of tensor)
def print_prinmom(prinmom):
    print('principal moments of inertia (amu * A^2)\n  ', end='')
    for p in range(3):
        print(' %12.6f' % (prinmom[p]), end='')
    print('\n')

# print molecule type
def print_moltype(moltype):
    print('this molecule is %s\n\n' % (moltype), end='')
    
# print rotational constants in mhz and wavenumber (cm^-1)
def print_rotfreq(freq, units):
    print('rotational frequencies (%s)\n  ' % (units), end='')
    for p in range(len(freq)):
        if (freq[p] < 10**(-1)):
            print(' %12.*e' % (3 - int(math.log(freq[p], 10.0)), freq[p]), end='')
        else:
            print(' %12.*f' % (6 - int(math.log(freq[p], 10.0)), freq[p]), end='')
    print('\n')
    
## MATH FUNCTIONS ##

# compare if two values are the same within specified threshold
def are_same(n1, n2, tol, minval):
    same = False
    nmax = max(abs(n1), abs(n2), abs(minval))
    comp = abs((n2 - n1) / nmax)
    if (comp <= 10**(-tol)):
        same = True
    return same

# calculate distance between two 3-d cartesian coordinates
def get_r12(coords1, coords2):
    r2 = 0.0
    for p in range(3):
        r2 += (coords2[p] - coords1[p])**2
    r = math.sqrt(r2)
    return r

# calculate unit vector between to 3-d cartesian coordinates
def get_u12(coords1, coords2):
    r12 = get_r12(coords1, coords2)
    u12 = [0.0 for p in range(3)]
    for p in range(3):
        u12[p] = (coords2[p] - coords1[p]) / r12
    return u12

# calculate dot product between two unit vectors
def get_udp(uvec1, uvec2):
    udp = 0.0
    for p in range(3):
        udp += uvec1[p] * uvec2[p]
    udp = max(min(udp, 1.0), -1.0)
    return udp

# calculate unit cross product between two unit vectors
def get_ucp(uvec1, uvec2):
    ucp = [0.0 for p in range(3)]
    cos_12 = get_udp(uvec1, uvec2)
    sin_12 = math.sqrt(1 - cos_12**2)
    ucp[0] = (uvec1[1]*uvec2[2] - uvec1[2]*uvec2[1]) / sin_12
    ucp[1] = (uvec1[2]*uvec2[0] - uvec1[0]*uvec2[2]) / sin_12
    ucp[2] = (uvec1[0]*uvec2[1] - uvec1[1]*uvec2[0]) / sin_12
    return ucp

# calculate angle between three 3-d cartesian coordinates
def get_a123(coords1, coords2, coords3):
    u21 = get_u12(coords2, coords1)
    u23 = get_u12(coords2, coords3)
    dp2123 = get_udp(u21, u23)
    a123 = rad2deg * math.acos(dp2123)
    return a123

# calculate torsion angle between four 3-d cartesian coordinates
def get_t1234(coords1, coords2, coords3, coords4):
    u21 = get_u12(coords2, coords1)
    u23 = get_u12(coords2, coords3)
    u32 = get_u12(coords3, coords2)
    u34 = get_u12(coords3, coords4)
    u21c23 = get_ucp(u21, u23)
    u32c34 = get_ucp(u32, u34)
    dp = get_udp(u21c23, u32c34)
    sign = 2 * float(get_udp(u21c23, u34) < 0) - 1
    t1234 = rad2deg * sign * math.acos(dp)
    return t1234

# calculate out-of-plane (improper torsion) angle between four 3-d cartesian coordinates
def get_o1234(coords1, coords2, coords3, coords4):
    u42 = get_u12(coords4, coords2)
    u43 = get_u12(coords4, coords3)
    u41 = get_u12(coords4, coords1)
    u42c43 = get_ucp(u42, u43)
    dp = get_udp(u42c43, u41)
    o1234 = rad2deg * math.asin(dp)
    return o1234

# translate coordinates by a defined vector and scale factor
def translate_coords(geom, vector, scale):
    coords = geom[1]
    n_atoms = len(coords)
    for i in range(n_atoms):
        for j in range(3):
            coords[i][j] += scale * vector[j]
    return geom

# calculate center of mass of a set of atoms
def get_com(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    com = [0.0 for p in range(3)]
    mass = 0.0
    for i in range(n_atoms):
        at_mass = at_masses[at_types[i]]
        mass += at_mass
        for j in range(3):
            com[j] += at_mass * coords[i][j]
    for p in range(3):
        com[p] /= mass
    return com

# calculate moment of inertia tensor for a set of atoms
def get_moi(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    moi = [[0.0 for q in range(3)] for p in range(3)]
    for i in range(n_atoms):
        at_mass = at_masses[at_types[i]]
        for p in range(3):
            for q in range(3):
                if (p == q):
                    r = (p+1) % 3
                    s = (p+2) % 3
                    moi[p][p] += at_mass * (coords[i][r]**2 + coords[i][s]**2)
                else:
                    moi[p][q] += -at_mass * coords[i][p] * coords[i][q]
    moi = numpy.matrix(moi)
    return moi

# calculate principal moments of inertia (eigenvalues of tensor)
def get_prinmom(moi):
    prinmom = numpy.linalg.eigvalsh(moi)
    return prinmom

# calculate rotational frequencies in MHz and wavenumbers (cm^-1)
def get_rotfreq(prinmom):
    freqcm1, freqmhz = [], []
    for p in range(3):
        iszero = are_same(prinmom[p], 0.0, mom_thresh, mom_min)
        degen = (p>0 and are_same(prinmom[p-1], prinmom[p], mom_thresh, mom_min))
        if (iszero or degen):
            continue
        mhz = h / (8 * math.pi**2 * prinmom[p])
        mhz *= (10**10)**2 * na * 10**3 * 10**(-6)
        cm1 = mhz / c * 10**6
        freqmhz.append(mhz)
        freqcm1.append(cm1)
    return freqmhz, freqcm1

# rotate molecule to inertial frame of principal moments
def get_inertial_coords(geom, moi):
    moi_eigvals, moi_eigvecs = numpy.linalg.eig(moi)
    coords = geom[1]
    coords = numpy.array(numpy.dot(coords, moi_eigvecs))
    geom[1] = coords
    moi = get_moi(geom)
    order = [0, 1, 2]
    for p in range(3):
        for q in range(p+1, 3):
            if (moi.item(p, p) < moi.item(q, q)):
                temp = order[p]
                order[p] = order[q]
                order[q] = temp
    moveaxes = numpy.zeros((3, 3))
    for p in range(3):
        moveaxes[p][order[p]] = 1.0
    coords = numpy.dot(coords, moveaxes)
    geom[1] = coords
    return geom

## TOPOLOGY FUNCTIONS ##

# build graph of which atoms are covalently bonded
def get_bond_graph(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    bond_graph = [[] for i in range(n_atoms)]
    for i in range(n_atoms):
        covrad1 = cov_rads[at_types[i]]
        for j in range(i+1, n_atoms):
            covrad2 = cov_rads[at_types[j]]
            thresh = bond_thresh * (covrad1 + covrad2)
            r12 = get_r12(coords[i], coords[j])
            if (r12 < thresh):
                bond_graph[i].append(j)
                bond_graph[j].append(i)
    return bond_graph

# determine atoms which are covalently bonded from bond graph
def get_bonds(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    bonds = []
    for i in range(n_atoms):      
        for a in range(len(bond_graph[i])):
            j = bond_graph[i][a]
            if (i < j):
                r12 = get_r12(coords[i], coords[j])
                bonds.append([i, j, r12])
    return bonds

# determine atoms which form a bond angle from bond graph
def get_angles(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    angles = []
    for j in range(n_atoms):
        n_jbonds = len(bond_graph[j])
        for a in range(n_jbonds):
            i = bond_graph[j][a]
            for b in range(a+1, n_jbonds):
                k = bond_graph[j][b]
                a123 = get_a123(coords[i], coords[j], coords[k])
                angles.append([i, j, k, a123])
    angles = sorted(angles, key=lambda angle:angle[0])
    return angles

# determine atoms which form torsion angles from bond graph
def get_torsions(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    torsions = []
    for j in range(n_atoms):
        n_jbonds = len(bond_graph[j])
        for a in range(n_jbonds):
            k = bond_graph[j][a]
            if (k < j):
                continue
            n_kbonds = len(bond_graph[k])
            for b in range(n_jbonds):
                i = bond_graph[j][b]
                if (i == k):
                    continue
                for c in range(n_kbonds):
                    l = bond_graph[k][c]
                    if (l == j or l == i):
                        continue
                    t1234 = get_t1234(coords[i], coords[j], coords[k], coords[l])
                    torsions.append([i, j, k, l, t1234])
    torsions = sorted(torsions, key=lambda torsion:torsion[0])
    return torsions

# determine atoms which form out-of-plane angles from bond graph
def get_outofplanes(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    outofplanes = []
    for l in range(n_atoms):
        n_lbonds = len(bond_graph[l])
        for a in range(n_lbonds):
            i = bond_graph[l][a]
            for b in range(n_lbonds):
                j = bond_graph[l][b]
                if (i == j):
                    continue
                for c in range(b+1, n_lbonds):
                    k = bond_graph[l][c]
                    if (i == k):
                        continue
                    o1234 = get_o1234(coords[i], coords[j], coords[k], coords[l])
                    outofplanes.append([i, j, k, l, o1234])
    outofplanes = sorted(outofplanes, key=lambda outofplane:outofplane[0])
    return outofplanes

# determine molecule type based on principal moments of inertia
def get_moltype(geom, pm):
    same12  = are_same(pm[0], pm[1], mom_thresh, mom_min)
    same13  = are_same(pm[0], pm[2], mom_thresh, mom_min)
    same23  = are_same(pm[1], pm[2], mom_thresh, mom_min)
    onezero = are_same(pm[0], 0.0,   mom_thresh, mom_min)
    allzero = are_same(pm[2], 0.0,   mom_thresh, mom_min)
    if (allzero):
        moltype = 'monatomic'
    elif (onezero):
        moltype = 'linear'
    elif (same13):
        moltype = 'a spherical top'
    elif (same12 or same23):
        moltype = 'a symmetric top'
    else:
        moltype = 'an asymmetric top'
    return moltype

## MAIN BLOCK ##

# read in geometry, determine bonded topology
xyz_file_name = get_inputs()
geom = get_geom(xyz_file_name)
bond_graph = get_bond_graph(geom)

# calculate bond lengths, angles, torsions, outofplanes, and center of mass
bonds = get_bonds(geom, bond_graph)
angles = get_angles(geom, bond_graph)
torsions = get_torsions(geom, bond_graph)
outofplanes = get_outofplanes(geom, bond_graph)
com = get_com(geom)

# print resulting values
print_geom(geom, 'initial geometry')
print_bonds(geom, bonds)
print_angles(geom, angles)
print_torsions(geom, torsions)
print_outofplanes(geom, outofplanes)
print_com(com)

# translate molecule to center of mass and print
com_geom = translate_coords(geom, com, -1.0)
print_geom(com_geom, 'center of mass translated geometry')

# calculate moment of inertia tensor, principal moments, molecule type, and rotational frequencies
moi = get_moi(com_geom)
prinmom = get_prinmom(moi)
moltype = get_moltype(com_geom, prinmom)
rotmhz, rotcm1 = get_rotfreq(prinmom)

# print resulting values
print_moi(moi)
print_prinmom(prinmom)
print_moltype(moltype)
print_rotfreq(rotmhz, 'MHz')
print_rotfreq(rotcm1, 'cm^-1')

# rotate molecule to inertial frame and print
moi_geom = get_inertial_coords(com_geom, moi)
print_geom(moi_geom, 'principal moment aligned geometry')

print('geometry analysis completed successfully\n')
# end of program

