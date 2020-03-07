import sys, math
import numpy as np
import numpy.linalg as la

## CONSTANTS ##

# conversion from radians to degrees and vice versa
rad2deg = 180.0 / math.pi
deg2rad = math.pi / 180.0

# cartesian indices
xyz = {0: 'X', 1: 'Y', 2: 'Z', 'X': 0, 'Y': 1, 'Z': 2}

# relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_masses = {    'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067,
  'F' : 18.9984, 'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
  'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
  'Be': 9.01218, 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
  'Si': 28.0855, 'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
  'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
  'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
  'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X' :  0.0000}

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
    
# input syntax and usage warnings
def get_input():
    if (not len(sys.argv) == 2):
        print('\nUsage: %s ZMAT_FILE\n' % sys.argv[0])
        print('  ZMAT_FILE: z-matrix file of target molecule\n')
        sys.exit()
    else:
        zmat_file_name = sys.argv[1]
        return zmat_file_name

## MATH FUNCTIONS ##

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

# get local axis system from 3 coordinates
def get_local_axes(coords1, coords2, coords3):
    u21 = get_u12(coords1, coords2)
    u23 = get_u12(coords2, coords3)
    if (abs(get_udp(u21, u23)) >= 1.0):
      print('\nError: Co-linear atoms in an internal coordinate definition')
      sys.exit()
    u23c21 = get_ucp(u23, u21)
    u21c23c21 = get_ucp(u21, u23c21)
    z = u21
    y = u21c23c21
    x = get_ucp(y, z)
    local_axes = [x, y, z]
    return local_axes

# calculate vector of bond in local axes of internal coordinates
def get_bond_vector(r, a, t):
    x = r * math.sin(a) * math.sin(t)
    y = r * math.sin(a) * math.cos(t)
    z = r * math.cos(a)
    bond_vector = [x, y, z]
    return bond_vector

## CLASSES ##

# atom class for atomic data
class atom:
    # constructor
    def __init__(self, at_type, rnum, anum, tnum, rval, aval, tval):
        self.attype = at_type
        self.rnum = rnum
        self.anum = anum
        self.tnum = tnum
        self.rval = rval
        self.aval = aval
        self.tval = tval
        self.coords = [None for j in range(3)]
        self.mass = at_masses[self.attype]

# molecule class for molecular data
class molecule:
    # constructor
    def __init__(self, zmat_file_name):
        self.zmat_file = zmat_file_name
        self.name = self.zmat_file.split('/')[-1].split('.')[0]
        self.get_zmat(self.zmat_file)

    # read in z-matrix from zmat file
    def get_zmat(self, zmat_file_name):
        zmat_array = get_file_string_array(zmat_file_name)
        self.n_atoms = int(zmat_array[0][0])
        self.comment = ' '.join(zmat_array[1])
        self.atoms = []
        for i in range(self.n_atoms):
            at_type = zmat_array[i+2][0]
            if (i >= 1):
                rnum = int(zmat_array[i+2][1]) - 1
                rval = float(zmat_array[i+2][2])
            else:
                rnum, rval = None, None
            if (i >= 2):
                anum = int(zmat_array[i+2][3]) - 1
                aval = deg2rad * float(zmat_array[i+2][4])
            else:
                anum, aval = None, None
            if (i >= 3):
                tnum = int(zmat_array[i+2][5]) - 1
                tval = deg2rad * float(zmat_array[i+2][6])
            else:
                tnum, tval = None, None
            self.atoms.append(atom(at_type, rnum, anum, tnum, rval, aval, tval))

    # print z-matrix to screen
    def print_zmat(self, comment):
        print('%i\n%s\n' % (self.n_atoms, comment), end='')
        for i in range(self.n_atoms):
            print('%-2s' % (self.atoms[i].attype), end='')
            if (i >= 1):
                print(' %3i' % (self.atoms[i].rnum+1), end='')
                print('%8.4f' % (self.atoms[i].rval), end='')
            if (i >= 2):
                print(' %3i' % (self.atoms[i].anum+1), end='')
                print('%8.3f' % (rad2deg * self.atoms[i].aval), end='')
            if (i >= 3):
                print(' %3i' % (self.atoms[i].tnum+1), end='')
                print('%8.3f' % (rad2deg * self.atoms[i].tval), end='')
            print('\n', end='')

    # print xyz coordinates to screen
    def print_coords(self, comment):
        print('%i\n%s\n' % (self.n_atoms, comment), end='')
        for i in range(self.n_atoms):
            print('%-2s' % (self.atoms[i].attype), end='')
            for j in range(3):
                print(' %12.6f' % (self.atoms[i].coords[j]), end='')
            print('\n', end='')

    # obtain cartesian xyz-coordinates from z-matrix values
    def zmat2xyz(self):
        if (self.n_atoms >= 1):
            self.atoms[0].coords = [0.0, 0.0, 0.0]
        if (self.n_atoms >= 2):
            self.atoms[1].coords = [0.0, 0.0, self.atoms[1].rval]
        if (self.n_atoms >= 3):
            r1,  r2  = self.atoms[1].rval, self.atoms[2].rval
            rn1, rn2 = self.atoms[1].rnum, self.atoms[2].rnum
            a1 = self.atoms[2].aval
            y = r2*math.sin(a1)
            z = self.atoms[rn2].coords[2] + (1-2*float(rn2==1))*r2*math.cos(a1)
            self.atoms[2].coords = [0.0, y, z]
        for i in range(3, self.n_atoms):
            atom = self.atoms[i]
            coords1 = self.atoms[atom.rnum].coords
            coords2 = self.atoms[atom.anum].coords
            coords3 = self.atoms[atom.tnum].coords
            self.atoms[i].local_axes = get_local_axes(coords1, coords2, coords3)
            bond_vector = get_bond_vector(atom.rval, atom.aval, atom.tval)
            disp_vector = np.array(np.dot(bond_vector, self.atoms[i].local_axes))
            for p in range(3):
                atom.coords[p] = self.atoms[atom.rnum].coords[p] + disp_vector[p]

## MAIN BLOCK ##

# get input arguments
zmat_file = get_input()

# read in z-matrix
mol = molecule(zmat_file)

# convert to xyz-coordinates
mol.zmat2xyz()

# print results
mol.print_coords(mol.comment)

# end of program

