import os, sys, math
import numpy as np
from mmlib import param, geomcalc, topology, energy, molecule

# fileio.py: functions for reading in and printing out molecular mechanics data

# create 2d array of strings from input file name
def get_file_string_array(infile_name):
    try:
        infile = open(infile_name, 'r')
    except IOError:
        print('Error: file (%s) does not exist!' % (infile_name))
        sys.exit()
    infile_data = infile.readlines()
    infile.close()
    infile_array = []
    for line in infile_data:
        infile_array.append(line.split())
    return infile_array

# get symbols, coordinates, and parameters from xyzq file
def get_geom(mol):
    infile_array = get_file_string_array(mol.infile)
    mol.n_atoms = int(infile_array[0][0])
    for i in range(mol.n_atoms):
        at_type = infile_array[i+2][0]
        at_coords = np.zeros(3)
        for j in range(3):
            at_coords[j] = float(infile_array[i+2][j+1])
        at_charge = float(infile_array[i+2][4])
        at_element = at_type[0].capitalize()
        if (at_type[-1].islower()): at_element += at_type[-1]
        at_mass = param.get_at_mass(at_element)
        at_ro, at_eps = param.get_vdw_param(at_type)
        new_atom = molecule.atom(at_type, at_coords, at_charge, at_ro,
            at_eps, at_mass)
        new_atom.set_covrad(param.get_cov_rad(at_element))
        mol.atoms.append(new_atom)
        mol.mass += at_mass

# parse atom record into an atom object
def get_atom(record):
    at_type = record[2]
    at_coords = np.zeros(3)
    for j in range(3):
        at_coords[j] = float(record[j+3])
    at_charge, at_ro, at_eps = (float(record[6]), float(record[7]),
        float(record[8]))
    at_element = at_type[0].capitalize()
    if (at_type[-1].islower()):
        at_element += at_type[-1]
    at_mass = param.get_at_mass(at_element)
    atom = molecule.atom(at_type, at_coords, at_charge, at_ro,
        at_eps, at_mass)
    atom.set_covrad(param.get_cov_rad(at_element))
    return atom

# parse bond record into a bond object
def get_bond(atoms, record):
    at1, at2 = int(record[1])-1, int(record[2])-1
    k_b, r_eq = float(record[3]), float(record[4])
    c1, c2 = atoms[at1].coords, atoms[at2].coords
    r_ij = geomcalc.get_r_ij(c1, c2)
    bond = molecule.bond(at1, at2, r_ij, r_eq, k_b)
    return bond

# parse angle record into an angle object
def get_angle(atoms, record):
    at1, at2, at3 = int(record[1])-1, int(record[2])-1, int(record[3])-1
    k_a, a_eq = float(record[4]), float(record[5])
    c1, c2, c3 = (atoms[at1].coords, atoms[at2].coords,
        atoms[at3].coords)
    a_ijk = geomcalc.get_a_ijk(c1, c2, c3)
    angle = molecule.angle(at1, at2, at3, a_ijk, a_eq, k_a)
    return angle

# parse torsion record into a torsion object
def get_torsion(atoms, record):
    at1, at2, at3, at4 = (int(record[1])-1, int(record[2])-1,
        int(record[3])-1, int(record[4])-1)
    v_n, gamma, nfold, paths = (float(record[5]), float(record[6]),
        int(record[7]), int(record[8]))
    c1, c2, c3, c4 = (atoms[at1].coords, atoms[at2].coords,
        atoms[at3].coords, atoms[at4].coords)
    t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4)
    torsion = molecule.torsion(at1, at2, at3, at4, t_ijkl,
        v_n, gamma, nfold, paths)
    return torsion

# parse outofplane record into an outofplane object
def get_outofplane(atoms, record):
    at1, at2, at3, at4 = (int(record[1])-1, int(record[2])-1,
        int(record[3])-1, int(record[4])-1)
    v_n, gamma, nfold = float(record[5]), float(record[6]), int(record[7])
    c1, c2, c3, c4 = (atoms[at1].coords, atoms[at2].coords,
        atoms[at3].coords, atoms[at4].coords)
    o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4) 
    outofplane = molecule.outofplane(at1, at2, at3, at4, o_ijkl,
        v_n, gamma, nfold)
    return outofplane

# get symbols, coordinates, parameters, and topology from prm file
def get_prm(mol):
    infile_array = get_file_string_array(mol.infile)
    for i in range(len(infile_array)):
        record = infile_array[i]
        rec_type = record[0]
        if (rec_type == 'ATOM'):
            mol.atoms.append(get_atom(record))
            mol.mass += mol.atoms[-1].mass
    mol.n_atoms = len(mol.atoms)
    for i in range(len(infile_array)):
        record = infile_array[i]
        rec_type = record[0]
        if (rec_type == 'BOND'):
            mol.bonds.append(get_bond(mol.atoms, record))
        elif (rec_type == 'ANGLE'):
            mol.angles.append(get_angle(mol.atoms, record))
        elif (rec_type == 'TORSION'):
            mol.torsions.append(get_torsion(mol.atoms, record))
        elif (rec_type == 'OUTOFPLANE'):
            mol.outofplanes.append(get_outofplane(mol.atoms, record))
    mol.n_bonds = len(mol.bonds)
    mol.n_angles = len(mol.angles)
    mol.n_torsions = len(mol.torsions)
    mol.n_outofplanes = len(mol.outofplanes)
    topology.get_nonints(mol)

# get simulation data from file
def get_sim_data(sim):
    infile_array = get_file_string_array(sim.infile)
    for q in range(len(infile_array)):
        if (len(infile_array[q]) < 2): continue
        kwarg = infile_array[q][0]
        kwargval = infile_array[q][1]
        kwargarr = infile_array[q][1:]
        if (kwarg == 'MOLECULE'):
            sim.mol = molecule.molecule(sim.indir + '/' + kwargval)
        elif (kwarg == 'TEMPERATURE'):
            sim.temp = float(kwargval)
        elif (kwarg == 'PRESSURE'):
            sim.press = float(kwargval)
        elif (kwarg == 'BOUNDARYSPRING'):
            sim.mol.k_box = float(kwargval)
        elif (kwarg == 'BOUNDARY'):
            sim.mol.bound = float(kwargval)
            sim.mol.get_volume()
        elif (kwarg == 'BOUNDARYTYPE'):
            sim.mol.boundtype = kwargval.lower()
            sim.mol.get_volume()
        elif (kwarg == 'ORIGIN'):
            sim.mol.origin = [float(kwargarr[i]) for i in range(3)]
        elif (kwarg == 'TOTALTIME'):
            sim.tottime = float(kwargval)
        elif (kwarg == 'TOTALCONFS'):
            sim.totconfs = int(kwargval)
        elif (kwarg == 'TIMESTEP'):
            sim.timestep = float(kwargval)
        elif (kwarg == 'GEOMTIME'):
            sim.geomtime = float(kwargval)
        elif (kwarg == 'GEOMCONF'):
            sim.geomconf = int(kwargval)
        elif (kwarg == 'GEOMOUT'):
            sim.geomout = sim.indir + '/' + kwargval
        elif (kwarg == 'ENERGYTIME'):
            sim.energytime = float(kwargval)
        elif (kwarg == 'ENERGYCONF'):
            sim.energyconf = int(kwargval)
        elif (kwarg == 'ENERGYOUT'):
            sim.energyout = sim.indir + '/' + kwargval
        elif (kwarg == 'STATUSTIME'):
            sim.statustime = float(kwargval)
        elif (kwarg == 'EQTIME'):
            sim.eqtime = float(kwargval)
        elif (kwarg == 'EQRATE'):
            sim.eqrate = float(kwargval)

# print atomic coordinates for a set of atoms
def print_coords(mol, comment):
    print('%i\n%s\n' % (mol.n_atoms, comment), end='')
    for i in range(mol.n_atoms):
        print('%-2s' % (mol.atoms[i].attype), end='')
        for j in range(3):
            print(' %12.6f' % (mol.atoms[i].coords[j]), end='')
        print('\n', end='')
    print('\n', end='')

# print atomic energy gradient for a set of atoms
def print_gradient(mol, grad_type):
    if   (grad_type == 'total'): grad = mol.g_total
    elif (grad_type == 'nonbonded'): grad = mol.g_nonbonded
    elif (grad_type == 'bonded'): grad = mol.g_bonded
    elif (grad_type == 'boundary'): grad = mol.g_bound
    elif (grad_type == 'vdw'): grad = mol.g_vdw
    elif (grad_type == 'elst'): grad = mol.g_elst
    elif (grad_type == 'bonds'): grad = mol.g_bonds
    elif (grad_type == 'angles'): grad = mol.g_angles
    elif (grad_type == 'torsions'): grad = mol.g_torsions
    elif (grad_type == 'outofplanes'): grad = mol.g_outofplanes
    else: print('Error: grad type (%s) not recognized!' % (grad_type))
    print('\n %s\n\n' % ('%s %s gradient' % (mol.grad_type,
        grad_type)), end='')
    for i in range(mol.n_atoms):
        print('%-2s' % (mol.atoms[i].attype), end='')
        for j in range(3):
            print(' %12.6f' % (grad[i][j]), end='')
        print('\n', end='')
    print('\n', end='')

# print a string in the middle of dash banner of specified length
def print_banner(string, length, newline1, newline2):
    n_dash1 = math.floor((length - len(string))/2) - 1
    n_dash2 = math.ceil((length - len(string))/2) - 1
    for i in range(newline1):
        print('')
    print_string = ' '
    for i in range(n_dash1):
        print_string += '-'
    print_string += string
    for i in range(n_dash2):
        print_string += '-'
    print(print_string, end='')
    for i in range(newline2):
        print('')

# print an array of strings padded by specified number of spaces
def print_padded(strings, spacings):
    print_string = ''
    for i in range(len(strings)):
        print_string += '%*s%s' % (spacings[i], '', strings[i])
    print(print_string)

# print header for a section of output
def print_header(header, n_banner, params, spaces):
    print_banner(header, n_banner, 1, 1)
    print_padded(params, spaces)
    print_banner('', n_banner, 0, 1)

# print geometry and non-bonded parameters for a set of atoms
def print_geom(mol):
    header = ' Molecular Geometry and Non-bonded Parameters '
    params = ['type', 'x', 'y', 'z', 'q', 'ro/2', 'eps']
    spaces = [6, 6, 11, 11, 10, 6, 4]
    n_banner = 73
    print_header(header, n_banner, params, spaces)
    for i in range(mol.n_atoms):
        print('%4i | %-2s' % (i+1, mol.atoms[i].attype), end='')
        for j in range(3):
            print('%12.6f' % (mol.atoms[i].coords[j]), end='')
        print('  %8.5f %7.4f %7.4f' % (mol.atoms[i].charge,
            mol.atoms[i].ro, mol.atoms[i].eps))

# print geometry and non-bonded parameters for a set of atoms
def print_geom_file(outfile, mol):
    outfile.write('# %s Atoms (at, type, x, y, z, q, ro, eps)\n' % (
        mol.n_atoms))
    for i in range(mol.n_atoms):
        outfile.write('ATOM %4i %-2s' % (i+1, mol.atoms[i].attype))
        for j in range(3):
            outfile.write(' %11.6f' % (mol.atoms[i].coords[j]))
        outfile.write(' %8.5f %7.4f %7.4f\n' % (mol.atoms[i].charge,
            mol.atoms[i].ro, mol.atoms[i].eps))

# print bond topology and bond parameters from an array
def print_bonds(mol):
    if (mol.n_bonds > 0):
        header, n_banner = ' Bond Length Data ', 57
        params = ['k_b', 'r_eq', 'r_ij', 'types', 'energy', 'atoms']
        spaces = [10, 5, 5, 3, 4, 1]
        print_header(header, n_banner, params, spaces)
    else:
        print('\n No Bonds Detected')
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        t1, t2 = mol.atoms[bond.at1].attype, mol.atoms[bond.at2].attype
        pstr = '%4i | %7.2f %8.4f %8.4f (%2s-%2s) %8.4f (%i-%i)' % (p+1,
            bond.k_b, bond.r_eq, bond.r_ij, t1, t2, bond.e, bond.at1+1,
            bond.at2+1)
        print(pstr)

# print bond topology and bond parameters to parameter file
def print_bonds_file(outfile, mol):
    outfile.write('# %i Bonds (At1, At2, K_b, R_eq)\n' % (mol.n_bonds))
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        outfile.write('BOND %4i %4i %7.2f %7.4f\n' % (bond.at1+1, bond.at2+1,
            bond.k_b, bond.r_eq))

# print bond angle topology and angle parameters from an array
def print_angles(mol):
    if (mol.n_angles > 0):
        header, n_banner = ' Bond Angle Data ', 58
        params = ['k_a', 'a_eq', 'a_ijk', 'types', 'energy', 'atoms']
        spaces = [9, 3, 4, 4, 5, 2]
        print_header(header, n_banner, params, spaces)
    else:
        print('\n No Bond Angles Detected')
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        t1, t2, t3 = (mol.atoms[ang.at1].attype, mol.atoms[ang.at2].attype,
            mol.atoms[ang.at3].attype)
        pstr = '%4i | %6.2f %7.3f %7.3f (%2s-%2s-%2s) %7.4f (%i-%i-%i)' % (
            p+1, ang.k_a, ang.a_eq, ang.a_ijk, t1, t2, t3, ang.e, ang.at1+1,
            ang.at2+1, ang.at3+1)
        print(pstr)

# print bond angle topology and angle parameters to parameter file
def print_angles_file(outfile, mol):
    outfile.write('# %i Angles (At1, At2, At3, K_a, A_eq)\n' % (mol.n_angles))
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        outfile.write('ANGLE %4i %4i %4i %7.4f %8.4f\n' % (ang.at1+1,
            ang.at2+1, ang.at3+1, ang.k_a, ang.a_eq))

# print torsion topology and torsion parameters from an array
def print_torsions(mol):
    if (mol.n_torsions > 0):
        header, n_banner = ' Torsion Angle Data ', 67
        params = ['vn/2', 'gamma', 't_ijkl n p', 'types', 'energy', 'atoms']
        spaces = [9, 2, 3, 5, 6, 3]
        print_header(header, n_banner, params, spaces)
    else:
        print('\n No Torsion Angles Detected')
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        t1, t2 = mol.atoms[tor.at1].attype, mol.atoms[tor.at2].attype
        t3, t4 = mol.atoms[tor.at3].attype, mol.atoms[tor.at4].attype
        pstr = '%4i | %6.2f %6.1f %8.3f %i %i (%2s-%2s-%2s-%2s)' % (p+1,
            tor.v_n, tor.gam, tor.t_ijkl, tor.n, tor.paths, t1, t2, t3, t4)
        pstr += ' %7.4f (%i-%i-%i-%i)' % (tor.e, tor.at1+1, tor.at2+2,
            tor.at3+1, tor.at4+1)
        print(pstr)

# print torsion topology and torsion parameters from an array
def print_torsions_file(outfile, mol):
    outfile.write('# %i Torsions (At1, At2, At3, At4,' % (mol.n_torsions))
    outfile.write('V_n, Gamma, N_f, paths)\n')
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        outfile.write('TORSION %4i %4i %4i %4i %6.2f %6.1f %i %i\n' % (
            tor.at1+1, tor.at2+1, tor.at3+1, tor.at4+1, tor.v_n, tor.gam,
            tor.n, tor.paths))

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes(mol):
    if (mol.n_outofplanes > 0):
        header, n_banner = ' Out-of-plane Angle Data ', 55
        params = ['vn/2', 'o_ijkl', 'types', 'energy', 'atoms']
        spaces = [8, 3, 5, 6, 3]
        print_header(header, n_banner, params, spaces)
    else:
        print('\n No Out-of-plane Angles Detected')
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        t1, t2 = mol.atoms[oop.at1].attype, mol.atoms[oop.at2].attype
        t3, t4 = mol.atoms[oop.at3].attype, mol.atoms[oop.at4].attype
        pstr = '%4i | %6.2f %7.3f (%2s-%2s-%2s-%2s) %7.4f (%i-%i-%i-%i)' % (
            p+1, oop.v_n, oop.o_ijkl, t1, t2, t3, t4, oop.e, oop.at1+1,
            oop.at2+1, oop.at3+1, oop.at4+1)
        print(pstr)

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes_file(outfile, mol):
    outfile.write('# %i Outofplanes (At1, At2, At3, At4,' % (mol.n_outofplanes))
    outfile.write(' V_n, Gamma, N_f)\n')
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        outfile.write('OUTOFPLANE %4i %4i %4i %4i %6.2f %6.1f %i\n' % (
            oop.at1+1, oop.at2+1, oop.at3+1, oop.at4+1, oop.v_n, oop.gam,
            oop.n_fold))

# print energy values to screen
def print_energy(mol):
    print_banner(' Energy Values ', 42, 1, 1)
    print('%10.4f kcal/mol Total Energy' % (mol.e_total))
    print('%10.4f kcal/mol Kinetic Energy' % (mol.e_kinetic))
    print('%10.4f kcal/mol Potential Energy' % (mol.e_potential))
    print('%10.4f kcal/mol Non-bonded Energy' % (mol.e_nonbonded))
    print('%10.4f kcal/mol Bonded Energy' % (mol.e_bonded))
    print('%10.4f kcal/mol Boundary Energy' % (mol.e_bound))
    print('%10.4f kcal/mol van der Waals Energy' % (mol.e_vdw))
    print('%10.4f kcal/mol Electrostatic Energy' % (mol.e_elst))
    print('%10.4f kcal/mol Bond Energy' % (mol.e_bonds))
    print('%10.4f kcal/mol Angle Energy' % (mol.e_angles))
    print('%10.4f kcal/mol Torsion Energy' % (mol.e_torsions))
    print('%10.4f kcal/mol Out-of-plane Energy' % (mol.e_outofplanes))

# check for proper input arguments and return result
def get_input():
    prog_name = sys.argv[0].split('/')[-1]
    if (not len(sys.argv) >= 2):
        print('\nUsage: %s INPUT_FILE\n' % (prog_name))
        print('  INPUT_FILE: ', end='')
        if (prog_name == 'mm.py'):
            print('xyzq or prm file for molecular mechanics\n')
        elif (prog_name == 'md.py'):
            print('simulation file for molecular dynamics\n')
        elif (prog_name == 'mmc.py'):
            print('simulation file for metropolis monte carlo\n')
        sys.exit()
    else:
        infile_name = sys.argv[1]
    return infile_name

# end of module

