import os, sys
import param, geomcalc, topology, energy, molecule
import numpy as np

# fileio.py: functions for reading in and printing out molecular mechanics data

# list of modules imported from mmlib
modules = ['__init__', 'fileio', 'param', 'geomcalc', 'topology', 'energy', 'molecule']

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

# get symbols, coordinates, and parameters from file
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
        new_atom = molecule.atom(at_type, at_coords, at_charge, at_ro, at_eps, at_mass)
        new_atom.set_covrad(param.get_cov_rad(at_element))
        mol.atoms.append(new_atom)

# get simulation data from file
def get_sim_data(sim):
    infile_array = get_file_string_array(sim.infile)
    for q in range(len(infile_array)):
        if (len(infile_array[q]) < 2): continue
        kwarg = infile_array[q][0]
        kwargval = infile_array[q][1]
        if (kwarg == 'MOLECULE'):
            sim.mol = molecule.molecule(kwargval)
        elif (kwarg == 'TEMPERATURE'):
            sim.temp = float(kwargval)
        elif (kwarg == 'TOTALTIME'):
            sim.tottime = float(kwargval)
        elif (kwarg == 'TIMESTEP'):
            sim.timestep = float(kwargval)
        elif (kwarg == 'GEOMTIME'):
            sim.geomtime = float(kwargval)
        elif (kwarg == 'GEOMOUT'):
            sim.geomout = kwargval
        elif (kwarg == 'ENERGYTIME'):
            sim.energytime = float(kwargval)
        elif (kwarg == 'ENERGYOUT'):
            sim.energyout = kwargval

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
    elif (grad_type == 'vdw'): grad = mol.g_vdw
    elif (grad_type == 'elst'): grad = mol.g_elst
    elif (grad_type == 'bonds'): grad = mol.g_bonds
    elif (grad_type == 'angles'): grad = mol.g_angles
    elif (grad_type == 'torsions'): grad = mol.g_torsions
    elif (grad_type == 'outofplanes'): grad = mol.g_outofplanes
    else: print 'Error: grad type (%s) not recognized!' % (grad_type))
    print('\n %s\n\n' % ('%s %s gradient' % (mol.grad_type, grad_type)), end='')
    for i in range(mol.n_atoms):
        print('%-2s' % (mol.atoms[i].attype), end='')
        for j in range(3):
            print(' %12.6f' % (grad[i][j]), end='')
        print('\n', end='')
    print('\n', end='')

# print geometry and non-bonded parameters for a set of atoms
def print_geom(mol):
    print('\n ------------- Molecular Geoemetry and Non-bonded Parameters --------------')
    print('     | Atom      X            Y            Z', end='')
    print('          Q       Ro/2    Eps')
    print(    ' --------------------------------------------------------------------------')
    for i in range(mol.n_atoms):
        print('%4i | %-2s' % (i+1, mol.atoms[i].attype), end='')
        for j in range(3):
            print('%12.6f' % (mol.atoms[i].coords[j]), end='')
        print('  %8.5f %7.4f %7.4f' % (mol.atoms[i].charge, mol.atoms[i].ro, mol.atoms[i].eps))

# print geometry and non-bonded parameters for a set of atoms
def print_geom_file(outfile, mol):
    outfile.write('# %s Atoms (At, Type, x, y, z, q, Ro, Eps)\n' % (mol.n_atoms))
    for i in range(mol.n_atoms):
        outfile.write('ATOM %4i %-2s' % (i+1, mol.atoms[i].attype))
        for j in range(3):
            outfile.write(' %11.6f' % (mol.atoms[i].coords[j]))
        outfile.write(' %8.5f %7.4f %7.4f\n' % (mol.atoms[i].charge, mol.atoms[i].ro, mol.atoms[i].eps))

# print bond topology and bond parameters from an array
def print_bonds(mol):
    if (mol.n_bonds > 0):
        print('\n --------------------- Bond Length Data ----------------------')
        print('     |    k_b      r_eq      r_ij     types   energy   atoms')
        print(  ' -----------------------------------------------------------)--'
    else:
        print('\n No Bonds Detected')
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        type1, type2 = mol.atoms[bond.at1].attype, mol.atoms[bond.at2].attype
        print('%4i | %7.2f  %8.4f  %8.4f  ' % (p+1, bond.k_b, bond.r_eq, bond.r_ij), end='')
        print('(%-2s-%-2s) %7.4f' % (type1, type2, bond.e), end='')
        print(' (%i-%i)' % (bond.at1+1, bond.at2+1))

# print bond topology and bond parameters to parameter file
def print_bonds_file(outfile, mol):
    outfile.write('# %i Bonds (At1, At2, K_b, R_eq)\n' % (mol.n_bonds))
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        outfile.write('BOND %4i %4i' % (bond.at1+1, bond.at2+1))
        outfile.write(' %7.2f %7.4f\n' % (bond.k_b, bond.r_eq))

# print bond angle topology and angle parameters from an array
def print_angles(mol):
    if (mol.n_angles > 0):
        print('\n ------------------- Bond Angle Data -----------------------')
        print('     |   k_a    a_eq     a_ijk    types     energy   atoms')
        print(  ' -----------------------------------------------------------')
    else:
        print('\n No Bond Angles Detected')
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        type1 = mol.atoms[ang.at1].attype
        type2 = mol.atoms[ang.at2].attype
        type3 = mol.atoms[ang.at3].attype
        print('%4i | %6.2f  %7.3f  %7.3f' % (p+1, ang.k_a, ang.a_eq, ang.a_ijk), end='')
        print('(%-2s-%-2s-%-2s)' % (type1, type2, type3), end='')
        print('%7.4f  (%i-%i-%i)' % (ang.e, ang.at1+1, ang.at2+1, ang.at3+1))

# print bond angle topology and angle parameters to parameter file
def print_angles_file(outfile, mol):
    outfile.write('# %i Angles (At1, At2, At3, K_a, A_eq)\n' % (mol.n_angles))
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        outfile.write('ANGLE %4i %4i %4i' % (ang.at1+1, ang.at2+1, ang.at3+1))
        outfile.write(' %7.4f %8.4f\n' % (ang.k_a, ang.a_eq))

# print torsion topology and torsion parameters from an array
def print_torsions(mol):
    if (mol.n_torsions > 0):
        print('\n ------------------- Torsion Angle Data -------------------------')
        print('     |   Vn/2  gamma   t_ijkl n p     types      energy   atoms')
        print(  ' ----------------------------------------------------------------')
    else:
        print('\n No Torsion Angles Detected')
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        type1, type2 = mol.atoms[tor.at1].attype, mol.atoms[tor.at2].attype
        type3, type4 = mol.atoms[tor.at3].attype, mol.atoms[tor.at4].attype
        print('%4i | %6.2f %6.1f' % (p+1, tor.v_n, tor.gam), end='')
        print('%8.3f %1i %1i' % (tor.t_ijkl, tor.n, tor.paths), end='')
        print('(%-2s-%-2s-%-2s-%-2s)' % (type1, type2, type3, type4), end='')
        print('%7.4f'  % (tor.e), end='')
        print(' (%i-%i-%i-%i)' % (tor.at1+1, tor.at2+1, tor.at3+1, tor.at4+1))

# print torsion topology and torsion parameters from an array
def print_torsions_file(outfile, mol):
    outfile.write('# %i Torsions (At1, At2, At3, At4,' % (mol.n_torsions))
    outfile.write('V_n, Gamma, N_f, paths)\n')
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        outfile.write('TORSION %4i %4i %4i %4i' % (tor.at1+1, tor.at2+1, tor.at3+1, tor.at4+1))
        outfile.write(' %6.2f %6.1f %1i %1i\n' % (tor.v_n, tor.gam, tor.n, tor.paths))

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes(mol):
    if (mol.n_outofplanes > 0):
        print('\n ---------------- Out-of-plane Angle Data ----------------')
        print('     |  Vn/2   o_ijkl      types      energy   atoms')
        print(  ' ---------------------------------------------------------')
    else:
        print('\n No Out-of-plane Angles Detected')
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        type1, type2 = mol.atoms[oop.at1].attype, mol.atoms[oop.at2].attype
        type3, type4 = mol.atoms[oop.at3].attype, mol.atoms[oop.at4].attype
        print('%4i | %6.2f %7.3f' % (p+1, oop.v_n, oop.o_ijkl), end='')
        print('  (%-2s-%-2s-%-2s-%-2s)' % (type1, type2, type3, type4), end='')
        print('%7.4f' % (oop.e), end='')
        print(' (%i-%i-%i-%i)' % (oop.at1+1, oop.at2+1, oop.at3+1, oop.at4+1))

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes_file(outfile, mol):
    outfile.write('# %i Outofplanes (At1, At2, At3, At4,' % (mol.n_outofplanes))
    outfile.write(' V_n, Gamma, N_f)\n')
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        outfile.write('OUTOFPLANE %4i %4i %4i %4i' % (oop.at1+1, oop.at2+1, oop.at3+1, oop.at4+1))
        outfile.write(' %6.2f %6.1f %1i\n' % (oop.v_n, oop.gam, oop.n_fold))

# print energy values to screen
def print_energy(mol):
    print('\n ------------ Energy Values -------------', end='')
    print('\n%10.4f kcal/mol Total Energy' % (mol.e_total), end='')
    print('\n%10.4f kcal/mol Kinetic Energy' % (mol.e_kinetic), end='')
    print('\n%10.4f kcal/mol Potential Energy' % (mol.e_potential), end='')
    print('\n%10.4f kcal/mol Non-bonded Energy' % (mol.e_nonbonded), end='')
    print('\n%10.4f kcal/mol Bonded Energy' % (mol.e_bonded), end='')
    print('\n%10.4f kcal/mol van der Waals Energy' % (mol.e_vdw), end='')
    print('\n%10.4f kcal/mol Electrostatic Energy' % (mol.e_elst), end='')
    print('\n%10.4f kcal/mol Bond Energy' % (mol.e_bonds), end='')
    print('\n%10.4f kcal/mol Angle Energy' % (mol.e_angles), end='')
    print('\n%10.4f kcal/mol Torsion Energy' % (mol.e_torsions), end='')
    print('\n%10.4f kcal/mol Out-of-plane Energy\n' % (mol.e_outofplanes), end='')
        
# check for proper input arguments and return result
def get_input():
    if (not len(sys.argv) > 1):
        print('\nUsage: mm.py INPUT_FILE\n')
        print('  INPUT_FILE: xyzq file for molecular mechanics input\n')
        sys.exit()
    else:
        infile_name = sys.argv[1]
    if ('debug' in sys.argv):
        for module in modules:
            os.system('rm %s.pyc' % (module))
    return infile_name
