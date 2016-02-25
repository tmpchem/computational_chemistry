import os, sys
import param, geomcalc, topology, energy, molecule

# fileio.py: functions for reading in and printing out molecular mechanics data

# list of modules imported from mmlib
modules = ['__init__', 'fileio', 'param', 'geomcalc', 'topology', 'energy', 'molecule']

# create 2d array of strings from input file name
def get_file_string_array(infile_name):
    try:
        infile = open(infile_name, 'r')
    except IOError:
        print 'Error: file (%s) does not exist!' % (infile_name)
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
        at_coords = []
        for j in range(3):
            at_coords.append(float(infile_array[i+2][j+1]))
        at_charge = float(infile_array[i+2][4])
        at_element = at_type[0].capitalize()
        at_mass = param.get_at_mass(at_element)
        at_ro, at_eps = param.get_vdw_param(at_type)
        new_atom = molecule.atom(at_type, at_coords, at_charge, at_ro, at_eps, at_mass)
        new_atom.set_covrad(param.get_cov_rad(at_element))
        mol.atoms.append(new_atom)

# print geometry and non-bonded parameters for a set of atoms
def print_geom(mol):
    print '\n ------------- Molecular Geoemetry and Non-bonded Parameters --------------'
    print '     | Atom      X            Y            Z',
    print '          Q       Ro/2    Eps'
    print     ' --------------------------------------------------------------------------'
    for i in range(mol.n_atoms):
        print '%4i | %-2s' % (i+1, mol.atoms[i].attype),
        for j in range(3):
            print '%12.6f' % (mol.atoms[i].coords[j]),
        print '  %8.5f %7.4f %7.4f' % (mol.atoms[i].charge, mol.atoms[i].ro, mol.atoms[i].eps)

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
        print '\n --------------------- Bond Length Data ----------------------'
        print '     |    k_b      r_eq      r_ij     types   energy   atoms'
        print   ' -------------------------------------------------------------'
    else:
        print '\n No Bonds Detected'
    for p in range(mol.n_bonds):
        at1, at2 = mol.bonds[p].at1, mol.bonds[p].at2
        print '%4i | %7.2f  %8.4f  %8.4f  ' % (p+1, mol.bonds[p].k_b, mol.bonds[p].r_eq, mol.bonds[p].r_ij),
        print '(%-2s-%-2s) %7.4f' % (mol.atoms[at1].attype, mol.atoms[at2].attype, mol.bonds[p].e_bond),
        print ' (%i-%i)' % (at1+1, at2+1)

# print bond topology and bond parameters to parameter file
def print_bonds_file(outfile, mol):
    outfile.write('# %i Bonds (At1, At2, K_b, R_eq)\n' % (mol.n_bonds))
    for p in range(mol.n_bonds):
        at1, at2 = mol.bonds[p].at1, mol.bonds[p].at2
        outfile.write('BOND %4i %4i %7.2f %7.4f\n' % (at1+1, at2+1, mol.bonds[p].k_b, mol.bonds[p].r_eq))

# print bond angle topology and angle parameters from an array
def print_angles(mol):
    if (mol.n_angles > 0):
        print '\n ------------------- Bond Angle Data ---------------------'
        print '     |   k_a    a_eq     a_ijk    types     energy   atoms'
        print   ' ---------------------------------------------------------'
    else:
        print '\n No Bond Angles Detected'
    for p in range(mol.n_angles):
        at1, at2, at3 = mol.angles[p].at1, mol.angles[p].at2, mol.angles[p].at3
        print '%4i | %6.2f  %7.3f  %7.3f' % (p+1, mol.angles[p].k_a, mol.angles[p].a_eq, mol.angles[p].a_ijk),
        print '(%-2s-%-2s-%-2s)' % (mol.atoms[at1].attype, mol.atoms[at2].attype, mol.atoms[at3].attype),
        print '%7.4f' % (mol.angles[p].e_angle),
        print ' (%i-%i-%i)' % (at1+1, at2+1, at3+1)

# print bond angle topology and angle parameters to parameter file
def print_angles_file(outfile, mol):
    outfile.write('# %i Angles (At1, At2, At3, K_a, A_eq)\n' % (mol.n_angles))
    for p in range(mol.n_angles):
        at1, at2, at3 = mol.angles[p].at1, mol.angles[p].at2, mol.angles[p].at3
        outfile.write('ANGLE %4i %4i %4i %7.4f %8.4f\n' % (at1+1, at2+1, at3+1, mol.angles[p].k_a, mol.angles[p].a_eq))

# print torsion topology and torsion parameters from an array
def print_torsions(mol):
    if (mol.n_torsions > 0):
        print '\n ------------------- Torsion Angle Data ----------------------'
        print '     |   Vn/2  gamma   t_ijkl n p     types      energy   atoms'
        print   ' -------------------------------------------------------------'
    else:
        print '\n No Torsion Angles Detected'
    for p in range(mol.n_torsions):
        at1, at2, at3, at4 = mol.torsions[p].at1, mol.torsions[p].at2, mol.torsions[p].at3, mol.torsions[p].at4
        print '%4i | %6.2f %6.1f %8.3f %1i %1i' % (p+1, mol.torsions[p].v_n, mol.torsions[p].gamma, mol.torsions[p].t_ijkl, mol.torsions[p].nfold, mol.torsions[p].paths),
        print '(%-2s-%-2s-%-2s-%-2s)' % (mol.atoms[at1].attype, mol.atoms[at2].attype, mol.atoms[at3].attype, mol.atoms[at4].attype),
        print '%7.4f' % (mol.torsions[p].e_torsion),
        print ' (%i-%i-%i-%i)' % (at1+1, at2+1, at3+1, at4+1)

# print torsion topology and torsion parameters from an array
def print_torsions_file(outfile, mol):
    outfile.write('# %i Torsions (At1, At2, At3, At4,' % (mol.n_torsions))
    outfile.write('V_n, Gamma, N_f, paths)\n')
    for p in range(mol.n_torsions):
        at1, at2, at3, at4 = mol.torsions[p].at1, mol.torsions[p].at2, mol.torsions[p].at3, mol.torsions[p].at4
        outfile.write('TORSION %4i %4i %4i %4i' % (at1+1, at2+1, at3+1, at4+1))
        outfile.write(' %6.2f %6.1f %1i %1i\n' % (mol.torsions[p].v_n, mol.torsions[p].gamma, mol.torsions[p].n_fold, mol.torsions[p].paths))

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes(mol):
    if (mol.n_outofplanes > 0):
        print '\n ---------------- Out-of-plane Angle Data ----------------'
        print '     |  Vn/2   o_ijkl      types      energy   atoms'
        print   ' ---------------------------------------------------------'
    else:
        print '\n No Out-of-plane Angles Detected'
    for p in range(mol.n_outofplanes):
        at1, at2, at3, at4 = mol.outofplanes[p].at1, mol.outofplanes[p].at2, mol.outofplanes[p].at3, mol.outofplanes[p].at4
        print '%4i | %6.2f %7.3f' % (p+1, mol.outofplanes[p].v_n, mol.outofplanes[p].o_ijkl),
        print '  (%-2s-%-2s-%-2s-%-2s)' % (mol.atoms[at1].attype, mol.atoms[at2].attype, mol.atoms[at3].attype, mol.atoms[at4].attype),
        print '%7.4f' % (mol.outofplanes[p].e_outofplane),
        print ' (%i-%i-%i-%i)' % (at1+1, at2+1, at3+1, at4+1)

# print out-of-plane angles and out-of-plane parameters from an array
def print_outofplanes_file(outfile, mol):
    outfile.write('# %i Outofplanes (At1, At2, At3, At4,' % (mol.n_outofplanes))
    outfile.write(' V_n, Gamma, N_f)\n')
    for p in range(mol.n_outofplanes):
        at1, at2, at3, at4 = mol.outofplanes[p].at1, mol.outofplanes[p].at2, mol.outofplanes[p].at3, mol.outofplanes[p].at4
        outfile.write('OUTOFPLANE %4i %4i %4i %4i' % (at1+1, at2+1, at3+1, at4+1))
        outfile.write(' %6.2f %6.1f %1i\n' % (mol.outofplanes[p].v_n, mol.outofplanes[p].gamma, mol.outofplanes[p].n_fold))

# print energy values to screen
def print_energy(mol):
    print '\n ------------ Energy Values -------------',
    print '\n%10.4f kcal/mol Total Energy' % (mol.e_total),
    print '\n%10.4f kcal/mol Non-bonded Energy' % (mol.e_nonbonded),
    print '\n%10.4f kcal/mol Bonded Energy' % (mol.e_bonded),
    print '\n%10.4f kcal/mol van der Waals Energy' % (mol.e_vdw),
    print '\n%10.4f kcal/mol Electrostatic Energy' % (mol.e_elst),
    print '\n%10.4f kcal/mol Bond Energy' % (mol.e_bonds),
    print '\n%10.4f kcal/mol Angle Energy' % (mol.e_angles),
    print '\n%10.4f kcal/mol Torsion Energy' % (mol.e_torsions),
    print '\n%10.4f kcal/mol Out-of-plane Energy\n' % (mol.e_outofplanes),
        
# check for proper input arguments and return result
def get_input():
    if (not len(sys.argv) > 1):
        print '\nUsage: mm.py INPUT_FILE'
        print '\nINPUT_FILE: xyzq file to be used for molecular mechanics input\n'
        sys.exit()
    else:
        infile_name = sys.argv[1]
    if ('debug' in sys.argv):
        for module in modules:
            os.system('rm %s.pyc' % (module))
    return infile_name
