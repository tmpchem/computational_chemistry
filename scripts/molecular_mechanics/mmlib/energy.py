import math
from mmlib import geomcalc

# energy.py: functions for calculating molecular mechanics energy of molecules 

# conversion of elst energy from ceu to kcal/mol
def ceu2kcal(): return 332.06375

# conversion of kinetic energy from amu*A^2/ps^2 to kcal/mol
def kin2kcal(): return 2.39005736 * 10**-3

# calculate bond length energy between bonded atoms
def get_e_bond(r_ij, r_eq, k_b):
    e_bond = k_b * (r_ij - r_eq)**2
    return e_bond

# calculate bond angle energy between bonded atoms
def get_e_angle(a_ijk, a_eq, k_a):
    e_angle = k_a * (geomcalc.deg2rad() * (a_ijk - a_eq) )**2
    return e_angle

# calculate torsion angle energy between bonded atoms
def get_e_torsion(t_ijkl, v_n, gamma, n_fold, paths):
    e_torsion = v_n * (1.0 + math.cos(geomcalc.deg2rad() * (n_fold * t_ijkl - gamma))) / paths
    return e_torsion

# calculate out-of-plane angle (improper torsion) energy bewteen bonded atoms
def get_e_outofplane(o_ijkl, v_n, gamma, n_fold):
    e_outofplane = v_n * (1.0 + math.cos(geomcalc.deg2rad() * (n_fold * o_ijkl - gamma)))
    return e_outofplane

# calculate van der waals interaction between atom pair
def get_e_vdw_ij(r_ij, eps_ij, ro_ij):
    r6_ij = (ro_ij / r_ij) ** 6
    e_vdw_ij = eps_ij * ( r6_ij**2 - 2.0 * r6_ij )
    return e_vdw_ij

# calculate electrostatics interaction between atom pair
def get_e_elst_ij(r_ij, q_i, q_j, epsilon):
    e_elst_ij = ceu2kcal() * ( q_i * q_j ) / ( epsilon * r_ij )
    return e_elst_ij

# calculate boundary energy of an atom
def get_e_bound_i(k_box, bounds, coords):
    e_bound_i = 0.0
    for j in range(3):
        scale = 1.0 if (abs(coords[j]) >= bounds[j]) else 0.0
        e_bound_i += scale * k_box * (abs(coords[j]) - bounds[j])**2
    return e_bound_i

# calculate kinetic energy of an atom
def get_e_kinetic_i(mass, vels):
    e_kin_i = 0.0
    for i in range(3):
        e_kin_i += kin2kcal() * 0.5 * mass * vels[i]**2
    return e_kin_i

# calculate non-bonded interactions between all atoms
def get_e_nonbonded(mol):
    mol.e_nonbonded, mol.e_vdw, mol.e_elst = 0.0, 0.0, 0.0
    # van der waals and electrostatic energy
    for i in range(mol.n_atoms):
        atom1 = mol.atoms[i]
        for j in range(i+1, mol.n_atoms):
            atom2 = mol.atoms[j]
            if (j in mol.nonints[i]): continue
            r_ij = geomcalc.get_r_ij(atom1.coords, atom2.coords)
            eps_ij = atom1.sreps * atom2.sreps
            ro_ij = atom1.ro + atom2.ro
            mol.e_elst += get_e_elst_ij(r_ij, atom1.charge, atom2.charge, mol.dielectric)
            mol.e_vdw += get_e_vdw_ij(r_ij, eps_ij, ro_ij)

# update bond values and calculate bond energies
def get_e_bonds(mol):
    mol.e_bonds = 0.0
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        c1 = mol.atoms[bond.at1].coords
        c2 = mol.atoms[bond.at2].coords
        bond.r_ij = geomcalc.get_r_ij(c1, c2)
        bond.e = get_e_bond(bond.r_ij, bond.r_eq, bond.k_b)
        mol.e_bonds += bond.e

# update angle values and calculate angle energies
def get_e_angles(mol):
    mol.e_angles = 0.0
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        c1 = mol.atoms[ang.at1].coords
        c2 = mol.atoms[ang.at2].coords
        c3 = mol.atoms[ang.at3].coords
        ang.a_ijk = geomcalc.get_a_ijk(c1, c2, c3)
        ang.e = get_e_angle(ang.a_ijk, ang.a_eq, ang.k_a)
        mol.e_angles += ang.e

# update torsion values and calculate torsion energies
def get_e_torsions(mol):
    mol.e_torsions = 0.0
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        c1 = mol.atoms[tor.at1].coords
        c2 = mol.atoms[tor.at2].coords
        c3 = mol.atoms[tor.at3].coords
        c4 = mol.atoms[tor.at4].coords
        tor.t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4)
        tor.e = get_e_torsion(tor.t_ijkl, tor.v_n, tor.gam, tor.n, tor.paths)
        mol.e_torsions += tor.e

# update outofplane values and calculate outofplane energies
def get_e_outofplanes(mol):
    mol.e_outofplanes = 0.0
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        c1 = mol.atoms[oop.at1].coords
        c2 = mol.atoms[oop.at2].coords
        c3 = mol.atoms[oop.at3].coords
        c4 = mol.atoms[oop.at4].coords
        oop.o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4)
        oop.e = get_e_outofplane(oop.o_ijkl, oop.v_n, oop.gam, oop.n)
        mol.e_outofplanes += oop.e

# update box boundary energies
def get_e_bound(mol):
    mol.e_bound = 0.0
    k_box = mol.k_box
    bound = mol.bound
    for i in range(mol.n_atoms):
        atom = mol.atoms[i]
        k_box = mol.k_box
        atom.e_bound = get_e_bound_i(k_box, bound, atom.coords)
        mol.e_bound += atom.e_bound

# calculate kinetic energy
def get_e_kinetic(mol, kintype):
    mol.e_kinetic = 0.0
    if (kintype == 'nokinetic'):
        pass
    elif (kintype == 'leapfrog'):
        for p in range(mol.n_atoms):
            mass = mol.atoms[p].mass
            vels = 0.5*(mol.atoms[p].vels + mol.atoms[p].pvels)
            e_kin = get_e_kinetic_i(mass, vels)
            mol.e_kinetic += e_kin
    else:
        for p in range(mol.n_atoms):
            mass = mol.atoms[p].mass
            vels = mol.atoms[p].vels
            e_kin = get_e_kinetic_i(mass, vels)
            mol.e_kinetic += e_kin

# update total system energy values
def get_e_totals(mol):
    mol.e_bonded  = mol.e_bonds + mol.e_angles
    mol.e_bonded += mol.e_torsions + mol.e_outofplanes
    mol.e_nonbonded = mol.e_vdw + mol.e_elst
    mol.e_potential = mol.e_bonded + mol.e_nonbonded + mol.e_bound
    mol.e_total = mol.e_kinetic + mol.e_potential

