import math, geomcalc

# energy.py: functions for calculating molecular mechanics energy of molecules 

# conversion of elst energy into kcal/mol
ceu2kcal = 332.06375

# calculate bond length energy between bonded atoms
def e_bond(r_ij, r_eq, k_b):
    e_bond = k_b * (r_ij - r_eq)**2
    return e_bond

# calculate bond angle energy between bonded atoms
def e_angle(a_ijk, a_eq, k_a):
    e_angle = k_a * (geomcalc.deg2rad() * (a_ijk - a_eq) )**2
    return e_angle

# calculate torsion angle energy between bonded atoms
def e_torsion(t_ijkl, v_n, gamma, n_fold, paths):
    e_torsion = v_n * (1 + math.cos(geomcalc.deg2rad() * (n_fold * t_ijkl - gamma))) / paths
    return e_torsion

# calculate out-of-plane angle (improper torsion) bewteen bonded atoms
def e_outofplane(o_ijkl, v_n, gamma, n_fold):
    e_outofplane = v_n * (1 + math.cos(geomcalc.deg2rad() * (n_fold * o_ijkl - gamma)))
    return e_outofplane

# calculate van der waals interaction between atom pair
def get_vdw_ij(r_ij, eps_ij, ro_ij):
    r6_ij = (ro_ij / r_ij) ** 6
    vdw_ij = eps_ij * ( r6_ij**2 - 2 * r6_ij )
    return vdw_ij

# calculate electrostatics interaction between atom pair
def get_elst_ij(r_ij, q_i, q_j, epsilon):
    elst_ij = ceu2kcal * ( q_i * q_j ) / ( epsilon * r_ij )
    return elst_ij

# calculate non-bonded interactions between all atoms
def get_e_nonbonded(mol):
    mol.e_nonbonded, mol.e_vdw, mol.e_elst = 0.0, 0.0, 0.0
    # van der waals and electrostatic energy
    for i in range(mol.n_atoms):
        for j in range(i+1, mol.n_atoms):
            if (j in mol.nonints[i]): continue
            r_ij = geomcalc.get_r_ij(mol.atoms[i].coords, mol.atoms[j].coords)
            eps_ij = mol.atoms[i].sreps * mol.atoms[j].sreps
            ro_ij = mol.atoms[i].ro + mol.atoms[j].ro
            mol.e_elst += get_elst_ij(r_ij, mol.atoms[i].charge, mol.atoms[j].charge, mol.dielectric)
            mol.e_vdw += get_vdw_ij(r_ij, eps_ij, ro_ij)

# update bond values and calculate bond energies
def get_e_bonds(mol):
    mol.e_bonds = 0.0
    for p in range(mol.n_bonds):
        at1, at2 = mol.bonds[p].at1, mol.bonds[p].at2
        mol.bonds[p].r_ij = geomcalc.get_r_ij(mol.atoms[at1].coords, mol.atoms[at2].coords)
        mol.bonds[p].e_bond = e_bond(mol.bonds[p].r_ij, mol.bonds[p].r_eq, mol.bonds[p].k_b)
        mol.e_bonds += mol.bonds[p].e_bond

# update angle values and calculate angle energies
def get_e_angles(mol):
    mol.e_angles = 0.0
    for p in range(mol.n_angles):
        at1, at2, at3 = mol.angles[p].at1, mol.angles[p].at2, mol.angles[p].at3
        mol.angles[p].a_ijk = geomcalc.get_a_ijk(mol.atoms[at1].coords, mol.atoms[at2].coords, mol.atoms[at3].coords)
        mol.angles[p].e_angle = e_angle(mol.angles[p].a_ijk, mol.angles[p].a_eq, mol.angles[p].k_a)
        mol.e_angles += mol.angles[p].e_angle

# update torsion values and calculate torsion energies
def get_e_torsions(mol):
    mol.e_torsions = 0.0
    for p in range(mol.n_torsions):
        at1, at2, at3, at4 = mol.torsions[p].at1, mol.torsions[p].at2, mol.torsions[p].at3, mol.torsions[p].at4
        mol.torsions[p].t_ijkl = geomcalc.get_t_ijkl(mol.atoms[at1].coords, mol.atoms[at2].coords, mol.atoms[at3].coords, mol.atoms[at4].coords)
        mol.torsions[p].e_torsion = e_torsion(mol.torsions[p].t_ijkl, mol.torsions[p].v_n, mol.torsions[p].gamma, mol.torsions[p].nfold, mol.torsions[p].paths)
        mol.e_torsions += mol.torsions[p].e_torsion

# update outofplane values and calculate outofplane energies
def get_e_outofplanes(mol):
    mol.e_outofplanes = 0.0
    for p in range(mol.n_outofplanes):
        at1, at2, at3, at4 = mol.outofplanes[p].at1, mol.outofplanes[p].at2, mol.outofplanes[p].at3, mol.outofplanes[p].at4
        mol.outofplanes[p].o_ijkl = geomcalc.get_o_ijkl(mol.atoms[at1].coords, mol.atoms[at2].coords, mol.atoms[at3].coords, mol.atoms[at4].coords)
        mol.outofplanes[p].e_outofplanes = e_outofplane(mol.outofplanes[p].o_ijkl, mol.outofplanes[p].v_n, mol.outofplanes[p].gamma, mol.outofplanes[p].nfold)
        mol.e_outofplanes += mol.outofplanes[p].e_outofplanes
        
# update aggregate energy values
def get_e_totals(mol):
    mol.e_bonded = mol.e_bonds + mol.e_angles + mol.e_torsions + mol.e_outofplanes
    mol.e_nonbonded = mol.e_vdw + mol.e_elst
    mol.e_total = mol.e_bonded + mol.e_nonbonded
