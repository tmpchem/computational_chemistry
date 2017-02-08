from mmlib import geomcalc, param, molecule

# topology.py: functions for determining the bonded structure of molecules

# threshold beyond average of covalent radii to determine bond cutoff
bond_thresh = 1.2

# build tree of which atoms are covalently bonded
def get_bond_tree(mol):
    mol.bond_tree = [[] for i in range(mol.n_atoms)]
    for i in range(mol.n_atoms):
        for j in range(i+1, mol.n_atoms):
            thresh = bond_thresh * (mol.atoms[i].covrad + mol.atoms[j].covrad)
            r2_12 = geomcalc.get_r2_ij(mol.atoms[i].coords, mol.atoms[j].coords)
            if (r2_12 < thresh**2):
                mol.bond_tree[i].append(j)
                mol.bond_tree[j].append(i)

# determine atoms which are covalently bonded from bond tree and get parameters
def get_bonds(mol):
    for i in range(mol.n_atoms):
        for a in range(len(mol.bond_tree[i])):
            j = mol.bond_tree[i][a]
            if (i < j):
                r_ij = geomcalc.get_r_ij(mol.atoms[i].coords, mol.atoms[j].coords)
                k_b, r_eq = param.get_bond_param(mol.atoms[i].attype, mol.atoms[j].attype)
                if (k_b > 0.0):
                    mol.bonds.append(molecule.bond(i, j, r_ij, r_eq, k_b))
    mol.bonds = sorted(mol.bonds, key=lambda bond:bond.at2)
    mol.bonds = sorted(mol.bonds, key=lambda bond:bond.at1)
    mol.n_bonds = len(mol.bonds)
    
# determine atoms which form a bond angle from bond tree and get parameters
def get_angles(mol):
    for j in range(mol.n_atoms):
        n_jbonds = len(mol.bond_tree[j])
        for a in range(n_jbonds):
            i = mol.bond_tree[j][a]
            for b in range(a+1, n_jbonds):
                k = mol.bond_tree[j][b]
                a_ijk = geomcalc.get_a_ijk(mol.atoms[i].coords, mol.atoms[j].coords, mol.atoms[k].coords)
                k_a, a_eq = param.get_angle_param(mol.atoms[i].attype, mol.atoms[j].attype, mol.atoms[k].attype)
                if (k_a > 0.0):
                    mol.angles.append(molecule.angle(i, j, k, a_ijk, a_eq, k_a))
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at3)
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at2)
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at1)
    mol.n_angles = len(mol.angles)
    
# determine atoms which form torsion angles from bond tree and get parameters
def get_torsions(mol):
    for j in range(mol.n_atoms):
        n_jbonds = len(mol.bond_tree[j])
        for a in range(n_jbonds):
            k = mol.bond_tree[j][a]
            if (k < j): continue
            n_kbonds = len(mol.bond_tree[k])
            for b in range(n_jbonds):
                i = mol.bond_tree[j][b]
                if (i == k): continue
                for c in range(n_kbonds):
                    l = mol.bond_tree[k][c]
                    if (l == j or l == i): continue
                    t_ijkl = geomcalc.get_t_ijkl(mol.atoms[i].coords, mol.atoms[j].coords, mol.atoms[k].coords, mol.atoms[l].coords)
                    v_n, gamma, nfold, paths = param.get_torsion_param(mol.atoms[i].attype, mol.atoms[j].attype, mol.atoms[k].attype, mol.atoms[l].attype)
                    if (v_n > 0.0):
                        mol.torsions.append(molecule.torsion(i, j, k, l, t_ijkl, v_n, gamma, nfold, paths))
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at4)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at3)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at2)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at1)
    mol.n_torsions = len(mol.torsions)

# determine atoms which form out-of-plane angles from bond tree and get parameters
def get_outofplanes(mol):
    for l in range(mol.n_atoms):
        n_lbonds = len(mol.bond_tree[l])
        for a in range(n_lbonds):
            i = mol.bond_tree[l][a]
            for b in range(n_lbonds):
                j = mol.bond_tree[l][b]
                if (i == j): continue
                for c in range(b+1, n_lbonds):
                    k = mol.bond_tree[l][c]
                    if (i == k): continue
                    o_ijkl = geomcalc.get_o_ijkl(mol.atoms[i].coords, mol.atoms[j].coords, mol.atoms[k].coords, mol.atoms[l].coords)
                    v_n, gamma, nfold = param.get_outofplane_param(mol.atoms[i].coords, mol.atoms[j].coords, mol.atoms[k].coords, mol.atoms[l].coords)
                    if (v_n > 0.0):
                        mol.outofplanes.append(molecule.outofplane(i, j, k, l, o_ijkl, v_n, gamma, nfold))
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at4)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at3)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at2)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at1)
    mol.n_outofplanes = len(mol.outofplanes)

# determine atoms which have bonded interactions and no non-covalent interaction
def get_nonints(mol):
    mol.nonints = [[] for i in range(mol.n_atoms)]
    for p in range(mol.n_bonds):
        mol.nonints[mol.bonds[p].at1].append(mol.bonds[p].at2)
        mol.nonints[mol.bonds[p].at2].append(mol.bonds[p].at1)
    for p in range(mol.n_angles):
        mol.nonints[mol.angles[p].at1].append(mol.angles[p].at3)
        mol.nonints[mol.angles[p].at3].append(mol.angles[p].at1)
    for p in range(mol.n_torsions):
        mol.nonints[mol.torsions[p].at1].append(mol.torsions[p].at4)
        mol.nonints[mol.torsions[p].at4].append(mol.torsions[p].at1)
    for i in range(mol.n_atoms):
        mol.nonints[i] = sorted(list(set(mol.nonints[i])))

