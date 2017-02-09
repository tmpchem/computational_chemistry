
"""Functions for determining bonded structure of molecular systems."""

from mmlib import geomcalc, param, molecule

def bond_thresh():
    """Threshold beyond covalent radii sum to determine bond cutoff"""
    return 1.2

def get_bond_graph(mol):
    """Build graph of which atoms are covalently bonded
    
    Search all atom pairs to find those within a threshold of the sum of
    the interatomic covalent radii. Append those found to the bond graph
    of each Atom object within the Molecule object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Atom objects
            with geometry and mm parameter data.
    """
    mol.bond_graph = [[] for i in range(mol.n_atoms)]
    bond_thresh = bond_thresh()
    for i in range(mol.n_atoms):
        for j in range(i+1, mol.n_atoms):
            thresh = bond_thresh * (mol.atoms[i].covrad + mol.atoms[j].covrad)
            r2_12 = geomcalc.get_r2_ij(mol.atoms[i].coords,
                mol.atoms[j].coords)
            if (r2_12 < thresh**2):
                mol.bond_graph[i].append(j)
                mol.bond_graph[j].append(i)

def get_bonds(mol):
    """Determine covalently bonded atoms from bond graph and get parameters.
    
    Search bond graph for bonded atom pairs, and parameter tables for mm
        parameters. Use to create a new Bond object and append to Molecule.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond graph of
            atom pairs within covalent radius cutoff threshold.
    """
    for i in range(mol.n_atoms):
        for a in range(len(mol.bond_graph[i])):
            j = mol.bond_graph[i][a]
            if (i < j):
                r_ij = geomcalc.get_r_ij(mol.atoms[i].coords,
                    mol.atoms[j].coords)
                k_b, r_eq = param.get_bond_param(mol.atoms[i].attype,
                    mol.atoms[j].attype)
                if (k_b > 0.0):
                    mol.bonds.append(molecule.Bond(i, j, r_ij, r_eq, k_b))
    mol.bonds = sorted(mol.bonds, key=lambda bond:bond.at2)
    mol.bonds = sorted(mol.bonds, key=lambda bond:bond.at1)
    mol.n_bonds = len(mol.bonds)
    
def get_angles(mol):
    """Determine bond angle atom triplets from bond graph and get parameters.
    
    Search bond graph for bond angle triplets, and parameter tables for mm
    parameters. Use to create a new Angle object and append to Molecule.
    Count as angle if ij and jk are bonded in triplet ijk.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond graph of
            atom pairs within covalent radius cutoff threshold.
    """
    for j in range(mol.n_atoms):
        n_jbonds = len(mol.bond_graph[j])
        for a in range(n_jbonds):
            i = mol.bond_graph[j][a]
            for b in range(a+1, n_jbonds):
                k = mol.bond_graph[j][b]
                a_ijk = geomcalc.get_a_ijk(mol.atoms[i].coords,
                    mol.atoms[j].coords, mol.atoms[k].coords)
                k_a, a_eq = param.get_angle_param(mol.atoms[i].attype,
                    mol.atoms[j].attype, mol.atoms[k].attype)
                if (k_a > 0.0):
                    mol.angles.append(molecule.Angle(i, j, k, a_ijk, a_eq,
                        k_a))
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at3)
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at2)
    mol.angles = sorted(mol.angles, key=lambda angle:angle.at1)
    mol.n_angles = len(mol.angles)
    
def get_torsions(mol):
    """Determine torsion angle atom quartets and parameters from bond graph.
    
    Search bond graph for torsion angle quartets, and parameter tables for
    mm parameters. Use to create a new Torsion object and append to Molecule.
    Count as torsion if ij and jk and kl are bonded in quartet ijkl.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond graph of
            atom pairs within covalent radius cutoff threshold.
    """
    for j in range(mol.n_atoms):
        n_jbonds = len(mol.bond_graph[j])
        for a in range(n_jbonds):
            k = mol.bond_graph[j][a]
            if (k < j): continue
            n_kbonds = len(mol.bond_graph[k])
            for b in range(n_jbonds):
                i = mol.bond_graph[j][b]
                if (i == k): continue
                for c in range(n_kbonds):
                    l = mol.bond_graph[k][c]
                    if (l == j or l == i): continue
                    t_ijkl = geomcalc.get_t_ijkl(mol.atoms[i].coords,
                        mol.atoms[j].coords, mol.atoms[k].coords,
                        mol.atoms[l].coords)
                    v_n, gamma, nfold, paths = param.get_torsion_param(
                        mol.atoms[i].attype, mol.atoms[j].attype,
                        mol.atoms[k].attype, mol.atoms[l].attype)
                    if (v_n > 0.0):
                        mol.torsions.append(molecule.Torsion(i, j, k, l,
                            t_ijkl, v_n, gamma, nfold, paths))
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at4)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at3)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at2)
    mol.torsions = sorted(mol.torsions, key=lambda torsion:torsion.at1)
    mol.n_torsions = len(mol.torsions)

def get_outofplanes(mol):
    """Determine outofplane atom quartets and parameters from bond graph.
    
    Search bond graph for outofplane angle quartets, and parameter tables for
    mm parameters. Use to create a new Outofplane object and append to
    Molecule. Count as outofplane if il and jl and kl are bonded in quartet
    ijkl.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond graph of
            atom pairs within covalent radius cutoff threshold.
    """
    for l in range(mol.n_atoms):
        n_lbonds = len(mol.bond_graph[l])
        for a in range(n_lbonds):
            i = mol.bond_graph[l][a]
            for b in range(n_lbonds):
                j = mol.bond_graph[l][b]
                if (i == j): continue
                for c in range(b+1, n_lbonds):
                    k = mol.bond_graph[l][c]
                    if (i == k): continue
                    o_ijkl = geomcalc.get_o_ijkl(mol.atoms[i].coords,
                        mol.atoms[j].coords, mol.atoms[k].coords,
                        mol.atoms[l].coords)
                    v_n = param.get_outofplane_param(mol.atoms[i].coords,
                        mol.atoms[j].coords, mol.atoms[k].coords,
                        mol.atoms[l].coords)
                    if (v_n > 0.0):
                        mol.outofplanes.append(molecule.Outofplane(i, j, k, l,
                            o_ijkl, v_n))
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at4)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at3)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at2)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda outofplane:outofplane.at1)
    mol.n_outofplanes = len(mol.outofplanes)

def get_nonints(mol):
    """Determine which atomic pairs have bonded interactions.
    
    If two atoms belong to a mutual bond, angle, and/or torsion, then
    they are 'bonded' and should not interact non-covalently. These
    pairs are appended to a list and ignored during computation of
    non-bonded interaction energies.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with arrays of
            Bond, Angle, and Torsion objects.
    """
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

# end of module

