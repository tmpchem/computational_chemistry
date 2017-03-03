
"""Functions for determining bonded structure of molecular systems.

Includes functions for determining bonds, angles, torsions, and
outofplanes for an mmlib.molecule.Molecule object through atoms which
are located within a threshold of the sum of interatomic covalent
radii.
"""

import math
from mmlib import geomcalc, param, molecule

def bond_threshold():
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
    mol.bond_graph = [{} for i in range(mol.n_atoms)]
    bond_thresh = bond_threshold()
    for i in range(mol.n_atoms):
        at1 = mol.atoms[i]
        for j in range(i+1, mol.n_atoms):
            at2 = mol.atoms[j]
            thresh = bond_thresh * (at1.covrad + at2.covrad)
            r2_12 = geomcalc.get_r2_ij(at1.coords, at2.coords)
            if (r2_12 < thresh**2):
                r_12 = math.sqrt(r2_12)
                mol.bond_graph[i][j] = r_12
                mol.bond_graph[j][i] = r_12

def get_bonds(mol):
    """Determine covalently bonded atoms from bond graph and get parameters.
    
    Search bond graph for bonded atom pairs, and parameter tables for mm
        parameters. Use to create a new Bond object and append to Molecule.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond graph of
            atom pairs within covalent radius cutoff threshold.
    """
    for i in range(mol.n_atoms):
        at1 = mol.atoms[i]
        for j in mol.bond_graph[i].keys():
            at2 = mol.atoms[j]
            if (i < j):
                r_ij = mol.bond_graph[i][j]
                k_b, r_eq = param.get_bond_param(at1.type, at2.type)
                if (k_b > 0.0):
                    mol.bonds.append(molecule.Bond(i, j, r_ij, r_eq, k_b))
    mol.bonds = sorted(mol.bonds, key=lambda b:b.at2)
    mol.bonds = sorted(mol.bonds, key=lambda b:b.at1)
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
        at2 = mol.atoms[j]
        for i in mol.bond_graph[j].keys():
            at1 = mol.atoms[i]
            r_ij = mol.bond_graph[i][j]
            for k in mol.bond_graph[j].keys():
                if (i >= k):
                    continue
                at3 = mol.atoms[k]
                r_jk = mol.bond_graph[j][k]
                a_ijk = geomcalc.get_a_ijk(at1.coords, at2.coords, at3.coords,
                    r_ij, r_jk)
                k_a, a_eq = param.get_angle_param(at1.type, at2.type,
                    at3.type)
                if (k_a > 0.0):
                    mol.angles.append(molecule.Angle(i, j, k, a_ijk, a_eq,
                        k_a))
    mol.angles = sorted(mol.angles, key=lambda a:a.at3)
    mol.angles = sorted(mol.angles, key=lambda a:a.at2)
    mol.angles = sorted(mol.angles, key=lambda a:a.at1)
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
        at2 = mol.atoms[j]
        for k in mol.bond_graph[j].keys():
            if (j >= k):
                continue
            at3 = mol.atoms[k]
            r_jk = mol.bond_graph[j][k]
            for i in mol.bond_graph[j].keys():
                if (i == j or i == k):
                    continue
                at1 = mol.atoms[i]
                r_ij = mol.bond_graph[i][j]
                for l in mol.bond_graph[k].keys():
                    if (i == l or j == l or k == l):
                        continue
                    at4 = mol.atoms[l]
                    r_kl = mol.bond_graph[k][l]
                    t_ijkl = geomcalc.get_t_ijkl(at1.coords, at2.coords,
                        at3.coords, at4.coords, r_ij, r_jk, r_kl)
                    params = param.get_torsion_param(at1.type, at2.type,
                        at3.type, at4.type)
                    for p in range(len(params)):
                        v_n, gamma, nfold, paths = params[p]
                        if (v_n > 0.0):
                            mol.torsions.append(molecule.Torsion(i, j, k, l,
                                t_ijkl, v_n, gamma, nfold, paths))
    mol.torsions = sorted(mol.torsions, key=lambda t:t.at4)
    mol.torsions = sorted(mol.torsions, key=lambda t:t.at3)
    mol.torsions = sorted(mol.torsions, key=lambda t:t.at2)
    mol.torsions = sorted(mol.torsions, key=lambda t:t.at1)
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
    for k in range(mol.n_atoms):
        at3 = mol.atoms[k]
        for l in mol.bond_graph[k].keys():
            if (k == l):
                continue
            at4 = mol.atoms[l]
            r34 = mol.bond_graph[k][l]
            for i in mol.bond_graph[k].keys():
                if (i == k or i == l):
                    continue
                at1 = mol.atoms[i]
                r31 = mol.bond_graph[k][i]
                for j in mol.bond_graph[k].keys():
                    if (j >= i or j == k or j == l):
                        continue
                    at2 = mol.atoms[j]
                    r32 = mol.bond_graph[k][j]
                    o_ijkl = geomcalc.get_o_ijkl(at1.coords, at2.coords,
                        at3.coords, at4.coords, r31, r32, r34)
                    v_n = param.get_outofplane_param(at1.type, at2.type,
                        at3.type, at4.type)
                    if (v_n > 0.0):
                        mol.outofplanes.append(molecule.Outofplane(i, j, k, l,
                            o_ijkl, v_n))
    mol.outofplanes = sorted(mol.outofplanes, key=lambda o:o.at4)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda o:o.at3)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda o:o.at2)
    mol.outofplanes = sorted(mol.outofplanes, key=lambda o:o.at1)
    mol.n_outofplanes = len(mol.outofplanes)

def get_nonints(mol):
    """Determine which atomic pairs have bonded interactions.
    
    If two atoms belong to a mutual bond, angle, and/or torsion, then
    they are 'bonded' and should not interact non-covalently. These
    pairs are appended to a list and ignored during computation of
    non-bonded interaction energies.
    
    It is not necessary to include outofplanes, since these atoms are
    already covalently linked through an Angle object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with arrays of
            Bond, Angle, and Torsion objects.
    """
    mol.nonints = [[] for i in range(mol.n_atoms)]
    for p in range(mol.n_bonds):
        at1 = mol.bonds[p].at1
        at2 = mol.bonds[p].at2
        mol.nonints[at1].append(at2)
        mol.nonints[at2].append(at1)
    for p in range(mol.n_angles):
        at1 = mol.angles[p].at1
        at3 = mol.angles[p].at3
        mol.nonints[at1].append(at3)
        mol.nonints[at3].append(at1)
    for p in range(mol.n_torsions):
        at1 = mol.torsions[p].at1
        at4 = mol.torsions[p].at4
        mol.nonints[at1].append(at4)
        mol.nonints[at4].append(at1)
    for i in range(mol.n_atoms):
        mol.nonints[i] = sorted(list(set(mol.nonints[i])))

def update_bonds(mol):
    """Update all bond lengths [Angstrom] within a molecule object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with bond data.
    """
    for p in range(mol.n_bonds):
        b = mol.bonds[p]
        c1 = mol.atoms[b.at1].coords
        c2 = mol.atoms[b.at2].coords
        b.r_ij = geomcalc.get_r_ij(c1, c2)
        mol.bond_graph[b.at1][b.at2] = b.r_ij
        mol.bond_graph[b.at2][b.at1] = b.r_ij

def update_angles(mol):
    """Update all bond angles [degrees] within a molecule object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with angle data.
    """
    for p in range(mol.n_angles):
        a = mol.angles[p]
        r12 = mol.bond_graph[a.at1][a.at2]
        r23 = mol.bond_graph[a.at2][a.at3]
        c1 = mol.atoms[a.at1].coords
        c2 = mol.atoms[a.at2].coords
        c3 = mol.atoms[a.at3].coords
        a.a_ijk = geomcalc.get_a_ijk(c1, c2, c3, r12, r23)

def update_torsions(mol):
    """Update all torsion angles [degrees] within a molecule object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with torsion data.
    """
    for p in range(mol.n_torsions):
        t = mol.torsions[p]
        r12 = mol.bond_graph[t.at1][t.at2]
        r23 = mol.bond_graph[t.at2][t.at3]
        r34 = mol.bond_graph[t.at3][t.at4]
        c1 = mol.atoms[t.at1].coords
        c2 = mol.atoms[t.at2].coords
        c3 = mol.atoms[t.at3].coords
        c4 = mol.atoms[t.at4].coords
        t.t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4, r12, r23, r34)

def update_outofplanes(mol):
    """Update all outofplane angles [degrees] within a molecule object.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with outofplane data.
    """
    for p in range(mol.n_outofplanes):
        o = mol.outofplanes[p]
        r31 = mol.bond_graph[o.at3][o.at1]
        r32 = mol.bond_graph[o.at3][o.at2]
        r34 = mol.bond_graph[o.at3][o.at4]
        c1 = mol.atoms[o.at1].coords
        c2 = mol.atoms[o.at2].coords
        c3 = mol.atoms[o.at3].coords
        c4 = mol.atoms[o.at4].coords
        o.o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4, r31, r32, r34)

# end of module

