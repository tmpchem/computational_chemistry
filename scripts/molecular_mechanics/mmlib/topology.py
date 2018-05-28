"""Functions for determining bonded structure of molecular systems.

Includes functions for determining bonds, angles, torsions, and outofplanes for
an mmlib.molecule.Molecule object through atoms which are located within a
threshold of the sum of interatomic covalent radii.
"""

import itertools
import math

from mmlib import constants as const
from mmlib import geomcalc
from mmlib import molecule
from mmlib import param

def GetBondGraph(mol):
  """Build graph of which atoms are covalently bonded and bond lengths.
  
  Search all atom pairs to find those within a threshold of the sum of the
  interatomic covalent radii. Append those found to the bond graph dictionary of
  each Atom object within the Molecule object.
  
  First atomic index is the array index. Second atomic index is a dictionary
  key. The value is the bond length [Angstrom] of the bond.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with Atom objects with
        geometry and mm parameter data.
  """
  mol.bond_graph = {i:{} for i in range(mol.n_atoms)}
  bond_thresh = const.BONDTHRESHOLD
  for i, j in itertools.combinations(range(mol.n_atoms), 2):
    at1, at2 = mol.atoms[i], mol.atoms[j]
    thresh = bond_thresh * (at1.covrad + at2.covrad)
    r2_12 = geomcalc.GetR2ij(at1.coords, at2.coords)
    if r2_12 < thresh**2:
      r_12 = math.sqrt(r2_12)
      mol.bond_graph[i][j] = r_12
      mol.bond_graph[j][i] = r_12


def GetBonds(mol):
  """Determine covalently bonded atoms from bond graph and get parameters.
  
  Search bond graph for bonded atom pairs, and parameter tables for mm
  parameters. Use to create a new Bond object and append to Molecule.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with bond graph of atom pairs
        within covalent radius cutoff threshold.
  """
  for i, at1 in enumerate(mol.atoms):
    for j in mol.bond_graph[i]:
      if i < j:
        continue
      at2 = mol.atoms[j]
      r_ij = mol.bond_graph[i][j]
      k_b, r_eq = param.GetBondParam(at1.type_, at2.type_)
      if k_b:
        mol.bonds.append(molecule.Bond(i, j, r_ij, r_eq, k_b))
  mol.bonds = sorted(mol.bonds, key=lambda b:(b.at1, b.at2))
  mol.n_bonds = len(mol.bonds)


def GetAngles(mol):
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
    for i in mol.bond_graph[j]:
      at1 = mol.atoms[i]
      r_ij = mol.bond_graph[i][j]
      for k in mol.bond_graph[j]:
        if i >= k:
          continue
        at3 = mol.atoms[k]
        r_jk = mol.bond_graph[j][k]
        a_ijk = geomcalc.GetAijk(at1.coords, at2.coords, at3.coords, r_ij, r_jk)
        k_a, a_eq = param.GetAngleParam(at1.type_, at2.type_, at3.type_)
        if k_a:
            mol.angles.append(molecule.Angle(i, j, k, a_ijk, a_eq, k_a))
  mol.angles = sorted(mol.angles, key=lambda a:(a.at1, a.at2, a.at3))
  mol.n_angles = len(mol.angles)

    
def GetTorsions(mol):
  """Determine torsion angle atom quartets and parameters from bond graph.
  
  Search bond graph for torsion angle quartets, and parameter tables for mm
  parameters. Use to create a new Torsion object and append to Molecule. Count
  as torsion if ij and jk and kl are bonded in quartet ijkl.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with bond graph of atom pairs
        within covalent radius cutoff threshold.
  """
  for j in range(mol.n_atoms):
    at2 = mol.atoms[j]
    for k in mol.bond_graph[j]:
      if j >= k:
        continue
      at3 = mol.atoms[k]
      r_jk = mol.bond_graph[j][k]
      for i in mol.bond_graph[j]:
        if i == j or i == k:
          continue
        at1 = mol.atoms[i]
        r_ij = mol.bond_graph[i][j]
        for l in mol.bond_graph[k]:
          if i == l or j == l or k == l:
            continue
          at4 = mol.atoms[l]
          r_kl = mol.bond_graph[k][l]
          t_ijkl = geomcalc.GetTijkl(at1.coords, at2.coords, at3.coords,
                                     at4.coords, r_ij, r_jk, r_kl)
          params = param.GetTorsionParam(at1.type_, at2.type_, at3.type_,
                                         at4.type_)
          for v_n, gamma, nfold, paths in params:
            if v_n:
              mol.torsions.append(molecule.Torsion(i, j, k, l, t_ijkl, v_n,
                                                   gamma, nfold, paths))
  mol.torsions = sorted(mol.torsions, key=lambda t:(t.at1, t.at2, t.at3, t.at4))
  mol.n_torsions = len(mol.torsions)


def GetOutofplanes(mol):
  """Determine outofplane atom quartets and parameters from bond graph.
  
  Search bond graph for outofplane angle quartets, and parameter tables for mm
  parameters. Use to create a new Outofplane object and append to Molecule.
  Count as outofplane if il and jl and kl are bonded in quartet ijkl.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with bond graph of atom pairs
        within covalent radius cutoff threshold.
  """
  for k in range(mol.n_atoms):
    at3 = mol.atoms[k]
    for l in mol.bond_graph[k]:
      if k == l:
        continue
      at4 = mol.atoms[l]
      r34 = mol.bond_graph[k][l]
      for i in mol.bond_graph[k]:
        if i == k or i == l:
          continue
        at1 = mol.atoms[i]
        r31 = mol.bond_graph[k][i]
        for j in mol.bond_graph[k]:
          if j >= i or j == k or j == l:
            continue
          at2 = mol.atoms[j]
          r32 = mol.bond_graph[k][j]
          o_ijkl = geomcalc.GetOijkl(at1.coords, at2.coords, at3.coords,
                                     at4.coords, r31, r32, r34)
          v_n = param.GetOutofplaneParam(at1.type_, at2.type_, at3.type_,
                                         at4.type_)
          if v_n:
            mol.outofplanes.append(molecule.Outofplane(i, j, k, l, o_ijkl, v_n))
  mol.outofplanes = sorted(mol.outofplanes,
                           key=lambda o:(o.at1, o.at2, o.at3, o.at4))
  mol.n_outofplanes = len(mol.outofplanes)


def GetNonints(mol):
  """Determine which atomic pairs have bonded interactions.
  
  If two atoms belong to a mutual bond, angle, and/or torsion, then they are
  'bonded' and should not interact non-covalently. These pairs are appended to a
  set and ignored during computation of non-bonded interaction energies. It is
  not necessary to include outofplanes, since these atoms are already covalently
  linked through an Angle object.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with arrays of Bond, Angle,
        and Torsion objects.
  """
  mol.nonints = set()
  for p in range(mol.n_bonds):
    b = mol.bonds[p]
    mol.nonints.add((b.at1, b.at2))
    mol.nonints.add((b.at2, b.at1))
  for p in range(mol.n_angles):
    a = mol.angles[p]
    mol.nonints.add((a.at1, a.at3))
    mol.nonints.add((a.at3, a.at1))
  for p in range(mol.n_torsions):
    t = mol.torsions[p]
    mol.nonints.add((t.at1, t.at4))
    mol.nonints.add((t.at4, t.at1))


def UpdateBonds(mol):
  """Update all bond lengths [Angstrom] within a molecule object.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with bond data.
  """
  for bond in mol.bonds:
    c1 = mol.atoms[bond.at1].coords
    c2 = mol.atoms[bond.at2].coords
    bond.r_ij = geomcalc.GetRij(c1, c2)
    mol.bond_graph[bond.at1][bond.at2] = bond.r_ij
    mol.bond_graph[bond.at2][bond.at1] = bond.r_ij


def UpdateAngles(mol):
  """Update all bond angles [degrees] within a molecule object.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with angle data.
  """
  for angle in mol.angles:
    r12 = mol.bond_graph[angle.at1][angle.at2]
    r23 = mol.bond_graph[angle.at2][angle.at3]
    c1 = mol.atoms[angle.at1].coords
    c2 = mol.atoms[angle.at2].coords
    c3 = mol.atoms[angle.at3].coords
    angle.a_ijk = geomcalc.GetAijk(c1, c2, c3, r12, r23)


def UpdateTorsions(mol):
  """Update all torsion angles [degrees] within a molecule object.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with torsion data.
  """
  for tors in mol.torsions:
    r12 = mol.bond_graph[torsion.at1][torsion.at2]
    r23 = mol.bond_graph[torsion.at2][torsion.at3]
    r34 = mol.bond_graph[torsion.at3][torsion.at4]
    c1 = mol.atoms[torsion.at1].coords
    c2 = mol.atoms[torsion.at2].coords
    c3 = mol.atoms[torsion.at3].coords
    c4 = mol.atoms[torsion.at4].coords
    torsion.t_ijkl = geomcalc.GetTijkl(c1, c2, c3, c4, r12, r23, r34)


def UpdateOutofplanes(mol):
  """Update all outofplane angles [degrees] within a molecule object.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with outofplane data.
  """
  for outofplane in mol.outofplanes:
    r31 = mol.bond_graph[outofplane.at3][outofplane.at1]
    r32 = mol.bond_graph[outofplane.at3][outofplane.at2]
    r34 = mol.bond_graph[outofplane.at3][outofplane.at4]
    c1 = mol.atoms[outofplane.at1].coords
    c2 = mol.atoms[outofplane.at2].coords
    c3 = mol.atoms[outofplane.at3].coords
    c4 = mol.atoms[outofplane.at4].coords
    outofplane.o_ijkl = geomcalc.GetOijkl(c1, c2, c3, c4, r31, r32, r34)
