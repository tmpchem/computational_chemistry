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

def GetBondGraph(atoms):
  """Build graph of which atoms are covalently bonded and bond lengths.
  
  Search all atom pairs to find those within a threshold of the sum of the
  interatomic covalent radii. Append those found to the bond graph dictionary of
  each Atom object within the Molecule object.
  
  First atomic index is the array index. Second atomic index is a dictionary
  key. The value is the bond length [Angstrom] of the bond.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of Atom objects with coordinate and MM
        parameter data.

  Returns:
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  bond_graph = {i:{} for i in range(len(atoms))}
  for i, j in itertools.combinations(range(len(atoms)), 2):
    at1, at2 = atoms[i], atoms[j]
    threshold = const.BONDTHRESHOLD * (at1.covrad + at2.covrad)
    r2_12 = geomcalc.GetR2ij(at1.coords, at2.coords)
    if r2_12 < threshold**2:
      r_12 = math.sqrt(r2_12)
      bond_graph[i][j] = bond_graph[j][i] = r_12
  return bond_graph


def GetBondGraphFromBonds(bonds, n_atoms):
  """Builds graph of which atoms are covalently bonded from list of bonds.

  Iterate over all bonds and add each atom pair to the dictionary of bond
  connectivity in the bond graph.

  Args:
    bonds (mmlib.molecule.Bond*): Array of molecule's Bond objects.
    n_atoms (int): Numer of atoms in molecule.

  Returns:
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  bond_graph = {i:{} for i in range(n_atoms)}
  for bond in bonds:
    bond_graph[bond.at1][bond.at2] = bond.r_ij
    bond_graph[bond.at2][bond.at1] = bond.r_ij
  return bond_graph


def GetBonds(atoms, bond_graph):
  """Determine covalently bonded atoms from bond graph and get parameters.
  
  Search bond graph for bonded atom pairs, and parameter tables for mm
  parameters. Use to create a new Bond object and append to Molecule.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int: float)): Dictionary of bond connectivity.

  Returns:
    bonds (mmlib.molecule.Bond*): Array of molecule's Bond objects.
  """
  bonds = []
  for i, at1 in enumerate(atoms):
    for j in bond_graph[i]:
      if i > j:
        continue
      at2 = atoms[j]
      r_ij = bond_graph[i][j]
      k_b, r_eq = param.GetBondParam(at1.type_, at2.type_)
      if k_b:
        bonds.append(molecule.Bond(i, j, r_ij, r_eq, k_b))
  bonds.sort(key=lambda b:(b.at1, b.at2))
  return bonds


def GetAngles(atoms, bond_graph):
  """Determine bond angle atom triplets from bond graph and get parameters.
  
  Search bond graph for bond angle triplets, and parameter tables for mm
  parameters. Use to create a new Angle object and append to Molecule.
  Count as angle if ij and jk are bonded in triplet ijk.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.

  Returns:
    angles (mmlib.molecule.Angle*): Array of molecule's Angle objects.
  """
  angles = []
  for j, at2 in enumerate(atoms):
    for i, k in itertools.combinations(bond_graph[j], 2):
      if i > k:
        continue
      at1, at3 = atoms[i], atoms[k]
      a_ijk = geomcalc.GetAijk(at1.coords, at2.coords, at3.coords)
      k_a, a_eq = param.GetAngleParam(at1.type_, at2.type_, at3.type_)
      if k_a:
        angles.append(molecule.Angle(i, j, k, a_ijk, a_eq, k_a))
  angles.sort(key=lambda a:(a.at1, a.at2, a.at3))
  return angles

    
def GetTorsions(atoms, bond_graph):
  """Determine torsion angle atom quartets and parameters from bond graph.
  
  Search bond graph for torsion angle quartets, and parameter tables for mm
  parameters. Use to create a new Torsion object and append to Molecule. Count
  as torsion if ij and jk and kl are bonded in quartet ijkl.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.

  Returns:
    torsions (mmlib.molecule.Torsion*): Array of molecule's Torsion objects.
  """
  torsions = []
  for j, at2 in enumerate(atoms):
    for i, k in itertools.permutations(bond_graph[j], 2):
      if j > k:
        continue
      at1, at3 = atoms[i], atoms[k]
      for l in bond_graph[k]:
        if l == i or l == j:
          continue
        at4 = atoms[l]
        t_ijkl = geomcalc.GetTijkl(
            at1.coords, at2.coords, at3.coords, at4.coords)
        params = param.GetTorsionParam(
            at1.type_, at2.type_, at3.type_, at4.type_)
        for v_n, gamma, nfold, paths in params:
          if v_n:
            torsions.append(
                molecule.Torsion(i, j, k, l, t_ijkl, v_n, gamma, nfold, paths))
  torsions.sort(key=lambda t:(t.at1, t.at2, t.at3, t.at4))
  return torsions

def GetOutofplanes(atoms, bond_graph):
  """Determine outofplane atom quartets and parameters from bond graph.
  
  Search bond graph for outofplane angle quartets, and parameter tables for mm
  parameters. Use to create a new Outofplane object and append to Molecule.
  Count as outofplane if il and jl and kl are bonded in quartet ijkl.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:foat)): Dictionary of bond connectivity.

  Returns:
    outofplanes (mmlib.molecule.Outofplane*): Array of molecule's Outofplane
        objects.
  """
  outofplanes = []
  for k, at3 in enumerate(atoms):
    for i, j, l in itertools.combinations(bond_graph[k], 3):
      at1, at2, at4 = atoms[i], atoms[j], atoms[l]
      # Three unique oop angles with this central atom (k) and atom quartet
      combos = [
        (min(i, j), max(i, j), k, l),
        (min(j, l), max(j, l), k, i),
        (min(l, i), max(l, i), k, j)]
      for combo in combos:
        o_ijkl = geomcalc.GetOijkl(*[atoms[x].coords for x in combo])
        v_n = param.GetOutofplaneParam(*[atoms[x].type_ for x in combo])
        if v_n:
          outofplanes.append(molecule.Outofplane(*combo, o_ijkl, v_n))
  outofplanes.sort(key=lambda o:(o.at1, o.at2, o.at3, o.at4))
  return outofplanes


def GetNonints(bonds, angles, torsions):
  """Determine which atomic pairs have bonded interactions.
  
  If two atoms belong to a mutual bond, angle, and/or torsion, then they are
  'bonded' and should not interact non-covalently. These pairs are appended to a
  set and ignored during computation of non-bonded interaction energies. It is
  not necessary to include outofplanes, since these atoms are already covalently
  linked through an Angle object.
  
  Args:
    bonds (mmlib.molecule.Bond*): Array of Bond objects.
    angles (mmlib.molecule.Angle*): Array of Angle objects.
    torsions (mmlib.molecule.Torsion*): Array of Torsion objects.

  Returns:
    nonints (set(int, int)): Set of atomic index pairs of non-interacting
        nonbonded atom pairs.
  """
  nonints = set()
  for bond in bonds:
    nonints.add((bond.at1, bond.at2))
    nonints.add((bond.at2, bond.at1))
  for angle in angles:
    nonints.add((angle.at1, angle.at3))
    nonints.add((angle.at3, angle.at1))
  for torsion in torsions:
    nonints.add((torsion.at1, torsion.at4))
    nonints.add((torsion.at4, torsion.at1))
  return nonints


def UpdateBonds(bonds, atoms, bond_graph):
  """Update all bond lengths [Angstrom] within a molecule object.
  
  Args:
    bonds (mmlib.molecule.Bond*): Array of molecule's Bond objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  for bond in bonds:
    c1 = atoms[bond.at1].coords
    c2 = atoms[bond.at2].coords
    bond.r_ij = geomcalc.GetRij(c1, c2)
    bond_graph[bond.at1][bond.at2] = bond.r_ij
    bond_graph[bond.at2][bond.at1] = bond.r_ij


def UpdateAngles(angles, atoms, bond_graph):
  """Update all bond angles [degrees] within a molecule object.
  
  Args:
    angles (mmlib.molecule.Angle*): Array of molecule's Angle objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  for angle in angles:
    r12 = bond_graph[angle.at1][angle.at2]
    r23 = bond_graph[angle.at2][angle.at3]
    c1 = atoms[angle.at1].coords
    c2 = atoms[angle.at2].coords
    c3 = atoms[angle.at3].coords
    angle.a_ijk = geomcalc.GetAijk(c1, c2, c3, r12, r23)


def UpdateTorsions(torsions, atoms, bond_graph):
  """Update all torsion angles [degrees] within a molecule object.
  
  Args:
    torsions (mmlib.molecule.Torsion*): Array of molecule's Torsion objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  for torsion in torsions:
    r12 = bond_graph[torsion.at1][torsion.at2]
    r23 = bond_graph[torsion.at2][torsion.at3]
    r34 = bond_graph[torsion.at3][torsion.at4]
    c1 = atoms[torsion.at1].coords
    c2 = atoms[torsion.at2].coords
    c3 = atoms[torsion.at3].coords
    c4 = atoms[torsion.at4].coords
    torsion.t_ijkl = geomcalc.GetTijkl(c1, c2, c3, c4, r12, r23, r34)


def UpdateOutofplanes(outofplanes, atoms, bond_graph):
  """Update all outofplane angles [degrees] within a molecule object.
  
  Args:
    outofplanes (mmlib.molecule.Outofplane*): Array of molecule's Outofplane
        objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  for outofplane in outofplanes:
    r31 = bond_graph[outofplane.at3][outofplane.at1]
    r32 = bond_graph[outofplane.at3][outofplane.at2]
    r34 = bond_graph[outofplane.at3][outofplane.at4]
    c1 = atoms[outofplane.at1].coords
    c2 = atoms[outofplane.at2].coords
    c3 = atoms[outofplane.at3].coords
    c4 = atoms[outofplane.at4].coords
    outofplane.o_ijkl = geomcalc.GetOijkl(c1, c2, c3, c4, r31, r32, r34)
