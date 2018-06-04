"""Functions for computing molecular mechanics energy gradient components.

Includes conversion factors, energy gradient magnitudes for individual
structural objects (atoms, bonds, angles, torsions, outofplanes), system
gradient components, and total energy gradient member data for
mmlib.molecule.Molecule objects.
"""

import itertools
import math
import numpy

from mmlib import constants as const
from mmlib import geomcalc
from mmlib import molecule

def GetGMagBond(r_ij, r_eq, k_b):
  """Calculate energy gradient magnitude of bond stretch.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    r_eq (float): Equilibrium bond length [Angstrom] of bond ij.
    k_b (float): Spring constant [kcal/(mol*A^2)] of bond ij.
  
  Returns:
      g_bond (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  return 2.0 * k_b * (r_ij - r_eq)


def GetGMagAngle(a_ijk, a_eq, k_a):
  """Calculate energy gradient magnitude of angle bend.
  
  Args:
    a_ijk (float): Angle [degrees] between atoms i, j, and k.
    a_eq (float): Equilibrium bond angle [degrees] of angle ijk.
    k_a (float): Spring constant [kcal/(mol*rad^2)] of angle ijk.
  
  Returns:
    g_angle (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  return 2.0 * k_a * (const.DEG2RAD * (a_ijk - a_eq) )


def GetGMagTorsion(t_ijkl, v_n, gamma, n_fold, paths):
  """Calculate energy gradient magnitude of torsion strain.
  
  Args:
    t_ijkl (float): Torsion [degrees] between atoms i, j, k, and l.
    v_n (float): Half-barrier height [kcal/mol] of torsion ijkl.
    gamma (float): Barrier offset [degrees] of torsion ijkl.
    n_fold (int): Barrier frequency of torsion ijkl.
    paths (int): Number of distinct paths in torsion ijkl.

  Returns:
    g_torsion (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  return -v_n * n_fold * math.sin(
      const.DEG2RAD * (n_fold * t_ijkl - gamma)) / paths


def GetGMagOutofplane(o_ijkl, v_n):
  """Calculate energy gradient magnitude of outofplane bend.
  
  Args:
    o_ijkl (float): Outofplane angle [degrees] between atoms i, j, k, and l.
    v_n (float): Half-barrier height [kcal/mol] of outofplane ijkl.
  
  Returns:
    g_outofplane (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  return -2.0 * v_n * math.sin(const.DEG2RAD * (2.0 * o_ijkl - 180.0))


def GetGMagVdwIJ(r_ij, eps_ij, ro_ij):
  """Calculate energy gradient magnitude of van der waals pair energy.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    eps_ij (float): Van der Waals epsilon [kcal/mol] between pair ij.
    ro_ij (float): Van der Waals radius [Angstrom] between pair ij.
  
  Returns:
    g_vdw_ij (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  rrel_ij = ro_ij / r_ij
  return -12.0 * (eps_ij / ro_ij) * (rrel_ij**13 - rrel_ij**7)


def GetGMagElstIJ(r_ij, q_i, q_j, epsilon):
  """Calculate energy gradient magnitude of electrostatic pair energy.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    q_i (float): Partial charge [e] of atom i.
    q_j (float): Partial charge [e] of atom j.
    epsilon (float): Dielectric constant of space (>= 1.0).
  
  Returns:
    e_elst_ij (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  return -const.CEU2KCAL * q_i * q_j / (epsilon * r_ij**2)


def GetGMagBoundI(k_box, bound, coord, origin, boundtype):
  """Calculate energy gradient magnitude of boundary energy.
  
  Args:
    k_box (float): Spring constant [kcal/(mol*A^2)] of boundary.
    bound (float): Distance from origin [Angstrom] of boundary.
    coords (float*): Array of cartesian coordinates [Angstrom] of atom.
    origin (float*): Array of cartesian coordiantes [Angstrom] of origin of
        simulation.
    boundtype (str): 'cube' or 'sphere', type of boundary condition.
  
  Returns:
    g_bound_i (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  g_bound_i = numpy.zeros(const.NUMDIM)
  if boundtype == 'cube':
    for j in range(const.NUMDIM):
      sign = 1.0 if (coord[j] - origin[j]) <= 0.0 else -1.0
      scale = float(abs(coord[j] - origin[j]))
      g_bound_i[j] = (-2.0 * sign * scale * k_box * (abs(coord[j]) - bound))
  elif boundtype == 'sphere':
    r_io = geomcalc.GetRij(origin, coord)
    u_io = geomcalc.GetUij(origin, coord, u_io)
    scale = float(r_io >= bound)
    g_bound_i = 2.0 * scale * k_box * (r_io - bound) * u_io
  return g_bound_i


def GetGDirInter(coords1, coords2, r_12=None):
  """Calculate direction of energy gradient between atom pair.
  
  Args:
    coords1 (float*): 3 Cartesian coordinates [Angstrom] of atom1.
    coords2 (float*): 3 Cartesian coordiantes [Angstrom] of atom2.
    r_ij (float): Distance between atom1 and atom2 (default None).
  
  Returns:
    gdir1 (float*), gdir2 (float*): unit vectors in the direction of max
        increasing inter-atomic distance.
  """
  gdir1 = geomcalc.GetUij(coords2, coords1, r_12)
  gdir2 = -1.0 * gdir1
  return gdir1, gdir2


def GetGDirAngle(coords1, coords2, coords3, r_21=None, r_23=None):
  """Calculate direction of energy gradients between bond angle atoms.
  
  Args:
    coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
    coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
    coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
    r_21 (float): Distance between atom2 and atom1 (default None).
    r_23 (float): Distance between atom2 and atom3 (default None).

  Returns:
    gdir1 (float*), gdir2 (float*), gdir3 (float*): vectors in the direction of
        max increasing bond angle.
  """
  u_21 = geomcalc.GetUij(coords2, coords1, r_21)
  u_23 = geomcalc.GetUij(coords2, coords3, r_23)
  cp = geomcalc.GetUcp(u_21, u_23)
  gdir1 = geomcalc.GetUcp(u_21, cp) / r_21
  gdir3 = geomcalc.GetUcp(cp, u_23) / r_23
  gdir2 = -1.0 * (gdir1 + gdir3)
  return gdir1, gdir2, gdir3


def GetGDirTorsion(coords1, coords2, coords3, coords4, r_12=None, r_23=None,
                   r_34=None):
  """Calculate direction of energy gradients between torsion atoms.
  
  Args:
    coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
    coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
    coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
    coords4 (float*): 3 cartesian coordinates [Angstrom] of atom4.
    r_12 (float): Distance between atom1 and atom2 (default None).
    r_23 (float): Distance between atom2 and atom3 (default None).
    r_34 (float): Distance between atom3 and atom4 (default None).
  
  Returns:
    gdir1 (float*), gdir2 (float*), gdir3 (float*), gdir4 (float*): Vectors in
        the direction of max increasing torsion angle.
  """
  u_21 = geomcalc.GetUij(coords2, coords1, r_12)
  u_34 = geomcalc.GetUij(coords3, coords4, r_34)
  u_23 = geomcalc.GetUij(coords2, coords3, r_23)
  u_32 = -1.0 * u_23
  a_123 = geomcalc.GetAijk(coords1, coords2, coords3, r_12, r_23)
  a_432 = geomcalc.GetAijk(coords4, coords3, coords2, r_34, r_23)
  s_123 = math.sin(const.DEG2RAD * a_123)
  s_432 = math.sin(const.DEG2RAD * a_432)
  c_123 = math.cos(const.DEG2RAD * a_123)
  c_432 = math.cos(const.DEG2RAD * a_432)
  gdir1 = geomcalc.GetUcp(u_21, u_23) / (r_12*s_123)
  gdir4 = geomcalc.GetUcp(u_34, u_32) / (r_34*s_432)
  gdir2 = (r_12/r_23*c_123 - 1.0)*gdir1 - (r_34/r_23*c_432)*gdir4
  gdir3 = (r_34/r_23*c_432 - 1.0)*gdir4 - (r_12/r_23*c_123)*gdir1
  return gdir1, gdir2, gdir3, gdir4


def GetGDirOutofplane(coords1, coords2, coords3, coords4, oop, r_31=None,
                      r_32=None, r_34=None):
  """Calculate direction of energy gradients between outofplane atoms.
  
  Args:
    coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
    coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
    coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
    coords4 (float*): 3 cartesian coordinates [Angstrom] of atom4.
    oop (float): Out-of-plane angles bewteen atoms 1, 2, 3, and 4.
    r_31 (float): Distance between atom3 and atom1 (default None).
    r_32 (float): Distance between atom3 and atom2 (default None).
    r_34 (float): Distance between atom3 and atom4 (default None).
  
  Returns:
    gdir1 (float*), gdir2 (float*), gdir3 (float*), gdir4 (float*): Vectors in
        the direction of max increasing outofplane angle.
  """
  u_31 = geomcalc.GetUij(coords3, coords1, r_31)
  u_32 = geomcalc.GetUij(coords3, coords2, r_32)
  u_34 = geomcalc.GetUij(coords3, coords4, r_34)
  cp_3234 = geomcalc.GetCp(u_32, u_34)
  cp_3431 = geomcalc.GetCp(u_34, u_31)
  cp_3132 = geomcalc.GetCp(u_31, u_32)
  a_132 = geomcalc.GetAijk(coords1, coords3, coords2)
  s_132 = math.sin(const.DEG2RAD * a_132)
  c_132 = math.cos(const.DEG2RAD * a_132)
  c_oop = math.cos(const.DEG2RAD * oop)
  t_oop = math.tan(const.DEG2RAD * oop)
  gdir1 = ((1.0 / r_31) * (cp_3234 / (c_oop*s_132)
      - (t_oop/s_132**2) * (u_31 - c_132*u_32)))
  gdir2 = ((1.0 / r_32) * (cp_3431 / (c_oop*s_132)
      - (t_oop/s_132**2) * (u_32 - c_132*u_31)))
  gdir4 = ((1.0 / r_34) * (cp_3132 / (c_oop*s_132) - (t_oop*u_34)))
  gdir3 = -1.0 * (gdir1 + gdir2 + gdir4)
  return gdir1, gdir2, gdir3, gdir4


def GetGBonds(g_bonds, bonds, atoms):
  """Calculate bond length energy gradients for all bonds.
  
  Args:
    g_bonds (float**): Nx3 array of molecule's atomic bond gradients
        [kcal/(mol*A)].
    bonds (mmlib.molecule.Bond*): Array of molecule's Bond objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
  """
  g_bonds.fill(0.0)
  for bond in bonds:
    bond.GetGradientMagnitude()
    c1 = atoms[bond.at1].coords
    c2 = atoms[bond.at2].coords
    dir1, dir2 = GetGDirInter(c1, c2, bond.r_ij)
    g_bonds[bond.at1] += bond.grad_mag * dir1
    g_bonds[bond.at2] += bond.grad_mag * dir2


def GetGAngles(g_angles, angles, atoms, bond_graph):
  """Calculate angle bend energy gradients for all angles.
  
  Args:
    g_angles (float**): Nx3 array of molecule's atomic angle gradients
        [kcal/(mol*A)].
    angles (mmlib.molecule.Angle*): Array of molecule's Angle objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  g_angles.fill(0.0)
  for angle in angles:
    angle.GetGradientMagnitude()
    c1 = atoms[angle.at1].coords
    c2 = atoms[angle.at2].coords
    c3 = atoms[angle.at3].coords
    r12 = bond_graph[angle.at1][angle.at2]
    r23 = bond_graph[angle.at2][angle.at3]
    dir1, dir2, dir3 = GetGDirAngle(c1, c2, c3, r12, r23)
    g_angles[angle.at1] += angle.grad_mag * dir1
    g_angles[angle.at2] += angle.grad_mag * dir2
    g_angles[angle.at3] += angle.grad_mag * dir3


def GetGTorsions(g_torsions, torsions, atoms, bond_graph):
  """Calculate torsion strain energy gradients for all torsions.
  
  Args:
    g_torsions (float**): Nx3 array of molecule's atomic torsion gradients
        [kcal/(mol*A)].
    torsions (mmlib.molecule.Angle*): Array of molecule's Torsion objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  g_torsions.fill(0.0)
  for torsion in torsions:
    torsion.GetGradientMagnitude()
    c1 = atoms[torsion.at1].coords
    c2 = atoms[torsion.at2].coords
    c3 = atoms[torsion.at3].coords
    c4 = atoms[torsion.at4].coords
    r12 = bond_graph[torsion.at1][torsion.at2]
    r23 = bond_graph[torsion.at2][torsion.at3]
    r34 = bond_graph[torsion.at3][torsion.at4]
    dir1, dir2, dir3, dir4 = GetGdirTorsion(c1, c2, c3, c4, r12, r23, r34)
    g_torsions[torsion.at1] += torsion.grad_mag * dir1
    g_torsions[torsion.at2] += torsion.grad_mag * dir2
    g_torsions[torsion.at3] += torsion.grad_mag * dir3
    g_torsions[torsion.at4] += torsion.grad_mag * dir4


def GetGOutofplanes(g_outofplanes, outofplanes, atoms, bond_graph):
  """Calculate outofplane bend energy gradients for all outofplanes.
  
  Args:
    g_outofplanes (float**): Nx3 array of molecule's atomic outofplane gradients
        [kcal/(mol*A)].
    outofplanes (mmlib.molecule.Angle*): Array of molecule's Outofplane objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    bond_graph (int:(int:float)): Dictionary of bond connectivity.
  """
  g_outofplanes.fill(0.0)
  for outofplane in outofplanes:
    outofplane.GetGradientMagnitude()
    c1 = atoms[outofplane.at1].coords
    c2 = atoms[outofplane.at2].coords
    c3 = atoms[outofplane.at3].coords
    c4 = atoms[outofplane.at4].coords
    r31 = bond_graph[outofplane.at3][outofplane.at1]
    r32 = bond_graph[outofplane.at3][outofplane.at2]
    r34 = bond_graph[outofplane.at3][outofplane.at4]
    dir1, dir2, dir3, dir4 = GetGdirOutofplane(
        c1, c2, c3, c4, outofplane.o_ijkl, r31, r32, r34)
    g_outofplanes[outofplane.at1] += outofplane.grad_mag * dir1
    g_outofplanes[outofplane.at2] += outofplane.grad_mag * dir2
    g_outofplanes[outofplane.at3] += outofplane.grad_mag * dir3
    g_outofplanes[outofplane.at4] += outofplane.grad_mag * dir4


def GetGNonbonded(g_vdw, g_elst, atoms, nonints, dielectric):
    """Calculate non-bonded energy gradients between all nonbonded atom pairs.
    
    Computes van der waals and electrostatic energy gradient [kcal/(mol*A)]
    components between all pairs of non-bonded atoms in a system.

    Args:
      g_vdw (float**): Nx3 array of molecule's van der waals gradients.
      g_elst (float**): Nx3 array of molecule's electrostatic gradients.
      atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
      nonints (set(int, int)): Set of atomic index pairs of atoms without
          nonbonded interactions due to covalent (near-)ajecency.
      dielectric (float): Dielectric constant of molecule.
    """
    g_vdw.fill(0.0)
    g_elst.fill(0.0)
    for i, j in itertools.combinations(range(len(atoms)), 2):
      if (i, j) in nonints:
        continue
      atom1, atom2 = atoms[i], atoms[j]
      r_ij = geomcalc.GetRij(atom1.coords, atom2.coords)
      dir1, dir2 = GetGDirInter(atom1.coords, atom2.coords, r_ij)
      eps_ij = atom1.sreps * atom2.sreps
      ro_ij = atom1.ro + atom2.ro
      g_elst_mag = GetGMagElstIJ(r_ij, atom1.charge, atom2.charge, dielectric)
      g_vdw_mag = GetGMagVdwIJ(r_ij, eps_ij, ro_ij)
      g_vdw[i] += g_vdw_mag * dir1
      g_vdw[j] += g_vdw_mag * dir2
      g_elst[i] += g_elst_mag * dir1
      g_elst[j] += g_elst_mag * dir2


def GetGBound(g_bound, atoms, k_box, boundary, origin, boundary_type):
  """Calculate boundary energy gradients for all atoms.
  
  Args:
    g_bound (float**): Nx3 array of molecule's boundary gradients.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    k_box (float): Spring constant [kcal/(mol*A^2)] of molecule boundary.
    boundary (float): Distance [Angstrom] from origin to molecule boundary.
    origin (float*): Cartesian coordinates of molecule origin.
    boundary_type (str): Molecular boundary type ('sphere' or 'cube').
  """
  g_bound.fill(0.0)
  for i, atom in enumerate(atoms):
    g_bound[i] += GetGBoundI(k_box, bound, atom.coords, origin, boundtype)


def GetVirial(g_total, atoms):
  """Clausius virial function for all atomic coordinates and forces.

  Args:
    g_total (float**): Nx3 Array of molecular energy gradient [kcal/(mol*A)].
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    virial (float): Clausius virial of molecule.
  """
  virial = 0.0
  for i, atom in enumerate(atoms):
    for j in range(const.NUMDIM):
      virial += atom.coords[j] * g_total[i][j]
  return virial


def GetPressure(n_atoms, temperature, virial, volume):
  """Update total pressure of a system of molecules with boundary.
  
  Args:
    n_atoms (int): Number of atoms in molecule.
    temperature (float): Temperature of molecule [Kelvin].
    virial (float): Clausius virial function [kcal/mol] of molecule.
    volume (float): Volume of system [Angstrom^3].

  Returns:
    pressure (float): Pressure of system [Pascals].
  """
  return const.KCALAMOL2PA * (
      n_atoms * const.KB * temperature + virial / const.NUMDIM) / volume


def GetGNumerical(mol):
  """Update total numerical energy gradient [kcal/(mol*A)] of all atoms.

  Args:
    mol (mmlib.molecule.Molecule): Molecule object with energy gradient
        component data [kcal/(mol*A)].
  """
  mol.g_bonds.fill(0.0)
  mol.g_angles.fill(0.0)
  mol.g_torsions.fill(0.0)
  mol.g_outofplanes.fill(0.0)
  mol.g_vdw.fill(0.0)
  mol.g_elst.fill(0.0)
  mol.g_bound.fill(0.0)

  for i in range(mol.n_atoms):
    for j in range(const.NUMDIM):
      q = mol.atoms[i].coords[j]

      # Displace in positive direction and compute energy.
      qp = q + 0.5 * const.NUMDISP
      mol.atoms[i].coords[j] = qp
      mol.UpdateInternals()
      mol.GetEnergy()
      ep_bond, ep_ang, ep_tor, ep_oop, ep_vdw, ep_elst, ep_bound = (
          mol.e_bonds, mol.e_angles, mol.e_torsions, mol.e_outofplanes,
          mol.e_vdw, mol.e_elst, mol.e_bound)

      # Displace in negative direction and compute energy.
      qm = q - 0.5 * const.NUMDISP
      mol.atoms[i].coords[j] = qm
      mol.UpdateInternals()
      mol.GetEnergy()
      em_bond, em_ang, em_tor, em_oop, em_vdw, e_elst, em_bound = (
          mol.e_bonds, mol.e_angles, mol.e_torsions, mol.e_outofplanes,
          mol.e_vdw, mol.e_elst, mol.e_bound)

      # Return to original coordinate and compute all gradient components.
      mol.atoms[i].coords[j] = q
      mol.g_bonds[i][j] = (ep_bond - em_bond) / disp
      mol.g_angles[i][j] = (ep_ang - em_ang) / disp
      mol.g_torsions[i][j] = (ep_tor - em_tor) / disp
      mol.g_outofplanes[i][j] = (ep_oop - em_oop) / disp
      mol.g_vdw[i][j] = (ep_vdw - em_vdw) / disp
      mol.g_elst[i][j] = (ep_elst - em_elst) / disp
      mol.g_bound[i][j] = (ep_bound - em_bound) / disp

  mol.UpdateInternals()
