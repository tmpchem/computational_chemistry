"""Functions for computing molecular mechanics energy components.

Includes conversion factors, energy functions for individual structural
objects (atoms, bonds, angles, torsions, outofplanes), system energy
components, and total energy member data for mmlib.molecule.Molecule
objects.
"""

import itertools
import math

from mmlib import constants as const
from mmlib import geomcalc

def GetEBond(r_ij, r_eq, k_b):
  """Calculate bond stretch energy between 2 bonded atoms.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    r_eq (float): Equilibrium bond length [Angstrom] of bond ij.
    k_b (float): Spring constant [kcal/(mol*A^2)] of bond ij.
      
  Returns:
    e_bond (float): Energy [kcal/mol] of bond ij.
  """
  return k_b * (r_ij - r_eq)**2


def GetEAngle(a_ijk, a_eq, k_a):
  """Calculate angle bend energy between 3 bonded atoms.
  
  Args:
    a_ijk (float): Angle [degrees] between atoms i, j, and k.
    a_eq (float): Equilibrium bond angle [degrees] of angle ijk.
    k_a (float): Spring constant [kcal/(mol*rad^2)] of angle ijk.
  
  Returns:
    e_angle (float): Energy [kcal/mol] of angle ijk.
  """
  return k_a * (const.DEG2RAD * (a_ijk - a_eq) )**2


def GetETorsion(t_ijkl, v_n, gamma, nfold, paths):
  """Calculate torsion strain energy between 4 bonded atoms.
  
  Args:
    t_ijkl (float): Torsion [degrees] between atoms i, j, k, and l.
    v_n (float): Half-barrier height [kcal/mol] of torsion ijkl.
    gamma (float): Barrier offset [degrees] of torsion ijkl.
    nfold (int): Barrier frequency of torsion ijkl.
    paths (int): Number of distinct paths in torsion ijkl.
  
  Returns:
    e_torsion (float): Energy [kcal/mol] of torsion ijkl.
  """
  return v_n * (1.0 + math.cos(const.DEG2RAD * (nfold*t_ijkl - gamma))) / paths


def GetEOutofplane(o_ijkl, v_n):
  """Calculate outofplane bend energy between 4 bonded atoms.
  
  Args:
    o_ijkl (float): Outofplane angle [degrees] between atoms i, j, k, and l.
    v_n (float): Half-barrier height [kcal/mol] of torsion ijkl.
  
  Returns:
    e_outofplane (float): Energy [kcal/mol] of outofplane ijkl.
  """
  return v_n * (1.0 + math.cos(const.DEG2RAD * (2.0 * o_ijkl - 180.0)))


def GetEVdwIJ(r_ij, eps_ij, ro_ij):
  """Calculate van der waals interaction energy between atom pair.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    eps_ij (float): Van der Waals epsilon [kcal/mol] between pair ij.
    ro_ij (float): Van der Waals radius [Angstrom] between pair ij.
  
  Returns:
    e_vdw_ij (float): Van der waals energy [kcal/mol] between pair ij.
  """
  r6_ij = (ro_ij / r_ij)**6
  return eps_ij * ( r6_ij**2 - 2.0 * r6_ij )


def GetEElstIJ(r_ij, q_i, q_j, epsilon):
  """Calculate electrostatic interaction energy between atom pair.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    q_i (float): Partial charge [e] of atom i.
    q_j (float): Partial charge [e] of atom j.
    epsilon (float): Dielectric constant of space (>= 1.0).
  
  Returns:
    e_elst_ij (float): Electrostatic energy [kcal/mol] between pair ij.
  """
  return const.CEU2KCAL * q_i * q_j / (epsilon * r_ij)


def GetEBoundI(k_box, bound, coords, origin, boundtype):
  """Calculate simulation boundary energy of an atom.
  
  Args:
    k_box (float): Spring constant [kcal/(mol*A^2)] of boundary.
    bound (float): Distance from origin [Angstrom] of boundary.
    coords (float*): Array of cartesian coordinates [Angstrom] of atom.
    origin (float*): Array of cartesian coordiantes [Angstrom] of origin of
        simulation.
    boundtype (str): 'cube' or 'sphere', type of boundary condition.
  
  Returns:
    e_bound_i (float): Boundary energy [kcal/mol] of atom.
  """
  e_bound_i = 0.0
  if (boundtype == 'cube'):
    for j in range(const.NUMDIM):
      scale = float(abs(coords[j] - origin[j]) >= bound)
      e_bound_i += scale * k_box * (abs(coords[j] - origin[j]) - bound)**2
  elif (boundtype == 'sphere'):
    r_io = geomcalc.GetRij(origin, coords)
    scale = float(r_io >= bound)
    e_bound_i += scale * k_box * (r_io - bound)**2
  return e_bound_i


def GetEKineticI(mass, vels):
  """Calculate kinetic energy of an atom
  
  Args:
    mass (float): Mass [g/mol] of atom.
    vels (float*): Array of velocities [Angstrom/ps] of atom.
  
  Returns:
    e_kin_i (float): Kinetic energy [kcal/mol] of atom.
  """
  e_kin_i = 0.0
  for i in range(const.NUMDIM):
    e_kin_i += mass * vels[i]**2
  return 0.5 * const.KIN2KCAL * e_kin_i


def GetEBonds(bonds):
  """Update bond length values and compute bond energy of system.
  
  For all bonds in a molecule, update the distance between atoms ij, and add the
      energy to the total molecular bond energy [kcal/mol].
  
  Args:
    bonds (mmlib.molecule.Bond*): Array of Bond objects containing atomic
        indices and molecular mechanics parameters.

  Returns:
    e_bonds (float): Bond energy [kcal/mol] of all bonds.
  """
  e_bonds = 0.0
  for bond in bonds:
    bond.GetEnergy()
    e_bonds += bond.energy
  return e_bonds


def GetEAngles(angles):
  """Update bond angle values and compute angle energy of system.
  
  For all angles in a molecule, update the angle between atoms ijk, and add the
  energy to the total molecular angle energy [kcal/mol].
  
  Args:
    angles (mmlib.molecule.Angle*): Array of Angle objects containing atomic
        indices and molecular mechanics parameters.

  Returns:
    e_angles (float): Angle energy [kcal/mol] of all bond angles.
  """
  e_angles = 0.0
  for angle in angles:
    angle.GetEnergy()
    e_angles += angle.energy
  return e_angles


def GetETorsions(torsions):
  """Update torsion angle values and compute torsion energy of system.
  
  For all torsions in a molecule, update the torsion between atoms ijkl, and add
  the energy to the total molecular torsion energy [kcal/mol].
  
  Args:
    torsions (mmlib.molecule.Torsion*): Array of Torsion objects containing
        atomic indices and molecular mechanics parameters.

  Returns:
    e_torsions (float): Torsion energy [kcal/mol] of all torsion angles.
  """
  e_torsions = 0.0
  for torsion in torsions:
    torsion.GetEnergy()
    e_torsions += torsion.energy
  return e_torsions


def GetEOutofplanes(outofplanes):
  """Update outofplane values and compute outofplane energy of system.
  
  For all outofplanes in molecule, update outofplane between atoms ijkl, and add
  the energy to the total molecular outofplane energy [kcal/mol].
  
  Args:
    outofplanes (mmlib.molecule.Outofplane*): Array of Outofplane objects
        containing atomic indices and molecular mechanics parameters.

  Return:
    e_outofplanes (float): Outofplane [kcal/mol] of all outofplane angles.
  """
  e_outofplanes = 0.0
  for outofplane in outofplanes:
    outofplane.GetEnergy()
    e_outofplanes += outofplane.energy
  return e_outofplanes


def GetENonbonded(atoms, nonints, dielectric):
  """Calculate non-bonded interaction energy between all atom pairs.
  
  Computes van der waals and electrostatic energy [kcal/mol] components
  between all pairs of non-bonded atoms in a system.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of Atom objects containing Cartesian
        coordinates and molecular mechanics parameters.
    nonints (set(int, int)): Set of atomic index pairs of atoms without
        nonbonded interactions due to covalent (near-)adjacency.
    dielectric (float): Dielectric constant of molecule.

  Returns:
    e_vdw (float): Van der waals energy [kcal/mol] of molecule.
    e_elst (float): Electrostatic energy [kcal/mol] of molecule.
  """
  e_vdw, e_elst = 0.0, 0.0
  for i, j in itertools.combinations(range(len(atoms)), 2):
    if (i, j) in nonints:
      continue
    atom1, atom2 = atoms[i], atoms[j]
    r_ij = geomcalc.GetRij(atom1.coords, atom2.coords)
    eps_ij = atom1.sreps * atom2.sreps
    ro_ij = atom1.ro + atom2.ro
    e_elst += GetEElstIJ(r_ij, atom1.charge, atom2.charge, dielectric)
    e_vdw += GetEVdwIJ(r_ij, eps_ij, ro_ij)
  return e_vdw, e_elst


def GetEBound(atoms, k_box, boundary, origin, boundary_type):
  """Compute total boundary energy of system.
  
  For all atoms in molecule, add the boundary energy to the total boundary
  energy [kcal/mol].
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of Atom objects containing Cartesian
        coordinates.
    k_box (float): Spring constant [kcal/(mol*A^2)] of molecule boundary.
    boundary (float): Distance [Angstrom] from origin to molecule boundary.
    origin (float*): Cartesian coordinates of molecule origin.
    boundary_type (str): Molecular boundary type ('sphere' or 'cube').

  Returns:
    e_bound (float): Boundary energy [kcal/mol] of molecule.
  """
  e_bound = 0.0
  for atom in atoms:
    e_bound += GetEBoundI(k_box, boundary, atom.coords, origin, boundary_type)
  return e_bound


def GetEKinetic(atoms, kintype=None):
  """Compute kinetic energy of all atoms in molecule.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of Atom objects containing (float) 3
        cartesian velocities [Angstrom/ps].
    kintype (str): Type of kinetic energy to be computed:
      'nokinetic': Do nothing (ke = 0.0).
      'leapfrog': Average of current and previous velocities.
      [default]: Use current velocities.

  Returns:
    e_kinetic (float): Kinetic energy [kcal/mol] of molecule.
  """
  if kintype == 'nokinetic':
    return 0.0
  
  e_kinetic = 0.0
  if kintype == 'leapfrog':
    for atom in atoms:
      vels = 0.5 * (atom.vels + atom.pvels)
      e_kinetic += GetEKineticI(atom.mass, vels)
  else:
    for atom in atoms:
      e_kinetic += GetEKineticI(atom.mass, atom.vels)
  return e_kinetic


def GetTemperature(e_kinetic, n_atoms):
  """Update kinetic temperature using current kinetic energy per atom.
  
  Args:
    e_kinetic (float): Kinetic energy [kcal/mol] of molecule.
    n_atoms (int): Number of atoms in molecule.

  Returns:
    temperature (float): Temperature [Kelvin] of molecule.
  """
  
  return (2.0 / const.NUMDIM) * e_kinetic / (const.KB * n_atoms)
