"""Functions for computing molecular mechanics energy gradient components.

Includes conversion factors, energy gradient magnitudes for individual
structural objects (atoms, bonds, angles, torsions, outofplanes), system
gradient components, and total energy gradient member data for
mmlib.molecule.Molecule objects.
"""

import numpy
import math

from mmlib import energy
from mmlib import geomcalc
from mmlib import molecule

def _NumDisp():
  """Displacement distance [Angstrom] for numerical gradient."""
  return 1.0 * 10**-6

def KcalAMol2Pa():
  """Conversion from [kcal*A^3/mol] to [Pa] for pressure."""
  return 69476.95

def GetGBond(r_ij, r_eq, k_b):
  """Calculate energy gradient magnitude of bond stretch.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    r_eq (float): Equilibrium bond length [Angstrom] of bond ij.
    k_b (float): Spring constant [kcal/(mol*A^2)] of bond ij.
  
  Returns:
      g_bond (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  g_bond = 2.0 * k_b * (r_ij - r_eq)
  return g_bond

def GetGAngle(a_ijk, a_eq, k_a):
  """Calculate energy gradient magnitude of angle bend.
  
  Args:
    a_ijk (float): Angle [degrees] between atoms i, j, and k.
    a_eq (float): Equilibrium bond angle [degrees] of angle ijk.
    k_a (float): Spring constant [kcal/(mol*rad^2)] of angle ijk.
  
  Returns:
    g_angle (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  g_angle = 2.0 * k_a * (geomcalc.Deg2Rad() * (a_ijk - a_eq) )
  return g_angle

def GetGTorsion(t_ijkl, v_n, gamma, n_fold, paths):
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
  g_torsion = -v_n * n_fold * math.sin(
      geomcalc.Deg2Rad() * (n_fold * t_ijkl - gamma)) / paths
  return g_torsion

def GetGOutofplane(o_ijkl, v_n):
  """Calculate energy gradient magnitude of outofplane bend.
  
  Args:
    o_ijkl (float): Outofplane angle [degrees] between atoms i, j, k, and l.
    v_n (float): Half-barrier height [kcal/mol] of outofplane ijkl.
  
  Returns:
    g_outofplane (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  g_outofplane = (
      -v_n * 2.0 * math.sin(geomcalc.Deg2Rad() * (2.0 * o_ijkl - 180.0)))
  return g_outofplane

def GetGVdwIJ(r_ij, eps_ij, ro_ij):
  """Calculate energy gradient magnitude of van der waals pair energy.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    eps_ij (float): Van der Waals epsilon [kcal/mol] between pair ij.
    ro_ij (float): Van der Waals radius [Angstrom] between pair ij.
  
  Returns:
    g_vdw_ij (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  rrel_ij = ro_ij / r_ij
  g_vdw_ij = -12.0 * (eps_ij / ro_ij) * (rrel_ij**13 - rrel_ij**7)
  return g_vdw_ij

def GetGElstIJ(r_ij, q_i, q_j, epsilon):
  """Calculate energy gradient magnitude of electrostatic pair energy.
  
  Args:
    r_ij (float): Distance [Angstrom] between atoms i and j.
    q_i (float): Partial charge [e] of atom i.
    q_j (float): Partial charge [e] of atom j.
    epsilon (float): Dielectric constant of space (>= 1.0).
  
  Returns:
    e_elst_ij (float): Magnitude of energy gradient [kcal/(mol*A)].
  """
  g_elst_ij = -energy.Ceu2Kcal() * ( q_i * q_j ) / ( epsilon * r_ij**2 )
  return g_elst_ij

def GetGBoundI(k_box, bound, coord, origin, boundtype):
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
  g_bound_i = numpy.zeros(3)
  if (boundtype == 'cube'):
    for j in range(3):
      sign = 1.0 if ((coord[j] - origin[j]) <= 0.0) else -1.0
      scale = 1.0 if (abs(coord[j] - origin[j]) >= bound) else 0.0
      g_bound_i[j] = (-2.0 * sign * scale * k_box * (abs(coord[j]) - bound))
  elif (boundtype == 'sphere'):
    r_io = geomcalc.GetRij(origin, coord)
    u_io = geomcalc.GetUij(origin, coord)
    scale = 1.0 if (r_io >= bound) else 0.0
    g_bound_i = 2.0 * scale * k_box * (r_io - bound) * u_io
  return g_bound_i

def GetGdirInter(coords1, coords2, r_12=None):
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

def GetGdirAngle(coords1, coords2, coords3, r_21=None, r_23=None):
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
  if (not r_21):
    r_21 = geomcalc.GetRij(coords2, coords1)
  if (not r_23):
    r_23 = geomcalc.GetRij(coords2, coords3)
  u_21 = geomcalc.GetUij(coords2, coords1, r_21)
  u_23 = geomcalc.GetUij(coords2, coords3, r_23)
  cp = geomcalc.GetUcp(u_21, u_23)
  gdir1 = geomcalc.GetUcp(u_21, cp) / r_21
  gdir3 = geomcalc.GetUcp(cp, u_23) / r_23
  gdir2 = -1.0 * (gdir1 + gdir3)
  return gdir1, gdir2, gdir3

def GetGdirTorsion(coords1, coords2, coords3, coords4, r_12=None,
                     r_23=None, r_34=None):
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
  if (not r_12):
    r_12 = geomcalc.GetRij(coords1, coords2)
  if (not r_23):
    r_23 = geomcalc.GetRij(coords2, coords3)
  if (not r_34):
    r_34 = geomcalc.GetRij(coords3, coords4)
  u_21 = geomcalc.GetUij(coords2, coords1, r_12)
  u_34 = geomcalc.GetUij(coords3, coords4, r_34)
  u_23 = geomcalc.GetUij(coords2, coords3, r_23)
  u_32 = -1.0 * u_23
  a_123 = geomcalc.GetAijk(coords1, coords2, coords3, r_12, r_23)
  a_432 = geomcalc.GetAijk(coords4, coords3, coords2, r_34, r_23)
  s_123 = math.sin(geomcalc.Deg2Rad() * a_123)
  s_432 = math.sin(geomcalc.Deg2Rad() * a_432)
  c_123 = math.cos(geomcalc.Deg2Rad() * a_123)
  c_432 = math.cos(geomcalc.Deg2Rad() * a_432)
  gdir1 = geomcalc.GetUcp(u_21, u_23) / (r_12*s_123)
  gdir4 = geomcalc.GetUcp(u_34, u_32) / (r_34*s_432)
  gdir2 = (r_12/r_23*c_123 - 1.0)*gdir1 - (r_34/r_23*c_432)*gdir4
  gdir3 = (r_34/r_23*c_432 - 1.0)*gdir4 - (r_12/r_23*c_123)*gdir1
  return gdir1, gdir2, gdir3, gdir4

def GetGdirOutofplane(coords1, coords2, coords3, coords4, oop, r_31=None,
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
  if (not r_31):
    r_31 = geomcalc.GetRij(coords3, coords1)
  if (not r_32):
    r_32 = geomcalc.GetRij(coords3, coords2)
  if (not r_34):
    r_34 = geomcalc.GetRij(coords3, coords4)
  u_31 = geomcalc.GetUij(coords3, coords1, r_31)
  u_32 = geomcalc.GetUij(coords3, coords2, r_32)
  u_34 = geomcalc.GetUij(coords3, coords4, r_34)
  cp_3234 = geomcalc.GetCp(u_32, u_34)
  cp_3431 = geomcalc.GetCp(u_34, u_31)
  cp_3132 = geomcalc.GetCp(u_31, u_32)
  a_132 = geomcalc.GetAijk(coords1, coords3, coords2)
  s_132 = math.sin(geomcalc.Deg2Rad() * a_132)
  c_132 = math.cos(geomcalc.Deg2Rad() * a_132)
  c_oop = math.cos(geomcalc.Deg2Rad() * oop)
  t_oop = math.tan(geomcalc.Deg2Rad() * oop)
  gdir1 = ((1.0/r_31)*(cp_3234/(c_oop*s_132)
      - (t_oop/s_132**2)*(u_31 - c_132*u_32)))
  gdir2 = ((1.0/r_32)*(cp_3431/(c_oop*s_132)
      - (t_oop/s_132**2)*(u_32 - c_132*u_31)))
  gdir4 = ((1.0/r_34)*(cp_3132/(c_oop*s_132)
      - (t_oop*u_34)))
  gdir3 = -1.0*(gdir1 + gdir2 + gdir4)
  return gdir1, gdir2, gdir3, gdir4

def GetGBonds(mol):
  """Calculate bond length energy gradients for all bonds.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with associated Bond objects
        with geometry and parameter data.
  """
  mol.g_bonds.fill(0.0)
  for p in range(mol.n_bonds):
    b = mol.bonds[p]
    c1 = mol.atoms[b.at1].coords
    c2 = mol.atoms[b.at2].coords
    b.grad = GetGBond(b.r_ij, b.r_eq, b.k_b)
    dir1, dir2 = GetGdirInter(c1, c2, b.r_ij)
    mol.g_bonds[b.at1] += b.grad * dir1
    mol.g_bonds[b.at2] += b.grad * dir2

def GetGAngles(mol):
  """Calculate angle bend energy gradients for all angles.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with associated Angle objects
        with geometry and parameter data.
  """
  mol.g_angles.fill(0.0)
  for p in range(mol.n_angles):
    a = mol.angles[p]
    c1 = mol.atoms[a.at1].coords
    c2 = mol.atoms[a.at2].coords
    c3 = mol.atoms[a.at3].coords
    r12 = mol.bond_graph[a.at1][a.at2]
    r23 = mol.bond_graph[a.at2][a.at3]
    a.grad = GetGAngle(a.a_ijk, a.a_eq, a.k_a)
    dir1, dir2, dir3 = GetGdirAngle(c1, c2, c3, r12, r23)
    mol.g_angles[a.at1] += a.grad * dir1
    mol.g_angles[a.at2] += a.grad * dir2
    mol.g_angles[a.at3] += a.grad * dir3

def GetGTorsions(mol):
  """Calculate torsion strain energy gradients for all torsions.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with associated Torsion
        objects with geometry and parameter data.
  """
  mol.g_torsions.fill(0.0)
  for p in range(mol.n_torsions):
    t = mol.torsions[p]
    c1 = mol.atoms[t.at1].coords
    c2 = mol.atoms[t.at2].coords
    c3 = mol.atoms[t.at3].coords
    c4 = mol.atoms[t.at4].coords
    r12 = mol.bond_graph[t.at1][t.at2]
    r23 = mol.bond_graph[t.at2][t.at3]
    r34 = mol.bond_graph[t.at3][t.at4]
    t.grad = GetGTorsion(t.t_ijkl, t.v_n, t.gam, t.n, t.paths)
    dir1, dir2, dir3, dir4 = GetGdirTorsion(c1, c2, c3, c4, r12, r23, r34)
    mol.g_torsions[t.at1] += t.grad * dir1
    mol.g_torsions[t.at2] += t.grad * dir2
    mol.g_torsions[t.at3] += t.grad * dir3
    mol.g_torsions[t.at4] += t.grad * dir4

def GetGOutofplanes(mol):
  """Calculate outofplane bend energy gradients for all outofplanes.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with associated Outofplane
        objects with geometry and parameter data.
  """
  mol.g_outofplanes.fill(0.0)
  for p in range(mol.n_outofplanes):
    o = mol.outofplanes[p]
    c1 = mol.atoms[o.at1].coords
    c2 = mol.atoms[o.at2].coords
    c3 = mol.atoms[o.at3].coords
    c4 = mol.atoms[o.at4].coords
    r31 = mol.bond_graph[o.at3][o.at1]
    r32 = mol.bond_graph[o.at3][o.at2]
    r34 = mol.bond_graph[o.at3][o.at4]
    o.grad = GetGOutofplane(o.o_ijkl, o.v_n)
    dir1, dir2, dir3, dir4 = GetGdirOutofplane(c1, c2, c3, c4, o.o_ijkl, r31,
                                               r32, r34)
    mol.g_outofplanes[o.at1] += o.grad * dir1
    mol.g_outofplanes[o.at2] += o.grad * dir2
    mol.g_outofplanes[o.at3] += o.grad * dir3
    mol.g_outofplanes[o.at4] += o.grad * dir4

def GetGNonbonded(mol):
    """Calculate vdw and elst energy gradients for all nonbonded atom pairs.
    
    Args:
      mol (mmlib.molecule.Molecule): Molecule object with associated Atom
          objects with geometry and parameter data.
    """
    mol.g_nonbonded.fill(0.0)
    mol.g_vdw.fill(0.0)
    mol.g_elst.fill(0.0)
    for i in range(mol.n_atoms):
      at1 = mol.atoms[i]
      for j in range(i+1, mol.n_atoms):
        if (not j in mol.nonints[i]):
          at2 = mol.atoms[j]
          r_ij = geomcalc.GetRij(at1.coords, at2.coords)
          dir1, dir2 = GetGdirInter(at1.coords, at2.coords, r_ij)
          eps_ij = at1.sreps * at2.sreps
          ro_ij = at1.ro + at2.ro
          g_elst = GetGElstIJ(r_ij, at1.charge, at2.charge, mol.dielectric)
          g_vdw = GetGVdwIJ(r_ij, eps_ij, ro_ij)
          mol.g_vdw[i] += g_vdw * dir1
          mol.g_vdw[j] += g_vdw * dir2
          mol.g_elst[i] += g_elst * dir1
          mol.g_elst[j] += g_elst * dir2

def GetGBound(mol):
  """Calculate boundary energy gradients for all atoms.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with boundary parameters and
        Atom objects with geometry data.
  """
  mol.g_bound.fill(0.0)
  k_box = mol.k_box
  bound = mol.bound
  origin = mol.origin
  boundtype = mol.boundtype
  for i in range(mol.n_atoms):
    coords = mol.atoms[i].coords
    mol.g_bound[i] = GetGBoundI(k_box, bound, coords, origin, boundtype)

def GetGTotals(mol):
  """Update total analytic energy gradient [kcal/(mol*A)] of all atoms.
  
  Fundamental components include bonds, angles, torsions, outofplanes, boundary,
  van der waals, and electrostatics. Pre-computed components sum to total.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with energy gradient
        component data [kcal/(mol*A)].
  """
  for i in range(mol.n_atoms):
    for j in range(3):
      mol.g_bonded[i][j] = (mol.g_bonds[i][j] + mol.g_outofplanes[i][j]
           + mol.g_torsions[i][j] + mol.g_angles[i][j])
      mol.g_nonbonded[i][j] = mol.g_vdw[i][j] + mol.g_elst[i][j]
      mol.g_total[i][j] = (
          mol.g_bonded[i][j] + mol.g_nonbonded[i][j] + mol.g_bound[i][j])
  mol.GetPressure()

def _GetVirial(mol):
  """Clausius virial function for all atoms, force,s and coordinates.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with coordinate and force
        data.
  """
  mol.virial = 0.0
  for i in range(mol.n_atoms):
    for j in range(3):
      mol.virial += -mol.atoms[i].coords[j] * mol.g_total[i][j]

def GetPressure(mol):
  """Update total pressure of a system of molecules with boundary.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with temperature, volume, and
        virial data.
  """
  _GetVirial(mol)
  pv = mol.n_atoms * energy.Kb() * mol.temp
  pv += mol.virial / (3.0 * mol.n_atoms)
  mol.press = KcalAMol2Pa() * pv / mol.vol

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
  disp = _NumDisp()
  for i in range(mol.n_atoms):
    for j in range(3):
      q = mol.atoms[i].coords[j]
      qp = q + 0.5*disp
      qm = q - 0.5*disp
      
      mol.atoms[i].coords[j] = qp
      mol.UpdateInternals()
      mol.GetEnergy('standard')
      ep_bond, ep_ang = mol.e_bonds, mol.e_angles
      ep_tor, ep_oop = mol.e_torsions, mol.e_outofplanes
      ep_vdw, ep_elst = mol.e_vdw, mol.e_elst
      ep_bound = mol.e_bound
      
      mol.atoms[i].coords[j] = qm
      mol.UpdateInternals()
      mol.GetEnergy('standard')
      em_bond, em_ang = mol.e_bonds, mol.e_angles
      em_tor, em_oop = mol.e_torsions, mol.e_outofplanes
      em_vdw, em_elst = mol.e_vdw, mol.e_elst
      em_bound = mol.e_bound
      
      mol.atoms[i].coords[j] = q
      mol.g_bonds[i][j] = (ep_bond - em_bond) / disp
      mol.g_angles[i][j] = (ep_ang - em_ang) / disp
      mol.g_torsions[i][j] = (ep_tor - em_tor) / disp
      mol.g_outofplanes[i][j] = (ep_oop - em_oop) / disp
      mol.g_vdw[i][j] = (ep_vdw - em_vdw) / disp
      mol.g_elst[i][j] = (ep_elst - em_elst) / disp
      mol.g_bound[i][j] = (ep_bound - em_bound) / disp
  mol.UpdateInternals()
