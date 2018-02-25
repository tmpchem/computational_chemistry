"""Functions for computing molecular geometry data.

Includes unit conversions, unit vectors, dot products, cross products, bond 
distances, bond angles, torsion angles, outofplane angles, and system volume.
"""

import numpy
import math

from mmlib import constants as const

def GetR2ij(coords_i, coords_j):
  """Calculate square of distance between two 3d cartesian points.
  
  Args:
    coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
    coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
  
  Returns:
    r2_ij (float): Square distance [Angstrom] between points i and j.
  """
  return ((coords_j[0] - coords_i[0])**2 + 
          (coords_j[1] - coords_i[1])**2 +
          (coords_j[2] - coords_i[2])**2)

def GetRij(coords_i, coords_j):
  """Calculate distance between two 3d cartesian points.
  
  Args:
    coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
    coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
  
  Returns:
    r_ij (float): Distance [Angstrom] between points i and j.
  """
  return math.sqrt((coords_j[0] - coords_i[0])**2 +
                   (coords_j[1] - coords_i[1])**2 +
                   (coords_j[2] - coords_i[2])**2)

def GetUij(coords_i, coords_j, r_ij=None):
  """Calculate 3d unit vector from cartesian points i to j.
  
  Gives zero vector if i and j are the same point.
  
  Args:
    coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
    coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
    r_ij (float): Distance from i to j [Angstrom], if provided.    

  Returns:
      u_ij (float*): 3 unit vector components from point i to j.
  """
  if not r_ij:
    r_ij = GetRij(coords_i, coords_j)
  u_ij = numpy.zeros(3)
  if r_ij > 0.0:
    u_ij[0] = (coords_j[0] - coords_i[0]) / r_ij 
    u_ij[1] = (coords_j[1] - coords_i[1]) / r_ij 
    u_ij[2] = (coords_j[2] - coords_i[2]) / r_ij 
  return u_ij

def GetUdp(uvec_i, uvec_j):
  """Calculate dot product between two 3d cartesian unit vectors.
  
  Args:
    coords_i (float*): 3 cartesian components of unit vector i.
    coords_j (float*): 3 cartesian components of unit vector j.

  Returns:
    udp (float): Dot product between unit vectors i and j.
  """
  udp  = uvec_i[0] * uvec_j[0]
  udp += uvec_i[1] * uvec_j[1]
  udp += uvec_i[2] * uvec_j[2]
  return max(-1.0, min(1.0, udp))

def GetUcp(uvec_i, uvec_j):
  """Calculate unit cross product between two 3d cartesian unit vectors.
  
  Args:
    coords_i (float*): 3 cartesian components of unit vector i.
    coords_j (float*): 3 cartesian components of unit vector j.

  Returns:
    ucp (float*): Normalized cross product between unit vectors i and j.
  """
  ucp = numpy.zeros(3)
  cos_ijk = GetUdp(uvec_i, uvec_j)
  sin_ijk = math.sqrt(1.0 - cos_ijk**2)
  if sin_ijk > 0.0:
    ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1]) / sin_ijk 
    ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2]) / sin_ijk 
    ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0]) / sin_ijk 
  return ucp

def GetCp(uvec_i, uvec_j):
  """Calculate cross product between two 3d cartesian unit vectors.

  Args:
    coords_i (float*): 3 cartesian components of unit vector i.
    coords_j (float*): 3 cartesian components of unit vector j.

  Returns:
    cp (float*): Cross product between unit vectors i and j.
  """
  ucp = numpy.zeros(3)
  cos_ijk = GetUdp(uvec_i, uvec_j)
  sin_ijk = math.sqrt(1 - cos_ijk**2)
  if sin_ijk > 0.0:
    ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1])
    ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2])
    ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0])
  return ucp

def GetAijk(coords_i, coords_j, coords_k, r_ij=None, r_jk=None):
  """Calculate angle between 3 3d cartesian points.
  
  Args:
    coords_i (float*): 3 cartesian components of point i.
    coords_j (float*): 3 cartesian components of point j.
    coords_k (float*): 3 cartesian components of point k.
    r_ij (float): Distance between i and j (default None).
    r_jk (float): Distance between j and k (default None).

  Returns:
    a_ijk (float): Angle [degrees] between unit vectors ji and jk.
  """
  u_ji = GetUij(coords_j, coords_i, r_ij)
  u_jk = GetUij(coords_j, coords_k, r_jk)
  dp_jijk = GetUdp(u_ji, u_jk)
  return const.RAD2DEG * math.acos(dp_jijk)

def GetTijkl(coords_i, coords_j, coords_k, coords_l, r_ij=None, r_jk=None,
             r_kl=None):
  """Calculate torsion angle between 4 3d cartesian points.
  
  Args:
    coords_i (float*): 3 cartesian components of point i.
    coords_j (float*): 3 cartesian components of point j.
    coords_k (float*): 3 cartesian components of point k.
    coords_l (float*): 3 cartesian components of point l.
    r_ij (float): Distance between i and j (default None).
    r_jk (float): Distance between j and k (default None).
    r_kl (float): Distance between k and l (default None).
  Returns:
    t_ijkl (float): Signed angle [degrees] between planes ijk and jkl.
  """
  u_ji = GetUij(coords_j, coords_i, r_ij)
  u_jk = GetUij(coords_j, coords_k, r_jk)
  u_kl = GetUij(coords_k, coords_l, r_kl)
  u_jijk =  GetUcp(u_ji, u_jk)
  u_kjkl = -GetUcp(u_jk, u_kl)
  dp_jijk_kjkl = GetUdp(u_jijk, u_kjkl)
  sign = 1.0 if GetUdp(u_jijk, u_kl) <= 0.0 else -1.0
  return const.RAD2DEG * sign * math.acos(dp_jijk_kjkl)

def GetOijkl(coords_i, coords_j, coords_k, coords_l, r_ki=None, r_kj=None,
             r_kl=None):
  """Calculate outofplane angle between 4 3d cartesian points.
  
  Args:
    coords_i (float*): 3 cartesian components of point i.
    coords_j (float*): 3 cartesian components of point j.
    coords_k (float*): 3 cartesian components of point k.
    coords_l (float*): 3 cartesian components of point l.
    r_ki (float): Distance between k and i (default None).
    r_kj (float): Distance between k and j (default None).
    r_kl (float): Distance between k and l (default None).

  Returns:
    o_ijkl (float): Signed angle between plane ijk and vector kl.
  """
  u_ki = GetUij(coords_k, coords_i, r_ki)
  u_kj = GetUij(coords_k, coords_j, r_kj)
  u_kl = GetUij(coords_k, coords_l, r_kl)
  u_kikj = GetUcp(u_ki, u_kj)
  dp_kikj_kl = GetUdp(u_kikj, u_kl)
  return const.RAD2DEG * math.asin(dp_kikj_kl)

def GetVolume(mol):
  """Calculate volume of molecular system based on boundary type
  
  Boundary may be 'cube' (V=l**3) or 'sphere' (V=4/3pi*l**3). Units
  assumed to be [Angstrom].

  Args:
    mol (mmlib.molecule.Molecule): Molecule object with boundary data.
  """
  if mol.boundtype == 'cube':
    mol.vol = 8.0 * mol.bound**3
  elif mol.boundtype == 'sphere':
    mol.vol = 4.0/3.0 * math.pi * mol.bound**3
  else:
    mol.vol = float('inf')
