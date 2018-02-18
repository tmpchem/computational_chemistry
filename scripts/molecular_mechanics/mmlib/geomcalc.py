"""Functions for computing molecular geometry data.

Includes unit conversions, unit vectors, dot products, cross products,
bond distances, bond angles, torsion angles, outofplane angles, and
system volume.
"""

import numpy
import math

def rad2deg():
    """Conversion from radians to degrees."""
    return 180.0 / math.pi

def deg2rad():
    """Conversion from degrees to radians."""
    return math.pi / 180.0

def get_r2_ij(coords_i, coords_j):
    """Calculate square of distance between two 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
        coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
    
    Returns:
        r2_ij (float): Square distance [Angstrom] between points i and j.
    """
    r2_ij  = (coords_i[0] - coords_j[0])**2
    r2_ij += (coords_i[1] - coords_j[1])**2
    r2_ij += (coords_i[2] - coords_j[2])**2
    return r2_ij

def get_r_ij(coords_i, coords_j):
    """Calculate distance between two 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
        coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
    
    Returns:
        r_ij (float): Distance [Angstrom] between points i and j.
    """
    r2_ij  = (coords_i[0] - coords_j[0])**2
    r2_ij += (coords_i[1] - coords_j[1])**2
    r2_ij += (coords_i[2] - coords_j[2])**2
    r_ij = math.sqrt(r2_ij)
    return r_ij

def get_u_ij(coords_i, coords_j, r_ij=None):
    """Calculate 3d unit vector from cartesian points i to j.
    
    Gives zero vector if i and j are the same point.
    
    Args:
        coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
        coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
        r_ij (float): Distance from i to j [Angstrom], if provided.    

    Returns:
        u_ij (float*): 3 unit vector components from point i to j.
    """
    if (not r_ij):
        r_ij = get_r_ij(coords_i, coords_j)
    u_ij = numpy.zeros(3)
    if (r_ij > 0.0):
        u_ij[0] = (coords_j[0] - coords_i[0]) / r_ij 
        u_ij[1] = (coords_j[1] - coords_i[1]) / r_ij 
        u_ij[2] = (coords_j[2] - coords_i[2]) / r_ij 
    return u_ij

def get_udp(uvec_i, uvec_j):
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
    if (udp > 1.0):
        udp = min(udp, 1.0)
    elif (udp < -1.0):
        udp = max(udp, -1.0)
    return udp

def get_ucp(uvec_i, uvec_j):
    """Calculate unit cross product between two 3d cartesian unit vectors.
    
    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.

    Returns:
        ucp (float*): Normalized cross product between unit vectors i and j.
    """
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1.0 - cos_ijk**2)
    if (sin_ijk > 0.0):
        ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1]) / sin_ijk 
        ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2]) / sin_ijk 
        ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0]) / sin_ijk 
    return ucp

def get_cp(uvec_i, uvec_j):
    """Calculate cross product between two 3d cartesian unit vectors.

    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.

    Returns:
        cp (float*): Cross product between unit vectors i and j.
    """
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1 - cos_ijk**2)
    if (sin_ijk > 0.0):
        ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1])
        ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2])
        ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0])
    return ucp

def get_a_ijk(coords_i, coords_j, coords_k, r_ij=None, r_jk=None):
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
    u_ji = get_u_ij(coords_j, coords_i, r_ij)
    u_jk = get_u_ij(coords_j, coords_k, r_jk)
    dp_jijk = get_udp(u_ji, u_jk)
    a_ijk = rad2deg() * math.acos(dp_jijk)
    return a_ijk

def get_t_ijkl(coords_i, coords_j, coords_k, coords_l, r_ij=None, r_jk=None,
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
    u_ji = get_u_ij(coords_j, coords_i, r_ij)
    u_jk = get_u_ij(coords_j, coords_k, r_jk)
    u_kl = get_u_ij(coords_k, coords_l, r_kl)
    u_jijk =  get_ucp(u_ji, u_jk)
    u_kjkl = -get_ucp(u_jk, u_kl)
    dp_jijk_kjkl = get_udp(u_jijk, u_kjkl)
    sign = 2.0*float(get_udp(u_jijk, u_kl) <= 0.0) - 1.0
    t_ijkl = rad2deg() * sign * math.acos(dp_jijk_kjkl)
    return t_ijkl

def get_o_ijkl(coords_i, coords_j, coords_k, coords_l, r_ki=None, r_kj=None,
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
    u_ki = get_u_ij(coords_k, coords_i, r_ki)
    u_kj = get_u_ij(coords_k, coords_j, r_kj)
    u_kl = get_u_ij(coords_k, coords_l, r_kl)
    u_kikj = get_ucp(u_ki, u_kj)
    dp_kikj_kl = get_udp(u_kikj, u_kl)
    o_ijkl = rad2deg() * math.asin(dp_kikj_kl)
    return o_ijkl

def get_volume(mol):
    """Calculate volume of molecular system based on boundary type
    
    Boundary may be `cube` (V=l**3) or `sphere` (V=4/3pi*l**3). Units
    assumed to be [Angstrom].

    Args:
        mol (mmlib.molecule.Molecule): Molecule object with boundary data.
    """
    if (mol.boundtype == 'cube'):
        mol.vol = 8.0 * mol.bound**3
    elif (mol.boundtype == 'sphere'):
        mol.vol = (4.0/3.0) * math.pi * mol.bound**3
    else:
        mol.vol = float('inf')
