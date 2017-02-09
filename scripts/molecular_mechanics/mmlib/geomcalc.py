
"""Functions for computing molecular geometry data."""

import math, numpy

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
    r2_ij = 0.0
    for m in range(3):
        r2_ij += (coords_i[m] - coords_j[m])**2
    return r2_ij

def get_r_ij(coords_i, coords_j):
    """Calculate distance between two 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
        coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
    
    Returns:
        r_ij (float): Distance [Angstrom] between points i and j.
    """
    r2_ij = 0.0
    for m in range(3):
        r2_ij += (coords_i[m] - coords_j[m])**2
    r_ij = math.sqrt(r2_ij)
    return r_ij

def get_u_ij(coords_i, coords_j):
    """Calculate unit vector from cartesian points i to j.
    
    Args:
        coords_i (float*): 3 cartesian coordinates [Angstrom] of point i.
        coords_j (float*): 3 cartesian coordinates [Angstrom] of point j.
    
    Returns:
        u_ij (float*): 3 unit vector components from point i to j.
    """
    r_ij = get_r_ij(coords_i, coords_j)
    u_ij = numpy.zeros(3)
    if (r_ij > 0.0):
        for m in range(3):
            u_ij[m] = (coords_j[m] - coords_i[m]) / r_ij
    return u_ij

def get_udp(uvec_i, uvec_j):
    """Calculate dot product between two 3d cartesian unit vectors.
    
    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.

    Returns:
        udp (float): Dot product between unit vectors i and j.
    """
    udp = 0.0
    for m in range(3):
        udp += uvec_i[m] * uvec_j[m]
    udp = min(udp, 1.0)
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
    udp = 0.0
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1 - cos_ijk**2)
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
    udp = 0.0
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1 - cos_ijk**2)
    if (sin_ijk > 0.0):
        ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1])
        ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2])
        ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0])
    return ucp

def get_a_ijk(coords_i, coords_j, coords_k):
    """Calculate angle between 3 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.
        coords_k (float*): 3 cartesian components of unit vector k.

    Returns:
        a_ijk (float): Angle [degrees] between unit vectors ji and jk.
    """
    udp = 0.0
    u_ji = get_u_ij(coords_j, coords_i)
    u_jk = get_u_ij(coords_j, coords_k)
    dp_jijk = get_udp(u_ji, u_jk)
    a_ijk = rad2deg() * math.acos(dp_jijk)
    return a_ijk

def get_t_ijkl(coords_i, coords_j, coords_k, coords_l):
    """Calculate torsion angle between 4 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.
        coords_k (float*): 3 cartesian components of unit vector k.
        coords_l (float*): 3 cartesian components of unit vector l.

    Returns:
        t_ijkl (float): Signed angle [degrees] between planes ijk and jkl.
    """
    udp = 0.0
    u_ji = get_u_ij(coords_j, coords_i)
    u_jk = get_u_ij(coords_j, coords_k)
    u_kj = get_u_ij(coords_k, coords_j)
    u_kl = get_u_ij(coords_k, coords_l)
    u_jijk = get_ucp(u_ji, u_jk)
    u_kjkl = get_ucp(u_kj, u_kl)
    dp_jijk_kjkl = get_udp(u_jijk, u_kjkl)
    sign = 2 * float(get_udp(u_jijk, u_kl) <= 0) - 1
    t_ijkl = rad2deg() * sign * math.acos(dp_jijk_kjkl)
    return t_ijkl

def get_o_ijkl(coords_i, coords_j, coords_k, coords_l):
    """Calculate outofplane angle between 4 3d cartesian points.
    
    Args:
        coords_i (float*): 3 cartesian components of unit vector i.
        coords_j (float*): 3 cartesian components of unit vector j.
        coords_k (float*): 3 cartesian components of unit vector k.
        coords_l (float*): 3 cartesian components of unit vector l.

    Returns:
        o_ijkl (float): Signed angle between plane ijk and vector kl.
    """
    udp = 0.0
    u_ki = get_u_ij(coords_k, coords_i)
    u_kj = get_u_ij(coords_k, coords_j)
    u_kl = get_u_ij(coords_k, coords_l)
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
    udp = 0.0
    if (mol.boundtype == 'cube'):
        mol.vol = 8.0 * mol.bound**3
    elif (mol.boundtype == 'sphere'):
        mol.vol = (4.0/3.0) * math.pi * mol.bound**3
    else:
        mol.vol = float('inf')

# end of module

