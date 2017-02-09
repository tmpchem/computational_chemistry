import math, numpy

# geomcalc.py: functions for calculating molecular coordinate data

# radians to degrees conversion
def rad2deg(): return 180.0 / math.pi
def deg2rad(): return math.pi / 180.0

# calculate square distance between two points
def get_r2_ij(coords_i, coords_j):
    r2_ij = 0.0
    for m in range(3):
        r2_ij += (coords_i[m] - coords_j[m])**2
    return r2_ij

# calculate distance between two points
def get_r_ij(coords_i, coords_j):
    r2_ij = 0.0
    for m in range(3):
        r2_ij += (coords_i[m] - coords_j[m])**2
    r_ij = math.sqrt(r2_ij)
    return r_ij

# calculate unit vector between two points
def get_u_ij(coords_i, coords_j):
    r_ij = get_r_ij(coords_i, coords_j)
    u_ij = numpy.zeros(3)
    if (r_ij > 0.0):
        for m in range(3):
            u_ij[m] = (coords_j[m] - coords_i[m]) / r_ij
    return u_ij

# calculate dot product between two unit vectors
def get_udp(uvec_i, uvec_j):
    udp = 0.0
    for m in range(3):
        udp += uvec_i[m] * uvec_j[m]
    udp = min(udp, 1.0)
    udp = max(udp, -1.0)
    return udp

# calculate unit cross product between two unit vectors
def get_ucp(uvec_i, uvec_j):
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1 - cos_ijk**2)
    if (sin_ijk > 0.0):
        ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1]) / sin_ijk 
        ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2]) / sin_ijk 
        ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0]) / sin_ijk 
    return ucp

# calculate cross product between two unit vectors
def get_cp(uvec_i, uvec_j):
    ucp = numpy.zeros(3)
    cos_ijk = get_udp(uvec_i, uvec_j)
    sin_ijk = math.sqrt(1 - cos_ijk**2)
    if (sin_ijk > 0.0):
        ucp[0] = (uvec_i[1]*uvec_j[2] - uvec_i[2]*uvec_j[1])
        ucp[1] = (uvec_i[2]*uvec_j[0] - uvec_i[0]*uvec_j[2])
        ucp[2] = (uvec_i[0]*uvec_j[1] - uvec_i[1]*uvec_j[0])
    return ucp

# calculate angle between three points
def get_a_ijk(coords_i, coords_j, coords_k):
    u_ji = get_u_ij(coords_j, coords_i)
    u_jk = get_u_ij(coords_j, coords_k)
    dp_jijk = get_udp(u_ji, u_jk)
    a_ijk = rad2deg() * math.acos(dp_jijk)
    return a_ijk

# calculate torsion angle between four points
def get_t_ijkl(coords_i, coords_j, coords_k, coords_l):
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

# calculate out-of-plane (improper torsion) angle between four points
def get_o_ijkl(coords_i, coords_j, coords_k, coords_l):
    u_ki = get_u_ij(coords_k, coords_i)
    u_kj = get_u_ij(coords_k, coords_j)
    u_kl = get_u_ij(coords_k, coords_l)
    u_kikj = get_ucp(u_ki, u_kj)
    dp_kikj_kl = get_udp(u_kikj, u_kl)
    o_ijkl = rad2deg() * math.asin(dp_kikj_kl)
    return o_ijkl

# determine volume of molecular system in Angstroms based on boundary type
def get_volume(mol):
    if (mol.boundtype == 'cube'):
        mol.vol = 8.0 * mol.bound**3
    elif (mol.boundtype == 'sphere'):
        mol.vol = (4.0/3.0) * math.pi * mol.bound**3
    else:
        mol.vol = float('inf')

