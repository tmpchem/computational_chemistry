import math
from mmlib import energy, geomcalc, molecule
import numpy as np

# numerical displacement distance (Angstrom)
num_disp = 1.0 * 10**-6

# conversion from kcal*A^3/mol to joules
kcalamol2j = 69476.95

# calculate bond length energy gradient magnitude between bonded atoms
def get_g_bond(r_ij, r_eq, k_b):
    g_bond = 2.0 * k_b * (r_ij - r_eq)
    return g_bond

# calculate bond angle energy gradient magnitude between bonded atoms
def get_g_angle(a_ijk, a_eq, k_a):
    g_angle = 2.0 * k_a * (geomcalc.deg2rad() * (a_ijk - a_eq) )
    return g_angle

# calculate torsion angle energy gradient magnitude between bonded atoms
def get_g_torsion(t_ijkl, v_n, gamma, n_fold, paths):
    g_torsion = (-v_n * n_fold * math.sin(geomcalc.deg2rad()
        * (n_fold * t_ijkl - gamma)) / paths)
    return g_torsion

# calculate out-of-plane angle (improper torsion) gradient bewteen bonded atoms
def get_g_outofplane(o_ijkl, v_n, gamma, n_fold):
    g_outofplane = (-v_n * n_fold * math.sin(geomcalc.deg2rad()
        * (n_fold * o_ijkl - gamma)))
    return g_outofplane

# calculate van der waals interaction gradient magnitude between atom pair
def get_g_vdw_ij(r_ij, eps_ij, ro_ij):
    rrel_ij = ro_ij / r_ij
    g_vdw_ij = -12.0 * (eps_ij / ro_ij) * (rrel_ij**13 - rrel_ij**7)
    return g_vdw_ij

# calculate electrostatics interaction gradient magnitude between atom pair
def get_g_elst_ij(r_ij, q_i, q_j, epsilon):
    g_elst_ij = -energy.ceu2kcal() * ( q_i * q_j ) / ( epsilon * r_ij**2 )
    return g_elst_ij

# calculate boundary repulsion gradient magnitude for a direction
def get_g_bound_i(k_box, bound, coord, origin, boundtype):
    g_bound_i = np.zeros(3)
    if (boundtype == 'prism'):
        for j in range(3):
            sign = 1.0 if ((coord[j] - origin[j]) <= 0.0) else -1.0
            scale = 1.0 if (abs(coord[j] - origin[j]) >= bound) else 0.0
            g_bound_i[j] = (-2.0 * sign * scale * k_box
                * (abs(coord[j]) - bound))
    elif (boundtype == 'sphere'):
        r_io = geomcalc.get_r_ij(origin, coord)
        u_io = geomcalc.get_u_ij(origin, coord)
        scale = 1.0 if (r_io >= bound) else 0.0
        g_bound_i = 2.0 * scale * k_box * (r_io - bound) * u_io
    elif (bound_type == 'none'):
        pass
    return g_bound_i

# calculate direction of energy gradient between a pair of atoms
def get_gdir_inter(coords1, coords2):
    gdir1 = geomcalc.get_u_ij(coords2, coords1)
    gdir2 = -1.0 * gdir1
    return gdir1, gdir2

# calculate bond angle energy gradient directions between bonded atoms
def get_gdir_angle(coords1, coords2, coords3):
    r_ji = geomcalc.get_r_ij(coords2, coords1)
    r_jk = geomcalc.get_r_ij(coords2, coords3)
    u_ji = geomcalc.get_u_ij(coords2, coords1)
    u_jk = geomcalc.get_u_ij(coords2, coords3)
    cp    = geomcalc.get_ucp(u_ji, u_jk)
    gdir1 = geomcalc.get_ucp(u_ji, cp) / r_ji
    gdir3 = geomcalc.get_ucp(cp, u_jk) / r_jk
    gdir2 = -1.0 * (gdir1 + gdir3)
    return gdir1, gdir2, gdir3

# calculate torsion angle energy gradient directions boetween bonded atoms
def get_gdir_torsion(coords1, coords2, coords3, coords4):
    r_kl = geomcalc.get_r_ij(coords3, coords4)
    r_ij = geomcalc.get_r_ij(coords1, coords2)
    r_jk = geomcalc.get_r_ij(coords2, coords3)
    u_ji = geomcalc.get_u_ij(coords2, coords1)
    u_jk = geomcalc.get_u_ij(coords2, coords3)
    u_kj = geomcalc.get_u_ij(coords3, coords2)
    u_kl = geomcalc.get_u_ij(coords3, coords4)
    a_ijk = geomcalc.get_a_ijk(coords1, coords2, coords3)
    a_lkj = geomcalc.get_a_ijk(coords4, coords3, coords2)
    s_ijk = math.sin(geomcalc.deg2rad() * a_ijk)
    s_lkj = math.sin(geomcalc.deg2rad() * a_lkj)
    c_ijk = math.cos(geomcalc.deg2rad() * a_ijk)
    c_lkj = math.cos(geomcalc.deg2rad() * a_lkj)
    gdir1 = geomcalc.get_ucp(u_ji, u_jk) / (r_ij*s_ijk)
    gdir4 = geomcalc.get_ucp(u_kl, u_kj) / (r_kl*s_lkj)
    gdir2 = gdir1*(r_ij/r_jk*c_ijk - 1.0) - gdir4*(r_kl/r_jk*c_lkj)
    gdir3 = gdir4*(r_kl/r_jk*c_lkj - 1.0) - gdir1*(r_ij/r_jk*c_ijk)
    return gdir1, gdir2, gdir3, gdir4

# calculate outofplane angle energy gradient directions between bonded atoms
def get_gdir_outofplane(coords1, coords2, coords3, coords4, oop):
    r_ki = geomcalc.get_r_ij(coords3, coords1)
    r_kj = geomcalc.get_r_ij(coords3, coords2)
    r_kl = geomcalc.get_r_ij(coords3, coords4)
    u_ki = geomcalc.get_u_ij(coords3, coords1)
    u_kj = geomcalc.get_u_ij(coords3, coords2)
    u_kl = geomcalc.get_u_ij(coords3, coords4)
    cp_kjkl = geomcalc.get_cp(u_kj, u_kl)
    cp_klki = geomcalc.get_cp(u_kl, u_ki)
    cp_kikj = geomcalc.get_cp(u_ki, u_kj)
    a_ikj = geomcalc.get_a_ijk(coords1, coords3, coords2)
    s_ikj = math.sin(geomcalc.deg2rad() * a_ikj)
    c_ikj = math.cos(geomcalc.deg2rad() * a_ikj)
    c_oop = math.cos(geomcalc.deg2rad() * oop)
    t_oop = math.tan(geomcalc.deg2rad() * oop)
    gdir1 = ((1.0/r_ki)*(cp_kjkl/(c_oop*s_ikj)
        - (t_oop/s_ikj**2)*(u_ki - c_ikj*u_kj)))
    gdir2 = ((1.0/r_kj)*(cp_klki/(c_oop*s_ikj)
        - (t_oop/s_ikj**2)*(u_kj - c_ikj*u_ki)))
    gdir4 = ((1.0/r_kl)*(cp_kikj/(c_oop*s_ikj)
        - (t_oop*u_kl)))
    gdir3 = -1.0*(gdir1 + gdir2 + gdir4)
    return gdir1, gdir2, gdir3, gdir4

# update bond lengths and calculate bond length energy gradients
def get_g_bonds(mol):
    mol.g_bonds = np.zeros((mol.n_atoms, 3))
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        c1 = mol.atoms[bond.at1].coords
        c2 = mol.atoms[bond.at2].coords
        bond.r_ij = geomcalc.get_r_ij(c1, c2)
        bond.g = get_g_bond(bond.r_ij, bond.r_eq, bond.k_b)
        dir1, dir2 = get_gdir_inter(c1, c2)
        mol.g_bonds[bond.at1] += bond.g * dir1
        mol.g_bonds[bond.at2] += bond.g * dir2

# update bond angles and calculate bond angle energy gradients
def get_g_angles(mol):
    mol.g_angles = np.zeros((mol.n_atoms, 3))
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        c1 = mol.atoms[ang.at1].coords
        c2 = mol.atoms[ang.at2].coords
        c3 = mol.atoms[ang.at3].coords
        ang.a_ijk = geomcalc.get_a_ijk(c1, c2, c3)
        ang.g = get_g_angle(ang.a_ijk, ang.a_eq, ang.k_a)
        dir1, dir2, dir3 = get_gdir_angle(c1, c2, c3)
        mol.g_angles[ang.at1] += ang.g * dir1
        mol.g_angles[ang.at2] += ang.g * dir2
        mol.g_angles[ang.at3] += ang.g * dir3

# update torsion angles and calculate torsion angle energy gradients
def get_g_torsions(mol):
    mol.g_torsions = np.zeros((mol.n_atoms, 3))
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        c1 = mol.atoms[tor.at1].coords
        c2 = mol.atoms[tor.at2].coords
        c3 = mol.atoms[tor.at3].coords
        c4 = mol.atoms[tor.at4].coords
        tor.t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4)
        tor.g = get_g_torsion(tor.t_ijkl, tor.v_n, tor.gam, tor.n, tor.paths)
        dir1, dir2, dir3, dir4 = get_gdir_torsion(c1, c2, c3, c4)
        mol.g_torsions[tor.at1] += tor.g * dir1
        mol.g_torsions[tor.at2] += tor.g * dir2
        mol.g_torsions[tor.at3] += tor.g * dir3
        mol.g_torsions[tor.at4] += tor.g * dir4

# update outofplane angles and calculate outofplane angle energy gradients
def get_g_outofplanes(mol):
    mol.g_outofplanes = np.zeros((mol.n_atoms, 3))
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        c1 = mol.atoms[oop.at1].coords
        c2 = mol.atoms[oop.at2].coords
        c3 = mol.atoms[oop.at3].coords
        c4 = mol.atoms[oop.at4].coords
        oop.o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4)
        oop.g = get_g_outofplane(oop.o_ijkl, oop.v_n, oop.gam, oop.n)
        dir1, dir2, dir3, dir4 = get_gdir_outofplane(c1,c2,c3,c4,oop.o_ijkl)
        mol.g_outofplanes[oop.at1] += oop.g * dir1
        mol.g_outofplanes[oop.at2] += oop.g * dir2
        mol.g_outofplanes[oop.at3] += oop.g * dir3
        mol.g_outofplanes[oop.at4] += oop.g * dir4

# calculate non-bonded interaction energy gradients between all atoms
def get_g_nonbonded(mol):
    mol.g_nonbonded = np.zeros((mol.n_atoms, 3))
    mol.g_vdw = np.zeros((mol.n_atoms, 3))
    mol.g_elst = np.zeros((mol.n_atoms, 3))
    # van der waals and electrostatic energy
    for i in range(mol.n_atoms):
        atom1 = mol.atoms[i]
        for j in range(i+1, mol.n_atoms):
            if (j in mol.nonints[i]): continue
            atom2 = mol.atoms[j]
            dir1, dir2 = get_gdir_inter(atom1.coords, atom2.coords)
            r_ij = geomcalc.get_r_ij(atom1.coords, atom2.coords)
            eps_ij = atom1.sreps * atom2.sreps
            ro_ij = atom1.ro + atom2.ro
            g_elst = get_g_elst_ij(r_ij, atom1.charge, atom2.charge,
                mol.dielectric)
            g_vdw = get_g_vdw_ij(r_ij, eps_ij, ro_ij)
            mol.g_vdw[i] += g_vdw * dir1
            mol.g_vdw[j] += g_vdw * dir2
            mol.g_elst[i] += g_elst * dir1
            mol.g_elst[j] += g_elst * dir2

# calculate boundary interaction energy gradients for all atoms
def get_g_bound(mol):
    mol.g_bound = np.zeros((mol.n_atoms, 3))
    k_box = mol.k_box
    bound = mol.bound
    origin = mol.origin
    boundtype = mol.boundtype
    for i in range(mol.n_atoms):
        coords = mol.atoms[i].coords
        mol.g_bound[i] = get_g_bound_i(k_box, bound, coords, origin, boundtype)

# update total system analytic energy gradient values
def get_g_totals(mol):
    for i in range(mol.n_atoms):
        for j in range(3):
            mol.g_bonded[i][j]  = mol.g_bonds[i][j] + mol.g_outofplanes[i][j]
            mol.g_bonded[i][j] += mol.g_torsions[i][j] + mol.g_angles[i][j]
            mol.g_nonbonded[i][j] = mol.g_vdw[i][j] + mol.g_elst[i][j]
            mol.g_total[i][j] = mol.g_bonded[i][j] + mol.g_nonbonded[i][j]
            mol.g_total[i][j] += mol.g_bound[i][j]
    mol.get_pressure()

# update clausius virial function
def get_virial(mol):
    mol.virial = 0.0
    for i in range(mol.n_atoms):
        for j in range(3):
            mol.virial += -mol.atoms[i].coords[j] * mol.g_total[i][j]

# update pressure of system
def get_pressure(mol):
    get_virial(mol)
    pv = mol.n_atoms * energy.kb() * mol.temp
    pv += mol.virial / (3.0 * mol.n_atoms)
    mol.press = kcalamol2j * pv / mol.vol

# update total system numerical energy gradient values
def get_g_numerical(mol):
    mol.g_bonds = np.zeros((mol.n_atoms, 3))
    mol.g_angles = np.zeros((mol.n_atoms, 3))
    mol.g_torsions = np.zeros((mol.n_atoms, 3))
    mol.g_outofplanes = np.zeros((mol.n_atoms, 3))
    mol.g_vdw = np.zeros((mol.n_atoms, 3))
    mol.g_elst = np.zeros((mol.n_atoms, 3))
    mol.g_bound = np.zeros((mol.n_atoms, 3))
    for i in range(mol.n_atoms):
        for j in range(3):
            q = mol.atoms[i].coords[j]
            qp = q + 0.5*num_disp
            qm = q - 0.5*num_disp
            mol.atoms[i].coords[j] = qp
            mol.get_energy()
            ep_bond, ep_ang = mol.e_bonds, mol.e_angles
            ep_tor, ep_oop = mol.e_torsions, mol.e_outofplanes
            ep_vdw, ep_elst = mol.e_vdw, mol.e_elst
            ep_bound = mol.e_bound
            mol.atoms[i].coords[j] = qm
            mol.get_energy()
            em_bond, em_ang = mol.e_bonds, mol.e_angles
            em_tor, em_oop = mol.e_torsions, mol.e_outofplanes
            em_vdw, em_elst = mol.e_vdw, mol.e_elst
            em_bound = mol.e_bound
            mol.atoms[i].coords[j] = q
            mol.g_bonds[i][j] = (ep_bond - em_bond) / num_disp
            mol.g_angles[i][j] = (ep_ang - em_ang) / num_disp
            mol.g_torsions[i][j] = (ep_tor - em_tor) / num_disp
            mol.g_outofplanes[i][j] = (ep_oop - em_oop) / num_disp
            mol.g_vdw[i][j] = (ep_vdw - em_vdw) / num_disp
            mol.g_elst[i][j] = (ep_elst - em_elst) / num_disp
            mol.g_bound[i][j] = (ep_bound - em_bound) / num_disp

