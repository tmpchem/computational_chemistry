
"""Functions for computing molecular mechanics energy gradient components."""

import math, numpy
from mmlib import energy, geomcalc, molecule

def num_disp():
    """Displacement distance [Angstrom] for numerical gradient"""
    return 1.0 * 10**-6

def kcalamol2pa():
    """Conversion from [kcal*A^3/mol] to [Pa] for pressure"""
    return 69476.95

def get_g_bond(r_ij, r_eq, k_b):
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

def get_g_angle(a_ijk, a_eq, k_a):
    """Calculate energy gradient magnitude of angle bend.
    
    Args:
        a_ijk (float): Angle [degrees] between atoms i, j, and k.
        a_eq (float): Equilibrium bond angle [degrees] of angle ijk.
        k_a (float): Spring constant [kcal/(mol*rad^2)] of angle ijk.
    
    Returns:
        g_angle (float): Magnitude of energy gradient [kcal/(mol*A)].
    """
    g_angle = 2.0 * k_a * (geomcalc.deg2rad() * (a_ijk - a_eq) )
    return g_angle

def get_g_torsion(t_ijkl, v_n, gamma, n_fold, paths):
    """Calculate energy gradient magnitude of torsion strain.
    
    Args:
        t_ijkl (float): Torsion [degrees] between atoms i, j, k, and l.
        v_n (float): Barrier height [kcal/mol] of torsion ijkl.
        gamma (float): Barrier offset [degrees] of torsion ijkl.
        n_fold (int): Barrier frequency of torsion ijkl.
        paths (int): Number of distinct paths in torsion ijkl.
    
    Returns:
        g_torsion (float): Magnitude of energy gradient [kcal/(mol*A)].
    """
    g_torsion = (-v_n * n_fold * math.sin(geomcalc.deg2rad()
        * (n_fold * t_ijkl - gamma)) / paths)
    return g_torsion

def get_g_outofplane(o_ijkl, v_n):
    """Calculate energy gradient magnitude of outofplane bend.
    
    Args:
        o_ijkl (float): Outofplane angle [degrees] between atoms
            i, j, k, and l.
        v_n (float): Barrier height [kcal/mol] of torsion ijkl.
    
    Returns:
        g_outofplane (float): Magnitude of energy gradient [kcal/(mol*A)].
    """
    g_outofplane = (-v_n * 2.0 * math.sin(geomcalc.deg2rad()
        * (2.0 * o_ijkl - 180.0)))
    return g_outofplane

def get_g_vdw_ij(r_ij, eps_ij, ro_ij):
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

def get_g_elst_ij(r_ij, q_i, q_j, epsilon):
    """Calculate energy gradient magnitude of electrostatic pair energy.
    
    Args:
        r_ij (float): Distance [Angstrom] between atoms i and j.
        q_i (float): Partial charge [e] of atom i.
        q_j (float): Partial charge [e] of atom j.
        epsilon (float): Dielectric constant of space (>= 1.0).
    
    Returns:
        e_elst_ij (float): Magnitude of energy gradient [kcal/(mol*A)].
    """
    g_elst_ij = -energy.ceu2kcal() * ( q_i * q_j ) / ( epsilon * r_ij**2 )
    return g_elst_ij

def get_g_bound_i(k_box, bound, coord, origin, boundtype):
    """Calculate energy gradient magnitude of boundary energy.
    
    Args:
        k_box (float): Spring constant [kcal/(mol*A^2)] of boundary.
        bound (float): Distance from origin [Angstrom] of boundary.
        coords (float*): Array of cartesian coordinates [Angstrom]
            of atom.
        origin (float*): Array of cartesian coordiantes [Angstrom]
            of origin of simulation.
        boundtype (str): `cube` or `sphere`, type of boundary condition.
    
    Returns:
        g_bound_i (float): Magnitude of energy gradient [kcal/(mol*A)].
    """
    g_bound_i = numpy.zeros(3)
    if (boundtype == 'cube'):
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
    elif (boundtype == 'none'):
        pass
    return g_bound_i

def get_gdir_inter(coords1, coords2):
    """Calculate direction of energy gradient between atom pair.
    
    Args:
        coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
        coords2 (float*): 3 cartesian coordiantes [Angstrom] of atom2.
    
    Returns:
        gdir1 (float*), gdir2 (float*): unit vectors in the direction
            of max increasing inter-atomic distance.
    """
    gdir1 = geomcalc.get_u_ij(coords2, coords1)
    gdir2 = -1.0 * gdir1
    return gdir1, gdir2

def get_gdir_angle(coords1, coords2, coords3):
    """Calculate direction of energy gradients between bond angle atoms.
    
    Args:
        coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
        coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
        coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
    
    Returns:
        gdir1 (float*), gdir2 (float*), gdir3 (float*): vectors in the
            direction of max increasing bond angle.
    """
    r_ji = geomcalc.get_r_ij(coords2, coords1)
    r_jk = geomcalc.get_r_ij(coords2, coords3)
    u_ji = geomcalc.get_u_ij(coords2, coords1)
    u_jk = geomcalc.get_u_ij(coords2, coords3)
    cp    = geomcalc.get_ucp(u_ji, u_jk)
    gdir1 = geomcalc.get_ucp(u_ji, cp) / r_ji
    gdir3 = geomcalc.get_ucp(cp, u_jk) / r_jk
    gdir2 = -1.0 * (gdir1 + gdir3)
    return gdir1, gdir2, gdir3

def get_gdir_torsion(coords1, coords2, coords3, coords4):
    """Calculate direction of energy gradients between torsion atoms.
    
    Args:
        coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
        coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
        coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
        coords4 (float*): 3 cartesian coordinates [Angstrom] of atom4.
    
    Returns:
        gdir1 (float*), gdir2 (float*), gdir3 (float*), gdir4 (float*):
            vectors in the direction of max increasing torsion angle.
    """
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

def get_gdir_outofplane(coords1, coords2, coords3, coords4, oop):
    """Calculate direction of energy gradients between outofplane atoms.
    
    Args:
        coords1 (float*): 3 cartesian coordinates [Angstrom] of atom1.
        coords2 (float*): 3 cartesian coordinates [Angstrom] of atom2.
        coords3 (float*): 3 cartesian coordinates [Angstrom] of atom3.
        coords4 (float*): 3 cartesian coordinates [Angstrom] of atom4.
    
    Returns:
        gdir1 (float*), gdir2 (float*), gdir3 (float*), gdir4 (float*):
            vectors in the direction of max increasing outofplane angle.
    """
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

def get_g_bonds(mol):
    """Calculate bond length energy gradients for all bonds.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with associated
            Bond objects with geometry and parameter data.
    """
    mol.g_bonds = numpy.zeros((mol.n_atoms, 3))
    for p in range(mol.n_bonds):
        b = mol.bonds[p]
        c1 = mol.atoms[b.at1].coords
        c2 = mol.atoms[b.at2].coords
        b.r_ij = geomcalc.get_r_ij(c1, c2)
        b.grad = get_g_bond(b.r_ij, b.r_eq, b.k_b)
        dir1, dir2 = get_gdir_inter(c1, c2)
        mol.g_bonds[b.at1] += b.grad * dir1
        mol.g_bonds[b.at2] += b.grad * dir2

def get_g_angles(mol):
    """Calculate angle bend energy gradients for all angles.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with associated
            Angle objects with geometry and parameter data.
    """
    mol.g_angles = numpy.zeros((mol.n_atoms, 3))
    for p in range(mol.n_angles):
        a = mol.angles[p]
        c1 = mol.atoms[a.at1].coords
        c2 = mol.atoms[a.at2].coords
        c3 = mol.atoms[a.at3].coords
        a.a_ijk = geomcalc.get_a_ijk(c1, c2, c3)
        a.grad = get_g_angle(a.a_ijk, a.a_eq, a.k_a)
        dir1, dir2, dir3 = get_gdir_angle(c1, c2, c3)
        mol.g_angles[a.at1] += a.grad * dir1
        mol.g_angles[a.at2] += a.grad * dir2
        mol.g_angles[a.at3] += a.grad * dir3

def get_g_torsions(mol):
    """Calculate torsion strain energy gradients for all torsions.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with associated
            Torsion objects with geometry and parameter data.
    """
    mol.g_torsions = numpy.zeros((mol.n_atoms, 3))
    for p in range(mol.n_torsions):
        t = mol.torsions[p]
        c1 = mol.atoms[t.at1].coords
        c2 = mol.atoms[t.at2].coords
        c3 = mol.atoms[t.at3].coords
        c4 = mol.atoms[t.at4].coords
        t.t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4)
        t.grad = get_g_torsion(t.t_ijkl, t.v_n, t.gam, t.n, t.paths)
        dir1, dir2, dir3, dir4 = get_gdir_torsion(c1, c2, c3, c4)
        mol.g_torsions[t.at1] += t.grad * dir1
        mol.g_torsions[t.at2] += t.grad * dir2
        mol.g_torsions[t.at3] += t.grad * dir3
        mol.g_torsions[t.at4] += t.grad * dir4

def get_g_outofplanes(mol):
    """Calculate outofplane bend energy gradients for all outofplanes.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with associated
            Outofplane objects with geometry and parameter data.
    """
    mol.g_outofplanes = numpy.zeros((mol.n_atoms, 3))
    for p in range(mol.n_outofplanes):
        o = mol.outofplanes[p]
        c1 = mol.atoms[o.at1].coords
        c2 = mol.atoms[o.at2].coords
        c3 = mol.atoms[o.at3].coords
        c4 = mol.atoms[o.at4].coords
        o.o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4)
        o.grad = get_g_outofplane(o.o_ijkl, o.v_n)
        dir1, dir2, dir3, dir4 = get_gdir_outofplane(c1, c2, c3, c4, o.o_ijkl)
        mol.g_outofplanes[o.at1] += o.grad * dir1
        mol.g_outofplanes[o.at2] += o.grad * dir2
        mol.g_outofplanes[o.at3] += o.grad * dir3
        mol.g_outofplanes[o.at4] += o.grad * dir4

def get_g_nonbonded(mol):
    """Calculate vdw and elst energy gradients for all nonbonded atom pairs.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with associated
            Atom objects with geometry and parameter data.
    """
    mol.g_nonbonded = numpy.zeros((mol.n_atoms, 3))
    mol.g_vdw = numpy.zeros((mol.n_atoms, 3))
    mol.g_elst = numpy.zeros((mol.n_atoms, 3))
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

def get_g_bound(mol):
    """Calculate boundary energy gradients for all atoms.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with boundary
            parameters and Atom objects with geometry data.
    """
    mol.g_bound = numpy.zeros((mol.n_atoms, 3))
    k_box = mol.k_box
    bound = mol.bound
    origin = mol.origin
    boundtype = mol.boundtype
    for i in range(mol.n_atoms):
        coords = mol.atoms[i].coords
        mol.g_bound[i] = get_g_bound_i(k_box, bound, coords, origin, boundtype)

def get_g_totals(mol):
    """Update total analytic energy gradient [kcal/(mol*A)] of all atoms.
    
    Fundamental components include bonds, angles, torsions, outofplanes,
    boundary, van der waals, and electrostatics. Pre-computed components
    sum to total.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with energy
            gradient component data [kcal/(mol*A)].
    """
    for i in range(mol.n_atoms):
        for j in range(3):
            mol.g_bonded[i][j]  = mol.g_bonds[i][j] + mol.g_outofplanes[i][j]
            mol.g_bonded[i][j] += mol.g_torsions[i][j] + mol.g_angles[i][j]
            mol.g_nonbonded[i][j] = mol.g_vdw[i][j] + mol.g_elst[i][j]
            mol.g_total[i][j] = mol.g_bonded[i][j] + mol.g_nonbonded[i][j]
            mol.g_total[i][j] += mol.g_bound[i][j]
    mol.get_pressure()

def get_virial(mol):
    """Clausius virial function for all atoms, force,s and coordinates.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with coordinate
            and force data.
    """
    mol.virial = 0.0
    for i in range(mol.n_atoms):
        for j in range(3):
            mol.virial += -mol.atoms[i].coords[j] * mol.g_total[i][j]

def get_pressure(mol):
    """Update total pressure of a system of molecules with boundary.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with temperature,
            volume, and virial data.
    """
    get_virial(mol)
    pv = mol.n_atoms * energy.kb() * mol.temp
    pv += mol.virial / (3.0 * mol.n_atoms)
    mol.press = kcalamol2pa() * pv / mol.vol

def get_g_numerical(mol):
    """Update total numerical energy gradient [kcal/(mol*A)] of all atoms.

    Args:
        mol (mmlib.molecule.Molecule): Molecule object with energy
            gradient component data [kcal/(mol*A)].
    """
    mol.g_bonds = numpy.zeros((mol.n_atoms, 3))
    mol.g_angles = numpy.zeros((mol.n_atoms, 3))
    mol.g_torsions = numpy.zeros((mol.n_atoms, 3))
    mol.g_outofplanes = numpy.zeros((mol.n_atoms, 3))
    mol.g_vdw = numpy.zeros((mol.n_atoms, 3))
    mol.g_elst = numpy.zeros((mol.n_atoms, 3))
    mol.g_bound = numpy.zeros((mol.n_atoms, 3))
    num_disp = num_disp()
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

# end of module

