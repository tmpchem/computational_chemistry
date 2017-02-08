
"""Functions for computing molecular mechanics energy components"""

import math
from mmlib import geomcalc

def ceu2kcal():
    """Conversion of electrostatic energy from [ceu] to [kcal/mol]"""
    return 332.06375

def kin2kcal():
    """Conversion of kinetic energy from [amu*A^2/ps^2] to [kcal/mol]"""
    return 0.00239005736

def kb():
    """Boltzmann constant [kcal/(mol*K)]"""
    return 0.001987204

def get_e_bond(r_ij, r_eq, k_b):
    """Calculate bond stretch energy between 2 bonded atoms.
    
    Args:
        r_ij (float): Distance [Angstrom] between atoms i and j.
        r_eq (float): Equilibrium bond length [Angstrom] of bond ij.
        k_b (float): Spring constant [kcal/(mol*A^2)] of bond ij.
        
    Returns:
        e_bond (float): Energy [kcal/mol] of bond ij.
    """
    e_bond = k_b * (r_ij - r_eq)**2
    return e_bond

def get_e_angle(a_ijk, a_eq, k_a):
    """Calculate angle bend energy between 3 bonded atoms.
    
    Args:
        a_ijk (float): Angle [degrees] between atoms i, j, and k.
        a_eq (float): Equilibrium bond angle [degrees] of angle ijk.
        k_a (float): Spring constant [kcal/(mol*rad^2)] of angle ijk.
    
    Returns:
        e_angle (float): Energy [kcal/mol] of angle ijk.
    """
    e_angle = k_a * (geomcalc.deg2rad() * (a_ijk - a_eq) )**2
    return e_angle

def get_e_torsion(t_ijkl, v_n, gamma, n_fold, paths):
    """Calculate torsion strain energy between 4 bonded atoms.
    
    Args:
        t_ijkl (float): Torsion [degrees] between atoms i, j, k, and l.
        v_n (float): Barrier height [kcal/mol] of torsion ijkl.
        gamma (float): Barrier offset [degrees] of torsion ijkl.
        n_fold (int): Barrier frequency of torsion ijkl.
        paths (int): Number of distinct paths in torsion ijkl.
    
    Returns:
        e_torsion (float): Energy [kcal/mol] of torsion ijkl.
    """
    e_torsion = (v_n * (1.0 + math.cos(geomcalc.deg2rad()
        * (n_fold * t_ijkl - gamma))) / paths)
    return e_torsion

def get_e_outofplane(o_ijkl, v_n):
    """Calculate outofplane bend energy between 4 bonded atoms.
    
    Args:
        o_ijkl (float): Outofplane angle [degrees] between atoms
            i, j, k, and l.
        v_n (float): Barrier height [kcal/mol] of torsion ijkl.
    
    Returns:
        e_outofplane (float): Energy [kcal/mol] of outofplane ijkl.
    """
    e_outofplane = (v_n * (1.0 + math.cos(geomcalc.deg2rad()
        * (2.0 * o_ijkl - 180.0))))
    return e_outofplane

def get_e_vdw_ij(r_ij, eps_ij, ro_ij):
    """Calculate van der waals interaction energy between atom pair.
    
    Args:
        r_ij (float): Distance [Angstrom] between atoms i and j.
        eps_ij (float): Van der Waals epsilon [kcal/mol] between pair ij.
        ro_ij (float): Van der Waals radius [Angstrom] between pair ij.
    
    Returns:
        e_vdw_ij (float): Van der waals energy [kcal/mol] between pair ij.
    """
    r6_ij = (ro_ij / r_ij) ** 6
    e_vdw_ij = eps_ij * ( r6_ij**2 - 2.0 * r6_ij )
    return e_vdw_ij

def get_e_elst_ij(r_ij, q_i, q_j, epsilon):
    """Calculate electrostatic interaction energy between atom pair.
    
    Args:
        r_ij (float): Distance [Angstrom] between atoms i and j.
        q_i (float): Partial charge [e] of atom i.
        q_j (float): Partial charge [e] of atom j.
        epsilon (float): Dielectric constant of space (>= 1.0).
    
    Returns:
        e_elst_ij (float): Electrostatic energy [kcal/mol] between pair ij.
    """
    e_elst_ij = ceu2kcal() * ( q_i * q_j ) / ( epsilon * r_ij )
    return e_elst_ij

def get_e_bound_i(k_box, bound, coords, origin, boundtype):
    """Calculate simulation boundary energy of an atom.
    
    Args:
        k_box (float): Spring constant [kcal/(mol*A^2)] of boundary.
        bound (float): Distance from origin [Angstrom] of boundary.
        coords (float*): Array of cartesian coordinates [Angstrom]
            of atom.
        origin (float*): Array of cartesian coordiantes [Angstrom]
            of origin of simulation.
        boundtype (str): `prism` or `sphere`, type of boundary condition.
    
    Returns:
        e_bound_i (float): Boundary energy [kcal/mol] of atom.
    """
    e_bound_i = 0.0
    if (boundtype == 'prism'):
        for j in range(3):
            scale = 1.0 if (abs(coords[j] - origin[j]) >= bound) else 0.0
            e_bound_i += (scale * k_box
                * (abs(coords[j] - origin[j]) - bound)**2)
    elif (boundtype == 'sphere'):
        r_io = geomcalc.get_r_ij(origin, coords)
        u_io = geomcalc.get_u_ij(origin, coords)
        scale = 1.0 if (r_io >= bound) else 0.0
        e_bound_i += scale * k_box * (r_io - bound)**2
    return e_bound_i

def get_e_kinetic_i(mass, vels):
    """Calculate kinetic energy of an atom
    
    Args:
        mass (float): Mass [g/mol] of atom.
        vels (float*): Array of velocities [Angstrom/ps] of atom.
    
    Returns:
        e_kin_i (float): Kinetic energy [kcal/mol] of atom.
    """
    e_kin_i = 0.0
    for i in range(3):
        e_kin_i += kin2kcal() * 0.5 * mass * vels[i]**2
    return e_kin_i

def get_e_nonbonded(mol):
    """Calculate non-bonded interaction energy between all atoms.
    
    Computes van der waals and electrostatic energy [kcal/mol] components
    between all pairs of non-bonded atoms in a system.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Atom objects
            containing cartesian coordinates and molecular mechanics
            parameters
    """
    mol.e_nonbonded, mol.e_vdw, mol.e_elst = 0.0, 0.0, 0.0
    for i in range(mol.n_atoms):
        atom1 = mol.atoms[i]
        for j in range(i+1, mol.n_atoms):
            atom2 = mol.atoms[j]
            if (j in mol.nonints[i]): continue
            r_ij = geomcalc.get_r_ij(atom1.coords, atom2.coords)
            eps_ij = atom1.sreps * atom2.sreps
            ro_ij = atom1.ro + atom2.ro
            mol.e_elst += get_e_elst_ij(r_ij, atom1.charge, atom2.charge,
                mol.dielectric)
            mol.e_vdw += get_e_vdw_ij(r_ij, eps_ij, ro_ij)

def get_e_bonds(mol):
    """Update bond length values and compute bond energy of system.
    
    For all bonds in a molecule, update the distance between atoms ij,
    and add the energy to the total molecular bond energy [kcal/mol].
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Bond objects
            containing atomic indices and molecular mechanics parameters.
    """
    mol.e_bonds = 0.0
    for p in range(mol.n_bonds):
        bond = mol.bonds[p]
        c1 = mol.atoms[bond.at1].coords
        c2 = mol.atoms[bond.at2].coords
        bond.r_ij = geomcalc.get_r_ij(c1, c2)
        bond.e = get_e_bond(bond.r_ij, bond.r_eq, bond.k_b)
        mol.e_bonds += bond.e

def get_e_angles(mol):
    """Update bond angle values and compute angle energy of system.
    
    For all angles in a molecule, update the angle between atoms ijk,
    and add the energy to the total molecular angle energy [kcal/mol].
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Angle objects
            containing atomic indices and molecular mechanics parameters.
    """
    mol.e_angles = 0.0
    for p in range(mol.n_angles):
        ang = mol.angles[p]
        c1 = mol.atoms[ang.at1].coords
        c2 = mol.atoms[ang.at2].coords
        c3 = mol.atoms[ang.at3].coords
        ang.a_ijk = geomcalc.get_a_ijk(c1, c2, c3)
        ang.e = get_e_angle(ang.a_ijk, ang.a_eq, ang.k_a)
        mol.e_angles += ang.e

def get_e_torsions(mol):
    """Update torsion angle values and compute torsion energy of system.
    
    For all torsions in a molecule, update the torsion between atoms ijkl,
    and add the energy to the total molecular torsion energy [kcal/mol].
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Torsion objects
            containing atomic indices and molecular mechanics parameters.
    """
    mol.e_torsions = 0.0
    for p in range(mol.n_torsions):
        tor = mol.torsions[p]
        c1 = mol.atoms[tor.at1].coords
        c2 = mol.atoms[tor.at2].coords
        c3 = mol.atoms[tor.at3].coords
        c4 = mol.atoms[tor.at4].coords
        tor.t_ijkl = geomcalc.get_t_ijkl(c1, c2, c3, c4)
        tor.e = get_e_torsion(tor.t_ijkl, tor.v_n, tor.gam, tor.n, tor.paths)
        mol.e_torsions += tor.e

def get_e_outofplanes(mol):
    """Update outofplane values and compute outofplane energy of system.
    
    For all outofplanes in molecule, update outofplane between atoms ijkl,
    and add the energy to the total molecular outofplane energy [kcal/mol].
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with Outofplane objects
            containing atomic indices and molecular mechanics parameters.
    """
    mol.e_outofplanes = 0.0
    for p in range(mol.n_outofplanes):
        oop = mol.outofplanes[p]
        c1 = mol.atoms[oop.at1].coords
        c2 = mol.atoms[oop.at2].coords
        c3 = mol.atoms[oop.at3].coords
        c4 = mol.atoms[oop.at4].coords
        oop.o_ijkl = geomcalc.get_o_ijkl(c1, c2, c3, c4)
        oop.e = get_e_outofplane(oop.o_ijkl, oop.v_n)
        mol.e_outofplanes += oop.e

def get_e_bound(mol):
    """Compute total boundary energy of system.
    
    For all atoms in molecule, add the boundary energy to the total
    boundary energy [kcal/mol].
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with boundary
            parameters and Atom objects containing cartesian coordinates.
    """
    mol.e_bound = 0.0
    k_box = mol.k_box
    bound = mol.bound
    origin = mol.origin
    boundtype = mol.boundtype
    for i in range(mol.n_atoms):
        atom = mol.atoms[i]
        k_box = mol.k_box
        atom.e_bound = get_e_bound_i(k_box, bound, atom.coords, origin,
            boundtype)
        mol.e_bound += atom.e_bound

def get_temperature(mol):
    """Update kinetic temperature using current kinetic energy per atom.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with current
            (float) kinetic energy [kcal/mol].
    """
    mol.temp = (2.0/3.0) * mol.e_kinetic / (kb() * mol.n_atoms)

def get_e_kinetic(mol, kintype):
    """Compute kinetic energy of all atoms in molecule.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule objects with Atom objects
            containing (float) 3 cartesian velocities [Angstrom/ps].
        kintype (str): Type of kinetic energy to be computed:
            `nokinetic`: Do nothing, ke = 0.0.
            `leapfrog`: Average of current and previous velocities.
            [all others]: Use current velocities.
    """
    mol.e_kinetic = 0.0
    if (kintype == 'nokinetic'):
        pass
    elif (kintype == 'leapfrog'):
        for p in range(mol.n_atoms):
            mass = mol.atoms[p].mass
            vels = 0.5*(mol.atoms[p].vels + mol.atoms[p].pvels)
            e_kin = get_e_kinetic_i(mass, vels)
            mol.e_kinetic += e_kin
    else:
        for p in range(mol.n_atoms):
            mass = mol.atoms[p].mass
            vels = mol.atoms[p].vels
            e_kin = get_e_kinetic_i(mass, vels)
            mol.e_kinetic += e_kin
    mol.get_temperature()

def get_e_totals(mol):
    """Update total energy [kcal/mol] of system from pre-computed components.
    
    Fundamental components include bonds, angles, torsions, outofplanes,
    boundary, van der waals, and electrostatics.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with energy component
            data [kcal/mol].
    """
    mol.e_bonded  = mol.e_bonds + mol.e_angles
    mol.e_bonded += mol.e_torsions + mol.e_outofplanes
    mol.e_nonbonded = mol.e_vdw + mol.e_elst
    mol.e_potential = mol.e_bonded + mol.e_nonbonded + mol.e_bound
    mol.e_total = mol.e_kinetic + mol.e_potential

# end of module

