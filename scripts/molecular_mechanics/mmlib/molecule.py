
"""Classes and functions for handling molecular system data."""

import os, sys, math, numpy
from mmlib import fileio, param, geomcalc, topology, energy, gradient

class Atom:
    """Atom class for atomic geometry and parameter data.
    
    Initialize attributes to corresponding specified argument
    values, look up in parameter tables, or set to zero.
    
    Args:
        at_type (str): AMBER94 mm atom type.
        at_coords (float*): 3 cartesian coordinates [Angstrom].
        at_charge (float): atomic partial charge [e].
        at_ro (float): atomic van der waals radius [Angstrom].
        at_eps (float): atomic van der waals epsilon [kcal/mol].
        at_mass (float): atomic mass [g/mol].
    
    Attributes:
        attype (str): AMBER94 mm atom type.
        coords (float*): 3 cartesian coordinates [Angstrom].
        ro (float): vdw radius [Angstrom].
        eps (float): vdw epsilon [kcal/mol].
        sreps (float): square root of vdw epsilon [(kcal/mol)^0.5].
        mass (float): atomic mass [g/mol].
        covard (float): covalent radius [Angstrom].
        vels (float*): 3 cartesian velocity components [Angstrom/ps].
        accs (float*): 3 cartesian acceleration components [A/(ps^2)].
        pcoords (float*): 3 previous `coords` [Angstrom].
        pvels (float*): 3 previous `vels` [Angstrom/ps].
        paccs (float*): 3 previous `accs` [Angstrom/(ps^2)].
    """
    def __init__(self, attype, coords, charge, ro, eps, mass):
        self.attype = attype
        self.element = attype[0].capitalize()
        if (attype[-1].islower()): self.element += attype[-1]
        self.coords = coords
        self.charge = charge
        self.ro = ro
        self.eps = eps
        self.sreps = math.sqrt(self.eps)
        self.mass = mass
        self.covrad = param.get_cov_rad(self.element)
        self.vels = numpy.zeros(3)
        self.accs = numpy.zeros(3)
        self.pcoords = numpy.zeros(3)
        self.pvels = numpy.zeros(3)
        self.paccs = numpy.zeros(3)
    def set_attype(self, attype):
        """Set new (str) atom type."""
        self.attype = attype
    def set_coords(self, coords):
        """Set new (float*) coodinates [Angstrom]."""
        self.coords = coords
    def set_coord(self, index, coord):
        """Set new (float) ith coordinate [Angstrom]."""
        self.coords[index] = coords
    def set_mass(self, mass):
        """Set new (float) atomic mass [g/mol]."""
        self.mass = mass
    def set_charge(self, charge):
        """Set new (float) partial charge [e]."""
        self.charge = charge
    def set_ro(self, ro):
        """Set new (float) vdw radius [Angstrom]."""
        self.ro = ro
    def set_eps(self, eps):
        """Set new (float) vdw epsilon [kcal/mol]."""
        self.eps = eps
    def set_element(self, element):
        """Set new (str) atomic element."""
        self.element = element
    def set_covrad(self, covrad):
        """Set new (float) covalent radius [Angstrom]."""
        self.covrad = covrad
    def set_vels(self, vels):
        """Set new (float*) velocities [Angstrom/ps]."""
        self.vels = vels
    def set_accs(self, accs):
        """Set new (float*) accelerations [Angstrom/(ps^2)]."""
        self.accs = accs

class Bond:
    """Bond class for bond geometry and parameter data.
    
    Initialize attributes to specified argument values. Change
    by calling appropriate `set_[param]` function.
    
    Args / Attributes:
        at1 (str): Atom1 atomic index in Molecule.
        at2 (str): Atom2 atomic index in Molecule.
        r_ij (float): Distance [Angstrom] between at1 and at2.
        r_eq (float): Equlibrium at1-at2 bond length [Angstrom].
        k_b (float): Bond spring constant [kcal/(mol*A^2)].
    Attributes:
        energy (float): Energy of bond [kcal/mol].
        grad (float): Energy gradient magnitude of bond [kcal/(mol*A)].
    """
    def __init__(self, at1, at2, r_ij, r_eq, k_b):
        self.at1 = at1
        self.at2 = at2
        self.r_ij = r_ij
        self.r_eq = r_eq
        self.k_b = k_b

        self.energy = energy.get_e_bond(r_ij, r_eq, k_b)
        self.grad = gradient.get_g_bond(r_ij, r_eq, k_b)

    def set_at1(self, at1):
        """Set new (int) at1."""
        self.at1 = at1
    def set_at2(self, at2):
        """Set new (int) at2."""
        self.at2 = at2
    def set_r_ij(self, r_ij):
        """Set new (float)  [Angstrom]."""
        self.r_ij = r_ij
    def set_r_eq(self, r_eq):
        """Set new (float)  [Angstrom]."""
        self.r_eq = r_eq
    def set_k_b(self, k_b):
        """Set new (float)  [kcal/(mol*A^2)]."""
        self.k_b = k_b

class Angle:
    """Angle class for angle geometry and parameter data.
    
    Initialize attributes to specified argument values. Change
    by calling appropriate `set_[param]` function.
    
    Args / Attributes:
        at1 (str): Atom1 atomic index in Molecule.
        at2 (str): Atom2 atomic index in Molecule.
        at3 (str): Atom3 atomic index in Molecule.
        a_ijk (float): Current at1-at2-at3 bond angle [degrees].
        a_eq (float): Equlibrium at1-at2-at3 bond angle [degrees].
        k_a (float): Angle spring constant [kcal/(mol*rad^2)].
    Attributes:
        energy (float): Energy of angle [kcal/mol].
        grad (float): Energy gradient magnitude of angle [kcal/(mol*A)].
    """
    def __init__(self, at1, at2, at3, a_ijk, a_eq, k_a):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.a_ijk = a_ijk
        self.a_eq = a_eq
        self.k_a = k_a

        self.energy = energy.get_e_angle(a_ijk, a_eq, k_a)
        self.grad = gradient.get_g_angle(a_ijk, a_eq, k_a)

    def set_at1(self, at1):
        """Set new (int) at1."""
        self.at1 = at1
    def set_at2(self, at2):
        """Set new (int) at2."""
        self.at2 = at2
    def set_at3(self, at3):
        """Set new (int) at3."""
        self.at3 = at3
    def set_a_ijk(self, a_ijk):
        """Set new (float) a_ijk [degrees]."""
        self.a_ijk = a_ijk
    def set_a_eq(self, a_eq):
        """Set new (float) a_eq [degrees]."""
        self.a_eq = a_eq
    def set_k_a(self, k_a):
        """Set new (float) k_a [kcal/(mol*rad^2)]."""
        self.k_a = k_a

class Torsion:
    """Torsion class for torsion geometry and parameter data.
    
    Initialize attributes to specified argument values. Change
    by calling appropriate `set_[param]` function.
    
    Args / Attributes:
        at1 (str): Atom1 atomic index in Molecule.
        at2 (str): Atom2 atomic index in Molecule.
        at3 (str): Atom3 atomic index in Molecule.
        at4 (str): Atom4 atomic index in Molecule.
        t_ijkl (float): Current at1-at2-at3-at4 torsion angle [degrees].
        v_n (float): Rotation barrier height [kcal/mol].
        gamma (float): Barrier offset [degrees].
        nfold (int): Frequency of barrier.
        paths (int): Unique paths through torsion.
    Attributes:
        energy (float): Energy of torsion [kcal/mol].
        grad (float): Energy gradient magnitude of torsion [kcal/(mol*A)].
    """
    def __init__(self, at1, at2, at3, at4, t_ijkl, v_n, gamma, nfold, paths):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.at4 = at4
        self.t_ijkl = t_ijkl
        self.v_n = v_n
        self.gam = gamma
        self.n = nfold
        self.paths = paths

        self.energy = energy.get_e_torsion(t_ijkl, v_n, gamma, nfold, paths)
        self.grad = gradient.get_g_torsion(t_ijkl, v_n, gamma, nfold, paths)

    def set_at1(self, at1):
        """Set new (int) at1."""
        self.at1 = at1
    def set_at2(self, at2):
        """Set new (int) at2."""
        self.at2 = at2
    def set_at3(self, at3):
        """Set new (int) at3."""
        self.at3 = at3
    def set_at4(self, at4):
        """Set new (int) at4."""
        self.at4 = at4
    def set_t_ijkl(self, t_ijkl):
        """Set new (float) t_ijkl [degrees]."""
        self.t_ijkl = t_ijkl
    def set_v_n(self, v_n):
        """Set new (float) v_n [kcal/mol]."""
        self.v_n = v_n
    def set_gamma(self, gamma):
        """Set new (float) gamma [degrees]."""
        self.gam = gamma
    def set_nfold(self, nfold):
        """Set new (int) nfold."""
        self.nfold = nfold
    def set_paths(self, paths):
        """Set new (int) paths."""
        self.paths = paths

class Outofplane:
    """Outofplane class for outofplane geometry and parameter data.
    
    Initialize attributes to specified argument values. Change
    by calling appropriate `set_[param]` function.
    
    Args / Attributes:
        at1 (str): Atom1 atomic index in Molecule.
        at2 (str): Atom2 atomic index in Molecule.
        at3 (str): Atom3 atomic index in Molecule.
        at4 (str): Atom4 atomic index in Molecule.
        o_ijkl (float): Current at1-at2-at3-at4 outofplane angle [degrees].
        v_n (float): Barrier height [kcal/mol].
    Attributes:
        energy (float): Energy of outofplane [kcal/mol].
        grad (float): Energy gradient magnitude of outofplane [kcal/(mol*A)].
    """
    def __init__(self, at1, at2, at3, at4, o_ijkl, v_n):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.at4 = at4
        self.o_ijkl = o_ijkl
        self.v_n = v_n
        
        self.energy = energy.get_e_outofplane(o_ijkl, v_n)
        self.grad = gradient.get_g_outofplane(o_ijkl, v_n)

    def set_at1(self, at1):
        """Set new (int) at1."""
        self.at1 = at1
    def set_at2(self, at2):
        """Set new (int) at2."""
        self.at2 = at2
    def set_at3(self, at3):
        """Set new (int) at3."""
        self.at3 = at3
    def set_at4(self, at4):
        """Set new (int) at4."""
        self.at4 = at4
    def set_o_ijkl(self, o_ijkl):
        """Set new (float) o_ijkl [degrees]."""
        self.o_ijkl = o_ijkl
    def set_v_n(self, v_n):
        """Set new (float) v_n [kcal/mol]."""
        self.v_n = v_n

class Molecule:
    """
    Molecule class for molecular geometry / topology / energy data.
    
    Contains `n_atoms` Atom objects, `n_bonds` Bond objects, `n_angles`
    Angle objects, `n_torsions` Torsion objects, and `n_outofplanes`
    Outofplane objects, with all their associated data.
    
    Also contains (float) energy and (float**) gradient total values
    and their components: kinetic, potential, bonded, non-bonded,
    boundary, van der waals, electrostatic, bonds, angles, torsions,
    and outofplanes.
    
    Also contains bulk properties like volume, pressure, temperature,
    boundaries, dielectric, origin and virial.
    
    Args:
        infile_name (str): xyzq or prm input file with molecule data.
    
    Attributes:
        infile (str): Absolute path to `infile_name`
        indir (str): Absolute path to infile directory.
        filetype (str): Input file format: `xyzq` or `prm`.
        name (str): Name of molecule from input file name.

        atoms (mmlib.molecule.Atom*): Array of Atom objects.
        bonds (mmlib.molecule.Bond*): Array of Bond objects.
        angles (mmlib.molecule.Angle*): Array of Atom objects.
        torsions (mmlib.molecule.Torsion*): Array of Torsion objects.
        outofplanes (mmlib.molecule.Outofplanes*): Array of Outofplane  
            objects.
        nonints (int**): Array of covalently bonded atomic indices.

        n_atoms (int): Number of atoms.
        n_bonds (int): Number of bonds.
        n_angles (int): Number of angles.
        n_torsions (int): Number of torsions.
        n_outofplanes (int): Number of outofplanes.

        dielectric (float): Dielectric constant. Default = 1.0 (free space).
        mass (float): Sum of all atomic masses.
        k_box (float): Spring constant [kcal/(mol*A^2)] of boundary potential.
        bound (float): (spherical / cubic) dimensions of system [Angstrom].
        boundtype (str): Type of boundary shape, `cube`, `sphere`, or `none`.
        origin (float*): Cartesian center of boundary object.
        vol (float): Volume of system [A^3].
        temp (float): Instantaneous kinetic temperature [K].
        press (float): Instantaneous kinetic pressure [Pa].
        virial (float): Instantaneous Clausius virial.
       
        e_bonds (float): Bond energy [kcal/mol].
        e_angles (float): Angle energy [kcal/mol].
        e_torsions (float): Torsion energy [kcal/mol].
        e_outofplanes (float): Outofplane energy [kcal/mol].
        e_vdw (float): Van der Waals energy [kcal/mol].
        e_elst (float): Electrostatic energy [kcal/mol].
        e_bound (float): Boundary energy [kcal/mol].
        e_bonded (float): Total bonded energy [kcal/mol].
        e_nonbonded (float): Total nonbonded energy [kcal/mol].
        e_potential (float): Potential energy [kcal/mol].
        e_kinetic (float): Kinetic energy [kcal/mol].
        e_total (float): Total energy [kcal/mol].
        
        g_bonds (float**): Bond energy gradient [kcal/(mol*A)].
        g_angles (float**): Angle energy gradient [kcal/(mol*A)].
        g_torsions (float**): Torsion energy gradient [kcal/(mol*A)].
        g_outofplanes (float**): Outofplane energy gradient [kcal/(mol*A)].
        g_vdw (float**): Van der Waals energy gradient [kcal/(mol*A)].
        g_elst (float**): Electrostatic energy gradient [kcal/(mol*A)].
        g_bound (float**): Boundary energy gradient [kcal/(mol*A)].
        g_bonded (float**): Total bonded energy gradient [kcal/(mol*A)].
        g_nonbonded (float**): Total nonbonded energy gradient [kcal/(mol*A)].
        g_potential (float**): Potential energy gradient [kcal/(mol*A)].
        g_kinetic (float**): Kinetic energy gradient [kcal/(mol*A)].
        g_total (float**): Total energy gradient [kcal/(mol*A)].
    """
    def __init__(self, infile_name):
        self.infile = os.path.realpath(infile_name)
        self.indir = os.path.dirname(self.infile)
        self.filetype = self.infile.split('.')[-1]
        self.name = '.'.join(self.infile.split('/')[-1].split('.')[:-1])
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.outofplanes = []
        self.nonints = []
        
        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_torsions = 0
        self.n_outofplanes = 0

        self.dielectric = 1.0
        self.mass = 0.0
        self.k_box = 0.0
        self.bound = 0.0
        self.boundtype = 'none'
        self.origin = numpy.zeros(3)
        self.vol = float('inf')
        self.temp = 0.0
        self.press = 0.0
        self.virial = 0.0

        self.e_bonds = 0.0
        self.e_angles = 0.0
        self.e_torsions = 0.0
        self.e_outofplanes = 0.0
        self.e_vdw = 0.0
        self.e_elst = 0.0
        self.e_bound = 0.0
        self.e_bonded = 0.0
        self.e_nonbonded = 0.0
        self.e_potential = 0.0
        self.e_kinetic = 0.0
        self.e_total = 0.0

        if (self.filetype == 'xyzq'):
            self.read_in_xyzq()
            self.get_topology()
        elif (self.filetype == 'prm'):
            self.read_in_prm()
            
        self.g_bonds = numpy.zeros((self.n_atoms, 3))
        self.g_angles = numpy.zeros((self.n_atoms, 3))
        self.g_torsions = numpy.zeros((self.n_atoms, 3))
        self.g_outofplanes = numpy.zeros((self.n_atoms, 3))
        self.g_vdw = numpy.zeros((self.n_atoms, 3))
        self.g_elst = numpy.zeros((self.n_atoms, 3))
        self.g_boundary = numpy.zeros((self.n_atoms, 3))
        self.g_bonded = numpy.zeros((self.n_atoms, 3))
        self.g_nonbonded = numpy.zeros((self.n_atoms, 3))
        self.g_total = numpy.zeros((self.n_atoms, 3))

    def read_in_xyzq(self):
        """Read in xyzq data from .xyzq input file"""
        fileio.get_geom(self)

    def read_in_prm(self):
        """Read in prm data from .prm input file"""
        fileio.get_prm(self)

    def get_topology(self):
        """Determine bonded topology of molecules from coordinates"""
        topology.get_bond_graph(self)
        topology.get_bonds(self)
        topology.get_angles(self)
        topology.get_torsions(self)
        topology.get_outofplanes(self)
        topology.get_nonints(self)

    def get_energy(self, kintype):
        """Calculate (float) energy [kcal/mol] and all energy components"""
        energy.get_e_bonds(self)
        energy.get_e_angles(self)
        energy.get_e_torsions(self)
        energy.get_e_outofplanes(self)
        energy.get_e_nonbonded(self)
        energy.get_e_bound(self)
        energy.get_e_kinetic(self, kintype)
        energy.get_e_totals(self)

    def get_gradient(self, grad_type):
        """Calculate analytical or numerical gradient of energy
        
        Args:
            grad_type (str): Type of gradient:
                `analytic`: exact, based on analytic derivatives.
                `numerical`: approximate, based on numerical derivatives.
        """
        self.grad_type = grad_type
        if (grad_type == 'analytic'):
          self.get_analytic_gradient()
        elif (grad_type == 'numerical'):
          self.get_numerical_gradient()
        else:
          print('Error: grad type (%s) not recognized!' % (grad_type))
          sys.exit()

    def get_analytic_gradient(self):
        """Calculate analytic (float**) gradient [kcal/(mol*A)] of energy"""
        gradient.get_g_bonds(self)
        gradient.get_g_angles(self)
        gradient.get_g_torsions(self)
        gradient.get_g_outofplanes(self)
        gradient.get_g_nonbonded(self)
        gradient.get_g_bound(self)
        gradient.get_g_totals(self)

    def get_numerical_gradient(self):
        """Calculate numerical (float**) gradient [kcal/(mol*A)] of energy"""
        gradient.get_g_numerical(self)
        gradient.get_g_totals(self)

    def get_temperature(self):
        """Calculate instantaneous kinetic temperature [K] of system."""
        energy.get_temperature(self)
    def get_pressure(self):
        """Calculate instantaneous kinetic pressure [Pa] of system."""
        gradient.get_pressure(self)
    def get_volume(self):
        """Caclculate approximate volume [A^3] of system."""
        geomcalc.get_volume(self)

    def print_data(self):
        """Print energy / geometry / topology data of molecule to screen."""
        fileio.print_energy(self)
        fileio.print_geom(self)
        fileio.print_bonds(self)
        fileio.print_angles(self)
        fileio.print_torsions(self)
        fileio.print_outofplanes(self)

    def print_energy(self):
        """Print energy and component data of molecule to screen."""
        fileio.print_energy(self)
    def print_geom(self):
        """Print geometry data of molecule to screen"""
        fileio.print_geom(self)
    def print_bonds(self):
        """Print bond data of molecule to screen"""
        fileio.print_bonds(self)
    def print_angles(self):
        """Print angle data of molecule to screen"""
        fileio.print_angles(self)
    def print_torsions(self):
        """Print torsion data of molecule to screen"""
        fileio.print_torsions(self)
    def print_outofplanes(self):
        """Print outofplane data of molecule to screen"""
        fileio.print_outofplanes(self)
    def print_gradient(self):
        """Print gradient data to screen"""
        fileio.print_gradient(self, 'total')

# end of module

