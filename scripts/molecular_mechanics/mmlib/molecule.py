import sys, math
import numpy as np
import fileio, param, geomcalc, topology, energy, gradient

# molecule.py: classes for handling molecular mechanics data

# atom class for atomic data
class atom:
    # constructor
    def __init__(self, at_type, at_coords, at_charge, at_ro, at_eps, at_mass):
        self.attype = at_type
        self.element = at_type[0].capitalize()
        if (at_type[-1].islower()): self.element += at_type[-1]
        self.coords = at_coords
        self.charge = at_charge
        self.ro = at_ro
        self.eps = at_eps
        self.sreps = math.sqrt(self.eps)
        self.mass = at_mass
        self.covrad = param.get_cov_rad(self.element)
        self.vels = np.zeros(3)
        self.accs = np.zeros(3)
        self.e_nonbonded = 0.0
        self.e_bonded = 0.0
        self.g_nonbonded = np.zeros(3)
        self.g_bonded = np.zeros(3)
    # set new atom type
    def set_attype(self, attype):
        self.attype = attype
    # set new xyz atomic coordinates (Angstrom)
    def set_coords(self, coords):
        self.coords = coords
    # set new single xyz atomic coordinate (Angstrom)
    def set_coord(self, index, coord):
        self.coords[index] = coords
    # set new atomic mass (a.m.u.)
    def set_mass(self, mass):
        self.mass = mass
    # set new atomic charge (electrons)
    def set_charge(self, charge):
        self.charge = charge
    # set new atomic ro van der waals parameter (Angstrom)
    def set_ro(self, ro):
        self.ro = ro
    # set new atomic epsilon van der waals parameter (kcal/mol)
    def set_eps(self, eps):
        self.eps = eps
    # set new atomic element
    def set_element(self, element):
        self.element = element
    # set new atomic covalent radius parameter (Angstrom)
    def set_covrad(self, covrad):
        self.covrad = covrad
    # set new atomic velocities (Angstrom / picosecond)
    def set_vels(self, vels):
        self.vels = vels
    # set new atomic acceleration (Angstrom / picosecond^2)
    def set_accs(self, accs):
        self.accs = accs

# bond class for bond data
class bond:
    # constructor
    def __init__(self, at1, at2, r_ij, r_eq, k_b):
        self.at1 = at1
        self.at2 = at2
        self.r_ij = r_ij
        self.r_eq = r_eq
        self.k_b = k_b
        self.e = 0.0
        self.g = 0.0
    # set new atomic index 1
    def set_at1(self, at1):
        self.at1 = at1
    # set new atomic index 2
    def set_at2(self, at2):
        self.at2 = at2
    # set new bond distance value (Angstrom)
    def set_r_ij(self, r_ij):
        self.r_ij = r_ij
    # set new equilibrium bond distance parameter (Angstrom)
    def set_r_eq(self, r_eq):
        self.r_eq = r_eq
    # set new bond spring constant parameter [kcal/(mol*Angstrom^2)]
    def set_k_b(self, k_b):
        self.k_b = k_b

# bond angle class for bond angle data
class angle:
    # constructor
    def __init__(self, at1, at2, at3, a_ijk, a_eq, k_a):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.a_ijk = a_ijk
        self.a_eq = a_eq
        self.k_a = k_a
        self.e = 0.0
        self.g = 0.0
    # set new atomic index 1
    def set_at1(self, at1):
        self.at1 = at1
    # set new atomic index 2
    def set_at2(self, at2):
        self.at2 = at2
    # set new atomic index 3
    def set_at3(self, at3):
        self.at3 = at3
    # set new bond angle value (degrees)
    def set_a_ijk(self, a_ijk):
        self.a_ijk = a_ijk
    # set new equilibrium bond angle parameter (degrees)
    def set_a_eq(self, a_eq):
        self.a_eq = a_eq
    # set new angle spring constant parameter [kcal/(mol*radian^2)]
    def set_k_a(self, k_a):
        self.k_a = k_a

# torsion class for torsion angle data
class torsion:
    # constructor
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
        self.e = 0.0
        self.g = 0.0
    # set new atomic index 1
    def set_at1(self, at1):
        self.at1 = at1
    # set new atomic index 2
    def set_at2(self, at2):
        self.at2 = at2
    # set new atomic index 3
    def set_at3(self, at3):
        self.at3 = at3
    # set new atomic index 4
    def set_at4(self, at4):
        self.at4 = at4
    # set new torsion angle value (degrees)
    def set_t_ijkl(self, t_ijkl):
        self.t_ijkl = t_ijkl
    # set new torsion magnitude parameter (kcal/mol)
    def set_v_n(self, v_n):
        self.v_n = v_n
    # set new torsion phase factor parameter (degrees)
    def set_gamma(self, gamma):
        self.gam = gamma
    # set new torsion angular frequency parameter (unitless)
    def set_nfold(self, nfold):
        self.nfold = nfold
    # set new torsion equivalent path degeneracy parameter (unitless)
    def set_paths(self, paths):
        self.paths = paths

# outofplane class for outofplane angle data
class outofplane:
    # constructor
    def __init__(self, at1, at2, at3, at4, o_ijkl, v_n, gamma, nfold):
        self.at1 = at1
        self.at2 = at2
        self.at3 = at3
        self.at4 = at4
        self.o_ijkl = o_ijkl
        self.v_n = v_n
        self.gam = gamma
        self.nfold = nfold
        self.e = 0.0
        self.g = 0.0
    # set new atomic index 1
    def set_at1(self, at1):
        self.at1 = at1
    # set new atomic index 2
    def set_at2(self, at2):
        self.at2 = at2
    # set new atomic index 3
    def set_at3(self, at3):
        self.at3 = at3
    # set new atomic index 4
    def set_at4(self, at4):
        self.at4 = at4
    # set new outofplane angle value (degrees)
    def set_o_ijkl(self, o_ijkl):
        self.o_ijkl = o_ijkl
    # set new outofplane magnitude parameter (kcal/mol)
    def set_v_n(self, v_n):
        self.v_n = v_n
    # set new outofplane phase factor parameter (degrees)
    def set_gamma(self, gamma):
        self.gam = gamma
    # set new outofplane angular frequency parameter (unitless)
    def set_nfold(self, nfold):
        self.nfold = nfold

# molecule class for molecular mechanics data
class molecule:
    # constructor
    def __init__(self, infile_name):
        self.infile = infile_name
        self.name = self.infile.split('/')[-1].split('.')[0]
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.torsions = []
        self.outofplanes = []
        self.nonints = []
        
        self.dielectric = 1.0
        
        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_torsions = 0
        self.n_outofplanes = 0
        
        self.e_bonds = 0.0
        self.e_angles = 0.0
        self.e_torsions = 0.0
        self.e_oufofplanes = 0.0
        self.e_bonded = 0.0
        self.e_vdw = 0.0
        self.e_elst = 0.0
        self.e_nonbonded = 0.0
        self.e_potential = 0.0
        self.e_kinetic = 0.0
        self.e_total = 0.0

        self.read_in_data()
        self.get_topology()

        self.g_bonds = np.zeros((self.n_atoms, 3))
        self.g_angles = np.zeros((self.n_atoms, 3))
        self.g_torsions = np.zeros((self.n_atoms, 3))
        self.g_outofplanes = np.zeros((self.n_atoms, 3))
        self.g_bonded = np.zeros((self.n_atoms, 3)) 
        self.g_vdw = np.zeros((self.n_atoms, 3)) 
        self.g_elst = np.zeros((self.n_atoms, 3)) 
        self.g_nonbonded = np.zeros((self.n_atoms, 3)) 
        self.g_total = np.zeros((self.n_atoms, 3)) 

    # read in data from input file
    def read_in_data(self):
        fileio.get_geom(self)

    # determine bonded topology of molecule
    def get_topology(self):
        topology.get_bond_tree(self)
        topology.get_bonds(self)
        topology.get_angles(self)
        topology.get_torsions(self)
        topology.get_outofplanes(self)
        topology.get_nonints(self)

    # calculate energy of molecule
    def get_energy(self):
        energy.get_e_bonds(self)
        energy.get_e_angles(self)
        energy.get_e_torsions(self)
        energy.get_e_outofplanes(self)
        energy.get_e_nonbonded(self)
        energy.get_e_kinetic(self)
        energy.get_e_totals(self)

    # calculate energy gradient of molecule
    def get_gradient(self, grad_type):
        self.grad_type = grad_type
        if (grad_type == 'analytic'):
          self.get_analytic_gradient()
        elif (grad_type == 'numerical'):
          self.get_numerical_gradient()
        else:
          print 'Error: grad type (%s) not recognized!' % (grad_type)
          sys.exit()

    # calculate analytic energy gradient of molecule
    def get_analytic_gradient(self):
        gradient.get_g_bonds(self)
        gradient.get_g_angles(self)
        gradient.get_g_torsions(self)
        gradient.get_g_outofplanes(self)
        gradient.get_g_nonbonded(self)
        gradient.get_g_totals(self)

    # calculate numerical energy gradient of molecule
    def get_numerical_gradient(self):
        gradient.get_g_numerical(self)
        gradient.get_g_totals(self)

    # print energy / topology to screen
    def print_data(self):
        fileio.print_energy(self)
        fileio.print_geom(self)
        fileio.print_bonds(self)
        fileio.print_angles(self)
        fileio.print_torsions(self)
        fileio.print_outofplanes(self)

    # print energy to screen
    def print_energy(self):
        fileio.print_energy(self)

    # print geometry to screen
    def print_geom(self):
        fileio.print_geom(self)

    # print bonds to screen
    def print_bonds(self):
        fileio.print_bonds(self)

    # print angles to screen
    def print_angles(self):
        fileio.print_angles(self)

    # print torsions to screen
    def print_torsions(self):
        fileio.print_torsions(self)

    # print outofplanes to screen
    def print_outofplanes(self):
        fileio.print_outofplanes(self)

    # print gradient to screen
    def print_gradient(self):
        fileio.print_gradient(self, 'total')

