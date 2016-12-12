import energy, gradient, molecule, fileio
import math
import numpy as np
import numpy.random

# dynamics.py: classes for handling molecular dynamics data

# gas constant (g * A^2 / (ps^2 * mol * K))
def rgas(): return 0.83144598

# (kcal / A * g) to (A / ps^2)
def acc_conv(): return 25.0 #418.4

# simulation class for molecular dynamics data
class simulation:
    # constructor
    def __init__(self, infile_name):
        self.infile = infile_name
        self.temp = 0.0
        self.tottime = 0.0
        self.time = 0.0
        self.timestep = 0.0
        self.bounds = []
        self.mol = []

        self.geomtime = 0.0
        self.geomout = ''
        self.energytime = 0.0
        self.energyout = ''

        self.read_in_data()

    # read in values from input file
    def read_in_data(self):
        fileio.get_sim_data(self)

    # initialize atomic velocities according to Maxwell-Boltzmann distribution
    def initialize_vels(self):
        if (self.temp > 0.0):
            for i in range(self.mol.n_atoms):
                sigma = math.sqrt(rgas() * self.temp / self.mol.atoms[i].mass)
                for j in range(3):
                    self.mol.atoms[i].vels[j] = np.random.normal(0.0, sigma)
            avgvel = 0.0
            for i in range(self.mol.n_atoms):
                for j in range(3):
                    avgvel += self.mol.atoms[i].vels[j]
            avgvel /= self.mol.n_atoms
            for i in range(self.mol.n_atoms):
                for j in range(3):
                    self.mol.atoms[i].vels[j] += -avgvel

    # molecular dynamics function
    def run_dynamics(self):
        self.open_output_files()
        self.mol.get_energy()
        self.mol.get_gradient('numerical')
        self.update_accs()
        self.initialize_vels()
        self.update_vels(0.5 * self.timestep)
        while (self.time < self.tottime):
            self.mol.print_energy()
            self.update_coords(self.timestep)
            self.update_vels(self.timestep)
            self.mol.get_energy()
            self.mol.get_gradient('numerical')
            self.update_accs()
            self.time += self.timestep
            self.check_print(self.timestep)

    # update accelerations (Angstrom / ps^2)
    def update_accs(self):
        for i in range(self.mol.n_atoms):
            mass = self.mol.atoms[i].mass
            for j in range(3):
                self.mol.atoms[i].accs[j] = -acc_conv() * self.mol.g_total[i][j] / mass

    # update velocities (Angstrom / ps)
    def update_vels(self, tstep):
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].vels[j] += self.mol.atoms[i].accs[j] * tstep

    # update coordinates (Angstrom)
    def update_coords(self, tstep):
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].coords[j] += self.mol.atoms[i].vels[j] * tstep

    # open files for output
    def open_output_files(self):
        self.gtime = 0.0
        self.etime = 0.0
        self.gfile = open(self.geomout, "w")
        self.efile = open(self.energyout, "w")
        self.efile.write('# energy of %s\n' % (self.mol.name))
        self.efile.write('# time e_tot e_kin e_pot e_non e_bon e_vdw e_els e_bon e_ang e_tor e_out\n')

    # check if printing trajectory or energy is necessary
    def check_print(self, timestep):
        self.etime += timestep
        self.gtime += timestep
        if (self.etime >= self.energytime):
            self.print_energy()
            self.etime = 0.0
        if (self.gtime >= self.geomtime):
            self.print_geom()
            self.gtime = 0.0

    # print geometry to trajectory file
    def print_geom(self):
        self.gfile.write('%i\n' % (self.mol.n_atoms))
        self.gfile.write('geometry at t = %.4f ps\n' % (self.time))
        for i in range(self.mol.n_atoms):
            self.gfile.write('%-2s' % (self.mol.atoms[i].element))
            for j in range(3):
                self.gfile.write('  %12.6f' % (self.mol.atoms[i].coords[j]))
            self.gfile.write('\n')

    # print energies to energy file
    def print_energy(self):
        self.efile.write(' %12.4f' % (self.time))
        self.efile.write(' %12.4f' % (self.mol.e_total))
        self.efile.write(' %12.4f' % (self.mol.e_kinetic))
        self.efile.write(' %12.4f' % (self.mol.e_potential))
        self.efile.write(' %12.4f' % (self.mol.e_nonbonded))
        self.efile.write(' %12.4f' % (self.mol.e_bonded))
        self.efile.write(' %12.4f' % (self.mol.e_vdw))
        self.efile.write(' %12.4f' % (self.mol.e_elst))
        self.efile.write(' %12.4f' % (self.mol.e_bonds))
        self.efile.write(' %12.4f' % (self.mol.e_angles))
        self.efile.write(' %12.4f' % (self.mol.e_torsions))
        self.efile.write(' %12.4f' % (self.mol.e_outofplanes))
        self.efile.write('\n')

