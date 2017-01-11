from mmlib import energy, gradient, molecule, fileio
import math, time
import numpy as np
import numpy.random

# dynamics.py: classes for handling molecular dynamics data

# gas constant (g * A^2 / (ps^2 * mol * K))
def rgas(): return 0.83144598

# boltzmann constant (kcal / (mol * K))
def kb(): return 0.001987204

# (kcal / A * g) to (A / ps^2)
def acc_conv(): return 418.4

# simulation class for molecular dynamics data
class simulation:
    # constructor
    def __init__(self, infile_name):
        self.infile = infile_name
        self.temp = 0.0
        self.tottime = 0.0
        self.time = 0.0
        self.timestep = 0.0
        self.mol = []

        self.geomtime = 0.01
        self.geomout = 'geom.xyz'
        self.energytime = 0.01
        self.energyout = 'energy.dat'
        self.statustime = float('inf')

        self.read_in_data()

    # read in values from input file
    def read_in_data(self):
        fileio.get_sim_data(self)

    # initialize atomic velocities according to Maxwell-Boltzmann distribution
    def initialize_vels(self):
        if (self.temp > 0.0):
            for i in range(self.mol.n_atoms):
                sigma = math.sqrt(2.0 * rgas() * self.temp / (3.0 * self.mol.atoms[i].mass))
                for j in range(3):
                    self.mol.atoms[i].vels[j] = np.random.normal(0.0, sigma)
            tempscale = 0.0
            while (abs(tempscale - 1.0) > 10**-10):
                self.mol.get_energy('standard')
                self.kintemp = self.mol.e_kinetic / kb()
                tempscale = math.sqrt(self.temp / self.kintemp)
                for i in range(self.mol.n_atoms):
                    for j in range(3):
                        self.mol.atoms[i].vels[j] *= tempscale
                vel_cm = np.zeros(3)
                self.mol.mass = 0.0
                for i in range(self.mol.n_atoms):
                    mass = self.mol.atoms[i].mass
                    self.mol.mass += mass
                    for j in range(3):
                        vel_cm[j] += mass * self.mol.atoms[i].vels[j]
                scale = 1.0/self.mol.mass
                for i in range(self.mol.n_atoms):
                    for j in range(3):
                        self.mol.atoms[i].vels[j] += -scale * vel_cm[j]

    # molecular dynamics function
    def run_dynamics(self, int_type):
        self.open_output_files()
        self.initialize_vels()
        # leapfrog integration
        if (int_type == 'leapfrog'):
            self.mol.get_gradient('analytic')
            self.update_accs()
            self.update_vels(0.5*self.timestep)
            while (self.time <= self.tottime):
                self.check_print(self.timestep)
                self.update_coords(self.timestep, 1.0, 0.0)
                self.mol.get_gradient('analytic')
                self.update_accs()
                self.update_vels(self.timestep)
                self.mol.get_energy('leapfrog')
                self.time += self.timestep
        # velocity-verlet integration
        elif (int_type == 'velverlet'):
            self.mol.get_gradient('analytic')
            self.update_accs()
            while (self.time <= self.tottime):
                self.mol.get_energy('standard')
                self.check_print(self.timestep)
                self.update_coords(self.timestep, 1.0, 0.5)
                self.mol.get_gradient('analytic')
                self.update_vels(0.5*self.timestep)
                self.update_accs()
                self.update_vels(0.5*self.timestep)
                self.time += self.timestep
        self.close_output_files()

    # update accelerations (Angstrom / ps^2)
    def update_accs(self):
        for i in range(self.mol.n_atoms):
            mass = self.mol.atoms[i].mass
            for j in range(3):
                self.mol.atoms[i].paccs[j] = self.mol.atoms[i].accs[j]
                self.mol.atoms[i].accs[j] = -acc_conv() * self.mol.g_total[i][j] / mass

    # update velocities (Angstrom / ps)
    def update_vels(self, tstep):
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pvels[j] = self.mol.atoms[i].vels[j]
                self.mol.atoms[i].vels[j] += self.mol.atoms[i].accs[j] * tstep

    # update coordinates (Angstrom)
    def update_coords(self, tstep, vconst, aconst):
        dt = vconst * tstep
        dt2 = aconst * tstep**2
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pcoords[j] = self.mol.atoms[i].coords[j]
                self.mol.atoms[i].coords[j] += self.mol.atoms[i].vels[j] * dt
                self.mol.atoms[i].coords[j] += self.mol.atoms[i].accs[j] * dt2

    # open files for output
    def open_output_files(self):
        self.gtime = 10**-10
        self.etime = 10**-10
        self.stime = time.time()
        self.starttime = self.stime
        self.print_status()
        self.gfile = open(self.geomout, "w")
        self.efile = open(self.energyout, "w")
        self.efile.write('# energy of %s\n' % (self.mol.name))
        self.efile.write('#       time     e_total       e_kin       e_pot   e_nonbond    e_bonded  e_boundary       e_vdw      e_elst      e_bond     e_angle      e_tors       e_oop\n')

    # close files for output
    def close_output_files(self):
        self.gfile.close()
        self.efile.close()

    # check if printing trajectory or energy is necessary
    def check_print(self, timestep):
        if (self.etime >= self.energytime):
            self.print_energy()
            self.etime = 10**-10
        if (self.gtime >= self.geomtime):
            self.print_geom()
            self.gtime = 10**-10
        #if (self.stime >= self.statustime):
        if (time.time() - self.stime > self.statustime):
            self.print_status()
            self.stime = time.time()
        self.etime += timestep
        self.gtime += timestep

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
        self.efile.write(' %11.4f' % (self.time))
        self.efile.write(' %11.4f' % (self.mol.e_total))
        self.efile.write(' %11.4f' % (self.mol.e_kinetic))
        self.efile.write(' %11.4f' % (self.mol.e_potential))
        self.efile.write(' %11.4f' % (self.mol.e_nonbonded))
        self.efile.write(' %11.4f' % (self.mol.e_bonded))
        self.efile.write(' %11.4f' % (self.mol.e_bound))
        self.efile.write(' %11.4f' % (self.mol.e_vdw))
        self.efile.write(' %11.4f' % (self.mol.e_elst))
        self.efile.write(' %11.4f' % (self.mol.e_bonds))
        self.efile.write(' %11.4f' % (self.mol.e_angles))
        self.efile.write(' %11.4f' % (self.mol.e_torsions))
        self.efile.write(' %11.4f' % (self.mol.e_outofplanes))
        self.efile.write('\n')


    # print status of simulation
    def print_status(self):
        print('%11.4f ps as of %s' % (self.time, time.strftime('%H:%M:%S')))

