from mmlib import energy, gradient, molecule, fileio
import os, sys, math, time
import numpy as np

# simulation.py: classes for handling molecular dynamics / monte carlo data

# gas constant (g * A^2 / (ps^2 * mol * K))
def rgas(): return 0.83144598

# (kcal / A * g) to (A / ps^2)
def acc_conv(): return 418.4

# simulation class for molecular dynamics data
class simulation:
    # constructor
    def __init__(self, infile_name, sim_type):
        self.infile = os.getcwd() + '/' + infile_name
        self.indir = '/'.join(self.infile.split('/')[:-1])
        self.simtype = sim_type
        self.temp = 0.0
        self.press = 0.0
        self.tottime = 0.0
        self.time = 0.0
        self.timestep = 1.0 * 10**-3
        self.time = 1.0 * 10**-10
        self.mol = []
        self.eqtime = 0.0
        self.eqrate = float('inf')

        self.totconfs = 0
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100

        self.geomtime = 0.01
        self.geomout = 'geom.xyz'
        self.energytime = 0.01
        self.energyout = 'energy.dat'
        self.statustime = float('inf')

        self.read_in_data()

    # read in values from input file
    def read_in_data(self):
        fileio.get_sim_data(self)

    # run simulation based on simulation type
    def run_simulation(self):
        if (self.simtype == 'md'):
            self.run_dynamics()
        elif (self.simtype == 'mmc'):
            self.run_mmc()
        else:
            print('Error: simulation type (%s) not recognized!' % (
                self.simtype))
            sys.exit()

    # initialize atomic velocities according to Maxwell-Boltzmann distribution
    def initialize_vels(self):
        if (self.temp > 0.0):
            self.etemp = self.temp
            self.epress = self.press
            for i in range(self.mol.n_atoms):
                sigma = (math.sqrt(2.0 * rgas() * self.temp
                    / (3.0 * self.mol.atoms[i].mass)))
                for j in range(3):
                    self.mol.atoms[i].vels[j] = np.random.normal(0.0, sigma)
            self.mol.get_energy('standard')
            vscale = math.sqrt(self.temp / self.mol.temp)
            for i in range(self.mol.n_atoms):
                for j in range(3):
                    self.mol.atoms[i].vels[j] *= vscale

    # adjust velocities to equilibrate energy to set temperature
    def equilibrate_temp(self):
        tscale = self.timestep / self.eqrate
        tweight = 10.0 * self.timestep
        self.etemp = (self.etemp + tweight * self.mol.temp)/(1.0 + tweight)
        velscale = 1.0 + tscale * (math.sqrt(self.temp / self.etemp) - 1.0)
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].vels[j] *= velscale

    # generate a random vector to displace 
    def get_rand_disp(self):
        self.rand_disp = np.zeros((self.mol.n_atoms, 3))
        for i in range(self.mol.n_atoms):
            for j in range(3):
                randval = np.random.normal(0.0, self.dispmag)
                self.rand_disp[i][j] = np.random.normal(0.0, self.dispmag)

    # molecular dynamics function
    def run_dynamics(self):
        self.open_output_files()
        self.initialize_vels()
        self.check_print_md(self.timestep)
        self.mol.get_gradient('analytic')
        self.update_accs()
        self.update_vels(0.5*self.timestep)
        while (self.time < self.tottime):
            self.update_coords(self.timestep, 1.0, 0.0)
            self.mol.get_gradient('analytic')
            self.update_accs()
            self.update_vels(self.timestep)
            self.mol.get_energy('leapfrog')
            if (self.time <= self.eqtime):
                self.equilibrate_temp()
            self.check_print_md(self.timestep)
            self.time += self.timestep
        self.check_print_md(self.timestep)
        self.close_output_files()

    # metropolis monte-carlo function
    def run_mmc(self):
        self.open_output_files()
        self.zero_vels()
        self.mol.get_energy('standard')
        penergy = self.mol.e_total
        while (self.conf < self.totconfs):
            self.get_rand_disp()
            self.disp_coords(self.rand_disp)
            self.mol.get_energy('standard')
            delta_e = self.mol.e_total - penergy
            bf = math.exp(min(1.0, -1.0*delta_e / (energy.kb()*self.temp)))
            if (bf >= np.random.random()):
                self.check_print_mc()
                self.conf += 1
                self.n_accept += 1
                penergy = self.mol.e_total
            else:
                self.disp_coords(-1.0*self.rand_disp)
                self.n_reject += 1
            self.check_disp()
        self.close_output_files()

    # update accelerations (Angstrom / ps^2)
    def update_accs(self):
        for i in range(self.mol.n_atoms):
            mass = self.mol.atoms[i].mass
            for j in range(3):
                self.mol.atoms[i].paccs[j] = self.mol.atoms[i].accs[j]
                self.mol.atoms[i].accs[j] = (-acc_conv() *
                    self.mol.g_total[i][j] / mass)

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

    # set velocities to zero
    def zero_vels(self):
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].vels[j] = 0.0

    # displace coordinates by given vector
    def disp_coords(self, disp_vector):
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pcoords[j] = self.mol.atoms[i].coords[j]
                self.mol.atoms[i].coords[j] += disp_vector[i][j]

    # change rms magnitude of displacement vector
    def changedisp(self):
        p_accept = float(self.n_accept) / float(self.n_reject + self.n_accept)
        self.n_accept, self.n_reject = 0, 0
        self.dispmag *= math.exp(2.0 * self.dispinc * (p_accept - 0.5))

    # open files for output
    def open_output_files(self):
        self.gfile = open(self.geomout, "w")
        self.efile = open(self.energyout, "w")
        self.print_energy_header()
        self.print_energy()
        self.print_geom()
        self.print_status()
        self.stime = time.time()
        if (self.simtype == 'md'):
            self.gtime = 10**-10
            self.etime = 10**-10
        elif (self.simtype == 'mmc'):
            self.gconf = 1
            self.econf = 1
            self.dconf = 1

    # close files for output
    def close_output_files(self):
        self.gfile.close()
        self.efile.close()
        self.print_status()

    # check if printing coordinates or energy is necessary (md)
    def check_print_md(self, timestep):
        if (self.etime >= self.energytime):
            self.print_energy()
            self.etime = 10**-10
        if (self.gtime >= self.geomtime):
            self.print_geom()
            self.gtime = 10**-10
        if (time.time() - self.stime > self.statustime):
            self.print_status()
            self.stime = time.time()
        self.etime += timestep
        self.gtime += timestep

    # check if printing coordinates (mc)
    def check_print_mc(self):
        if (self.econf >= self.energyconf):
            self.print_energy()
            self.econf = 0
        if (self.gconf >= self.geomconf):
            self.print_geom()
            self.gconf = 0
        if (time.time() - self.stime > self.statustime):
            self.print_status()
            self.stime = time.time()
        self.econf += 1
        self.gconf += 1

    # check if changing magnitude of random displacement vector
    def check_disp(self):
        if (self.dconf >= self.dispconf):
            self.changedisp()
            self.dconf = 0
        self.dconf += 1

    # print geometry to trajectory file
    def print_geom(self):
        g, m = self.gfile, self.mol
        g.write('%i\n' % (m.n_atoms))
        if (self.simtype == 'md'):
            g.write('geometry at t = %.4f ps\n' % (self.time))
        elif (self.simtype == 'mmc'):
            g.write('geometry at conf %i\n' % (self.conf+1))
        for i in range(m.n_atoms):
            g.write('%-2s' % (m.atoms[i].element))
            for j in range(3):
                g.write('  %12.6f' % (m.atoms[i].coords[j]))
            g.write('\n')

    # print header to energy file
    def print_energy_header(self):
        e = self.efile
        e.write('# energy of %s' % (self.mol.name))
        if (self.simtype == 'md'):
            e.write(' ( %.4f ps of eq)\n         time' % (self.eqtime))
            e.write('    e_total         e_kin       e_pot  ')
        elif (self.simtype == 'mmc'):
            e.write('\n#      conf       e_total')
        e.write('e_nonbond    e_bonded    e_boundary       e_vdw      e_elst')
        e.write('      e_bond     e_angle   e_tors        e_oop')
        if (self.simtype == 'md'):
            e.write(' temperature    pressure\n')

    # write term to file
    def print_val(self, totstr, decstr, val):
        self.efile.write(' %*.*f' % (totstr, decstr, val))

    # write energy sub-terms to file
    def print_e_terms(self, totstr, decstr):
        m = self.mol
        eterms = [m.e_kinetic, m.e_potential, m.e_nonbonded, m.e_bonded,
            m.e_bound, m.e_vdw, m.e_elst, m.e_bonds, m.e_angles,
            m.e_torsions, m.e_outofplanes]
        if (self.simtype == 'mmc'):
            eterms = eterms[2:]
        for i in range(len(eterms)):
            self.print_val(11, 4, eterms[i])

    # print energies to energy file
    def print_energy(self):
        if (self.simtype == 'md'):
            self.print_val(11, 4, self.time)
        elif (self.simtype == 'mmc'):
            self.print_val(9, 0, self.conf+1)
        self.print_val(13, 6, self.mol.e_total)
        self.print_e_terms(11, 4)
        if (self.simtype == 'md'):
            self.print_val(11, 5, self.mol.temp)
            self.print_val(11, 5, self.mol.press)
        self.efile.write('\n')

    # print status of simulation
    def print_status(self):
        if (self.simtype == 'md'):
            print('%.4f/%.4f ps' % (self.time, self.tottime), end='')
        elif (self.simtype == 'mmc'):
            print('%i/%i confs' % (self.conf, self.totconfs), end='')
        print(' as of %s' % (time.strftime('%H:%M:%S')))

