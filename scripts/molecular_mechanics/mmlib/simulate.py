
"""Classes and functions for handling molecular simulation data."""

import os, sys, math, time, numpy
from mmlib import energy, fileio

def rgas():
    """Gas constant in unit of [amu*A^2/(ps^2*K)]."""
    return 0.83144598

def acc_conv():
    """Conversion of acceleration from [kcal/(A*g) to [A/(ps^2)]."""
    return 418.4

class Simulation:
    """Simulation class for molecular simulation data.
    
    May be used for molecular dynamics or Monte Carlo. Contains
    attributes for handling system propogation and data output.
    
    Many attributes are set by default, but may be overridden by
    keywords in the input file. Mandatory values for input can be
    found in docstring for mmlib.fileio.get_sim_data function.
    
    Args:
        infile_name (str): Path of simulation input file. May be
            relative or absolute, though absolute is safer.
        sim_type (str): Type of simulation.
            `md`: Molecular dynamics.
            `mmc`: Metropolis Monte-Carlo.
    
    Attributes:
        infile (str): Input file (see Args).
        indir (str): Directory of input file
        simtype (str): Type of simulation.
        mol (mmlib.molecule.Molecule): Molecule object from input file.
        temp (float): Desired temperature [K].
        press (float): Desired pressure [bar].
        geomout (str): Geometry printing output file path.
        energyout (str): Energy printing output file path.
        statustime (float): Clock time between printing status to
            standard output [s].

        tottime (float): Total time [ps].
        time (float): Current time [ps].
        timstep (float): Time propogation increment [ps].
        eqtime (float): Total time for thermal equilibration [ps].
        eqrate (float): Rate of thermal equilibration [ps].
        energytime (float): Time between energy printing.
        geomtime (float): Time between geometry printing.

        totconfs (int): Total number of configurations.
        conf (int): Current configuration number.
        dispmag (float): Magnitude of average random displacement [Angstrom].
        dispinc (float): Rate constant for `dispmag` adjustment.
        n_accept (int): Number of accepted MMC trials.
        n_reject (int): Number of rejected MMC trials.
        dispconf (int): Number of configurations between adjusting
            `dispmag` value.
        energyconf (int): Number of configurations between energy printing.
        geomconf (int): Number of configurations between geometry printing.
    """
    def __init__(self, infile_name, sim_type):
        self.infile = os.getcwd() + '/' + infile_name
        self.indir = '/'.join(self.infile.split('/')[:-1])
        self.simtype = sim_type
        self.mol = []
        self.temp = 298.15
        self.press = 1.0
        self.geomout = 'geom.xyz'
        self.energyout = 'energy.dat'
        self.statustime = 60.0

        self.tottime = 1.0
        self.timestep = 1.0 * 10**-3
        self.time = 1.0 * 10**-10
        self.eqtime = 0.0
        self.eqrate = 2.0
        self.energytime = 0.01
        self.geomtime = 0.01

        self.totconfs = 1000
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100

        self.read_in_data()

    def read_in_data(self):
        """Read in simulation data from input file."""
        fileio.get_sim_data(self)
        self.temp += 1.0 * 10**-10

    def run_simulation(self):
        """Run simulation depending on simulation type."""
        if (self.simtype == 'md'):
            self.run_dynamics()
        elif (self.simtype == 'mmc'):
            self.run_mmc()
        else:
            print('Error: simulation type (%s) not recognized!' % (
                self.simtype))
            sys.exit()

    def initialize_vels(self):
        """Initialize atomic velocities depending on temperature.
        
        Selects random velocities [Angstrom/ps] from Gaussian distribution
        centered at zero with sigma according to Maxwell-Boltzmann
        distribution. Then rescales velocities to match specified desired
        temperature.
        """
        if (self.temp > 0.0):
            self.etemp = self.temp
            for i in range(self.mol.n_atoms):
                sigma = (math.sqrt(2.0 * rgas() * self.temp
                    / (3.0 * self.mol.atoms[i].mass)))
                for j in range(3):
                    self.mol.atoms[i].vels[j] = numpy.random.normal(0.0, sigma)
            self.mol.get_energy('standard')
            vscale = math.sqrt(self.temp / self.mol.temp)
            for i in range(self.mol.n_atoms):
                for j in range(3):
                    self.mol.atoms[i].vels[j] *= vscale

    def equilibrate_temp(self):
        """Adjust velocities to equilibrate energy to set temperature.
        
        Computes exponential moving average of kinetic temperature and
        compares it to desired temperature. Velocities are then scaled
        depending on this ratio and the equilibration rate parameters.
        """
        tscale = self.timestep / self.eqrate
        tweight = 10.0 * self.timestep
        self.etemp = (self.etemp + tweight * self.mol.temp)/(1.0 + tweight)
        velscale = 1.0 + tscale * (math.sqrt(self.temp / self.etemp) - 1.0)
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].vels[j] *= velscale

    def get_rand_disp(self):
        """Generate random displacment vector for coordinates.
        
        Random trial displacements for MMC are selected from a Gaussian
        distribution of mu = 0.0 and sigma = `dispmag` attribute for all
        3N atomic coordinates.
        """
        self.rand_disp = numpy.zeros((self.mol.n_atoms, 3))
        for i in range(self.mol.n_atoms):
            for j in range(3):
                randval = numpy.random.normal(0.0, self.dispmag)
                self.rand_disp[i][j] = numpy.random.normal(0.0, self.dispmag)

    def run_dynamics(self):
        """Run molecular dynamics according to simulation parameters.
        
        For every timestep, compute the potential energy and gradient
        of the system. Use velocities to propogate the coordinates and
        forces to propogate velocities in time. Print molecular geometry
        and/or energy data to output files as desired. Run until total
        time reached.
        """
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
            if (self.time < self.eqtime):
                self.equilibrate_temp()
            self.check_print_md(self.timestep)
            self.time += self.timestep
        self.check_print_md(self.timestep)
        self.close_output_files()

    def run_mmc(self):
        """Run Metropolis Monte-Carlo according to simulation parameters.
        
        For every configuration, compute the potential energy and compare
        to the previous step. If the relative Boltzmann factor is above
        a random number, accept or else reject. When desired, alter the
        magnitude of random displacement to seek 50% acceptance. Print
        molecular geometry and/or energy data to output files as desired.
        Run until total configurations reached.
        """
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
            if (bf >= numpy.random.random()):
                self.check_print_mc()
                self.conf += 1
                self.n_accept += 1
                penergy = self.mol.e_total
            else:
                self.disp_coords(-1.0*self.rand_disp)
                self.n_reject += 1
            self.check_disp()
        self.close_output_files()

    def update_accs(self):
        """Update accelerations of atoms [Angstrom/(ps^2)].
        
        Force is the negative gradient of the potential energy. Find
        accelerations by dividing the forces by the atomic masses.
        """
        for i in range(self.mol.n_atoms):
            mass = self.mol.atoms[i].mass
            for j in range(3):
                self.mol.atoms[i].paccs[j] = self.mol.atoms[i].accs[j]
                self.mol.atoms[i].accs[j] = (-acc_conv()
                    * self.mol.g_total[i][j] / mass)

    def update_vels(self, tstep):
        """Update velocities of atoms [Angstrom/ps].
        
        Acceleration is the derivative of velocity with respect to time.
        Find change by multiplying the acceleration by the timestep.
        
        Args:
            tstep (float): time propogation increment [ps].
        """
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pvels[j] = self.mol.atoms[i].vels[j]
                self.mol.atoms[i].vels[j] += self.mol.atoms[i].accs[j] * tstep

    def update_coords(self, tstep, vconst, aconst):
        """Update coordinates of atoms [Angstrom].
        
        Velocity is the derivative of position with respect to time.
        Find displacement by multiplying the velocity by the timestep.
        If desired may propogate using acceleration multiplied by square
        of timestep.
        
        Args:
            tstep (float): time propogation increment [ps].
            vconst (float): fraction of increment to propogate by velocity.
            aconst (float): fraction of increment^2 to propogate by
                acceleration.
        """
        dt = vconst * tstep
        dt2 = aconst * tstep**2
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pcoords[j] = self.mol.atoms[i].coords[j]
                self.mol.atoms[i].coords[j] += self.mol.atoms[i].vels[j] * dt
                self.mol.atoms[i].coords[j] += self.mol.atoms[i].accs[j] * dt2

    def zero_vels(self):
        """Set all 3N atomic velocity components to zero."""
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].vels[j] = 0.0

    def disp_coords(self, disp_vector):
        """Displace all 3N atomic coordinates by specified vector.
        
        Args:
            disp_vector (float**): Nx3 atomic displacement array [Angstrom].
        """
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].pcoords[j] = self.mol.atoms[i].coords[j]
                self.mol.atoms[i].coords[j] += disp_vector[i][j]

    def changedisp(self):
        """Change root-mean-square magnitude of displacement vector.
        
        The MMC random displacement vector has mu = 0.0, and sigma =
        `dispmag` chosen to best approach 50% acceptance ratio. Increase
        `dispmag` when `p_accept` > 0.5 and vice versa.
        """
        p_accept = float(self.n_accept) / float(self.n_reject + self.n_accept)
        self.n_accept, self.n_reject = 0, 0
        self.dispmag *= math.exp(2.0 * self.dispinc * (p_accept - 0.5))

    def open_output_files(self):
        """Open output files for energy and geometry data printing."""
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

    def close_output_files(self):
        """Close output files for energy and geometry data printing."""
        self.print_status()
        self.gfile.close()
        self.efile.close()

    def check_print_md(self, timestep):
        """Check if printing of various md data is needed at current time."""
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

    def check_print_mc(self):
        """Check if printing of various mc data is need at current time."""
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

    def check_disp(self):
        """Check if changing magnitude of random displacment vector needed."""
        if (self.dconf >= self.dispconf):
            self.changedisp()
            self.dconf = 0
        self.dconf += 1

    def print_geom(self):
        """Print xyz-format geometry of system to trajectory file."""
        if (self.simtype == 'md'):
            pstr = 'geometry at t = %.4f ps' % (self.time)
        elif (self.simtype == 'mmc'):
            pstr = 'geometry at conf %i' % (self.conf+1)
        fileio.print_coords_file(self.gfile, self.mol, pstr)

    def print_energy_header(self):
        """Print header of energy output columns to file."""
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

    def print_val(self, totstr, decstr, val):
        """Write specified file to energy output file in indicated format.
        
        Args:
            totstr (int): total number of characters in float print.
            decstr (int): number of post-decimal characters in float print.
            val (float): energy value to be printed to file.
        """
        self.efile.write(' %*.*f' % (totstr, decstr, val))

    def print_e_terms(self, totstr, decstr):
        """Write energy terms at current configuration to energy file.
        
        Args:
            totstr (int): total number of characters in float print.
            decstr (int): number of post-decimal characters in float print.
        """
        m = self.mol
        eterms = [m.e_kinetic, m.e_potential, m.e_nonbonded, m.e_bonded,
            m.e_bound, m.e_vdw, m.e_elst, m.e_bonds, m.e_angles,
            m.e_torsions, m.e_outofplanes]
        if (self.simtype == 'mmc'):
            eterms = eterms[2:]
        for i in range(len(eterms)):
            self.print_val(11, 4, eterms[i])

    def print_energy(self):
        """Print energy data to energy output file, depending on simtype"""
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

    def print_status(self):
        """Print completion progress of simulation to screen"""
        if (self.simtype == 'md'):
            print('%.4f/%.4f ps' % (self.time, self.tottime), end='')
        elif (self.simtype == 'mmc'):
            print('%i/%i confs' % (self.conf, self.totconfs), end='')
        print(' as of %s' % (time.strftime('%H:%M:%S')))
        self.gfile.flush()
        self.efile.flush()
        sys.stdout.flush()

# end of module

