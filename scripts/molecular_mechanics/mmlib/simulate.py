"""Classes and functions for handling molecular simulation data.

Includes unit conversions and classes to do molecular dynamics and Metropolis
Monte-Carlo simulations of mmlib.molecule.Molecule objects.
"""

import math
import numpy
import os
import sys
import time

from mmlib import constants as const
from mmlib import fileio

class Simulation:
  """Simulation base class for molecular simulation data.

  Derives into MolecularDynamics and MonteCarlo. Contains attributes for
  handling system propogation and data output.

  Many attributes are set by default, but may be overridden by keywords in the
  input file. Mandatory values for input can be found in docstring for
  mmlib.fileio.get_sim_data function.
  
  Args:
    infile_name (str): Path of simulation input file. May be relative or
        absolute, though absolute is safer.
  
  Attributes:
    infile (str): Absolute input file path (see Args).
    indir (str): Directory of input file
    simtype (str): Type of simulation (set through constructors of derived
        classes). Values:
          'md': Molecular dynamics.
          'mc': Metropolis Monte-Carlo.
    mol (mmlib.molecule.Molecule): Molecule object from input file.
    temperature (float): Desired temperature [K].
    pressure (float): Desired pressure [bar].
    geomout (str): Geometry printing output file path.
    energyout (str): Energy printing output file path.
    statustime (float): Clock time between printing status to standard output
        [s].
    random_seed (int): Random number generator seed for initial velocity
        assignment / displacement (random uint32 by default).
    eprintdig (int): Post-decimal digits for energy output.
    eprintchar (int): Total characters for energy output.
    gprintdig (int): Post-decimal digits for geometry output.
    gprintchar (int): Total characters for geometry output.
  """
  def __init__(self, infile_name):
    self.infile = os.path.realpath(infile_name)
    self.indir = os.path.dirname(self.infile)
    self.mol = []
    self.temperature = 298.15
    self.pressure = 1.0
    self.geomout = 'geom.xyz'
    self.energyout = 'energy.dat'
    self.statustime = 60.0
    self.random_seed = numpy.random.randint(2**31)
    self.eprintdig = 3
    self.eprintchar = 10
    self.gprintdig = 3
    self.gprintchar = 7

    self.ReadInData()

  def ReadInData(self):
    """Read in simulation data from input file."""
    fileio.GetSimData(self)
    # zero divison error workaround
    self.temperature += 1.0E-20

  def _OpenOutputFiles(self):
    """Open output files for energy and geometry data printing."""
    self.gfile = open(self.geomout, "w")
    self.efile = open(self.energyout, "w")
    self._PrintEnergyHeader()
    self.stime = time.time()
    if self.simtype == 'md':
      self.gtime = 10E-10
      self.etime = 10E-10
    elif self.simtype == 'mc':
      self.gconf = 0
      self.econf = 0
      self.dconf = 0

  def _CloseOutputFiles(self):
    """Close output files for energy and geometry data printing."""
    self._PrintStatus()
    self.gfile.close()
    self.efile.close()

  def _FlushBuffers(self):
    """Flush buffers to output files and screen output."""
    self.gfile.flush()
    self.efile.flush()
    sys.stdout.flush()

  def _PrintGeom(self):
    """Print xyz-format geometry of system to trajectory file."""
    if self.simtype == 'md':
      comment = '%.4f ps' % (self.time)
    elif self.simtype == 'mc':
      comment = 'conf %i' % (self.conf)

    self.gfile.write(fileio.GetPrintCoordsXyzString(
        self.mol.atoms, comment, self.gprintchar, self.gprintdig))

  def _PrintVal(self, totstr, decstr, val, ptype='f', n_space=1):
    """Write specified file to energy output file in indicated format.
    
    Args:
      totstr (int): Total number of characters in float print.
      decstr (int): Number of post-decimal characters in float print.
      val (float): Energy value [kcal/mol] to be printed to file.
      ptype (char): Type of number to print to output file:
        * 'f': Printf floating point.
        * 'e': Printf exponential.
      n_space (int): Leading number of spaces before printing value.
    """
    if ptype == 'f':
      self.efile.write('%*s%*.*f' % (n_space, '', totstr, decstr, val))
    elif ptype == 'e':
      self.efile.write('%*s%*.*e' % (n_space, '', totstr, decstr, val))

  def _PrintETerms(self, totstr, decstr, ptype):
    """Write energy terms at current configuration to energy file.
    
    Args:
      totstr (int): total number of characters in float print.
      decstr (int): number of post-decimal characters in float print.
    """
    m = self.mol
    eterms = [
        m.e_kinetic, m.e_potential, m.e_nonbonded, m.e_bonded, m.e_bound,
        m.e_vdw, m.e_elst, m.e_bonds, m.e_angles, m.e_torsions, m.e_outofplanes]
    if self.simtype == 'mc':
      eterms = eterms[2:]
    for i in range(len(eterms)):
      self._PrintVal(totstr, decstr, eterms[i], ptype)

  def _PrintEnergy(self):
    """Print energy data to energy output file, depending on simtype"""
    if self.simtype == 'md':
      self._PrintVal(self.tprintchar, self.tprintdig, self.time, 'f', 0)
    elif self.simtype == 'mc':
      self._PrintVal(self.cprintchar, 0, self.conf, 'f', 0)
    self._PrintVal(self.eprintchar+2, self.eprintdig+2, self.mol.e_total, 'e')
    self._PrintETerms(self.eprintchar, self.eprintdig, 'e')
    self.efile.write('\n')

  def _PrintStatus(self):
    """Print completion progress of simulation to screen"""
    if self.simtype == 'md':
      print('%.*f/%.*f ps' % (self.tprintdig, self.time, self.tprintdig, 
          self.tottime), end='')
    elif self.simtype == 'mc':
        print('%i/%i confs' % (self.conf, self.totconf), end='')
    print(' as of %s' % (time.strftime('%H:%M:%S')))
    self._FlushBuffers()


class MolecularDynamics(Simulation):
  """Molecular dynamics class derived from Simulation base class.
  
  Molecular dynamics takes an initial geometry, randomly assigns velocities, and
  propogates system forward in time for a specified amount of time through
  numerical integration of Newton's equations.
  
  Args:
    infile_name (str): Path of simulation input file. May be relative or
        absolute, though absolute is safer.
      
  Attributes:
    tottime (float): Total time [ps].
    time (float): Current time [ps].
    timstep (float): Time propogation increment [ps].
    eqtime (float): Total time for thermal equilibration [ps].
    eqrate (float): Rate of thermal equilibration [ps].
    energytime (float): Time between energy printing.
    geomtime (float): Time between geometry printing.
    tprintdig (int): Post-decimal digits for time printing.
    tprintchar (int): Total characters for time printing.
  """
  def __init__(self, infile_name):
    self.simtype = 'md'

    self.tottime = 1.0
    self.timestep = 1.0E-3
    self.time = 1.0E-10
    self.eqtime = 0.0
    self.eqrate = 2.0
    self.energytime = 1.0E-2
    self.geomtime = 1.0E-2

    Simulation.__init__(self, infile_name)

    self.tprintdig = 4
    self.tprintchar = 7

  def Run(self):
    """Run molecular dynamics according to simulation parameters.
    
    For every timestep, compute the potential energy and gradient of the system.
    Use velocities to propogate the coordinates and forces to propogate velocities
    in time. Print molecular geometry and/or energy data to output files as
    desired. Run until total time reached.
    """
    self._OpenOutputFiles()
    self._InitializeVels()
    self.mol.GetEnergy()
    self.mol.GetGradient()
    self._UpdateAccs()
    self._CheckPrint(0.0, print_all=True)
    self._UpdateVels(0.5*self.timestep)
    while self.time < self.tottime:
      self._UpdateCoords(self.timestep)
      self.mol.GetGradient()
      self._UpdateAccs()
      self._UpdateVels(self.timestep)
      self.mol.GetEnergy('leapfrog')
      if self.time < self.eqtime:
        self._EquilibrateTemp()
      self._CheckPrint(self.timestep)
      self.time += self.timestep
    self._CheckPrint(self.timestep)
    self._CloseOutputFiles()

  def _InitializeVels(self):
    """Initialize atomic velocities depending on temperature.
    
    Selects random velocities [Angstrom/ps] from Gaussian distribution
    centered at zero with sigma according to Maxwell-Boltzmann
    distribution. Then rescales velocities to match specified desired
    temperature.
    """
    if self.temperature:
      self.etemp = self.temperature
      numpy.random.seed(self.random_seed)
      sigma_base = math.sqrt(2.0 * const.RGAS * self.temperature / const.NUMDIM)
      for atom in self.mol.atoms:
        sigma = sigma_base * atom.mass**(-0.5)
        atom.vels = numpy.random.normal(0.0, sigma, const.NUMDIM)

      self.mol.GetEnergy()
      self.mol.GetTemperature()

      vscale = math.sqrt(self.temperature / self.mol.temperature)
      for atom in self.mol.atoms:
        atom.vels *= vscale

  def _EquilibrateTemp(self):
    """Adjust velocities to equilibrate energy to set temperature.
    
    Computes exponential moving average of kinetic temperature and compares it
    to desired temperature. Velocities are then scaled depending on this ratio
    and the equilibration rate parameters.
    """
    tscale = self.timestep / max(self.timestep, self.eqrate)
    tweight = 10.0 * self.timestep
    self.etemp = (self.etemp + tweight * self.mol.temperature)/(1.0 + tweight)
    velscale = 1.0 + tscale * (math.sqrt(self.temperature / self.etemp) - 1.0)
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        self.mol.atoms[i].vels[j] *= velscale

  def _UpdateAccs(self):
    """Update accelerations of atoms [Angstrom/(ps^2)].
    
    Force is the negative gradient of the potential energy. Find
    accelerations by dividing the forces by the atomic masses.
    """
    for i in range(self.mol.n_atoms):
      mass = self.mol.atoms[i].mass
      for j in range(const.NUMDIM):
        self.mol.atoms[i].paccs[j] = self.mol.atoms[i].accs[j]
        self.mol.atoms[i].accs[j] = (
            -const.ACCCONV * self.mol.g_total[i][j] / mass)

  def _UpdateVels(self, dt):
    """Update velocities of atoms [Angstrom/ps].
    
    Acceleration is the derivative of velocity with respect to time.
    Find change by multiplying the acceleration by the timestep.
    
    Args:
      dt (float): time propogation increment [ps].
    """
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        self.mol.atoms[i].pvels[j] = self.mol.atoms[i].vels[j]
        self.mol.atoms[i].vels[j] += self.mol.atoms[i].accs[j] * dt

  def _UpdateCoords(self, dt):
    """Update coordinates of atoms [Angstrom].
    
    Velocity is the derivative of position with respect to time. Find
    displacement by multiplying the velocity by the timestep.
    
    Args:
      dt (float): time propogation increment [ps].
    """
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        self.mol.atoms[i].coords[j] += self.mol.atoms[i].vels[j] * dt
    self.mol.UpdateInternals()

  def _CheckPrint(self, timestep, print_all=False):
    """Check if printing of various md data is needed at current time.
    
    Args:
      timestep (float): Simulation time [ps] between previous check.
      print_all (bool): Print regardless of time status.
    """
    if print_all or self.etime >= self.energytime:
      self._PrintEnergy()
      self.etime = 1.0E-10
    if print_all or self.gtime >= self.geomtime:
      self._PrintGeom()
      self.gtime = 1.0E-10
    if print_all or (time.time() - self.stime) > self.statustime:
      self._PrintStatus()
      self.stime = time.time()
    self.etime += timestep
    self.gtime += timestep

  def _PrintEnergyHeader(self):
    """Print header of energy output columns to file."""
    e = self.efile
    e.write('#\n# INPUTFILE %s' % (self.infile))
    e.write('\n#\n# -- INPUT DATA --\n#')
    e.write('\n# MOLFILE %s' % (self.mol.infile))
    e.write('\n# ENERGYOUT %s' % (self.energyout))
    e.write('\n# GEOMOUT %s' % (self.geomout))
    e.write('\n# RANDOMSEED %i' % (self.random_seed))
    e.write('\n# TEMPERATURE %.6f K' % (self.temperature))
    e.write('\n# BOUNDARY %.6f A' % (self.mol.boundary))
    e.write('\n# BOUNDARYSPRING %.6f kcal/(mol*A^2)' % (self.mol.k_box))
    e.write('\n# BOUNDARYTYPE %s' % (self.mol.boundary_type))
    e.write('\n# STATUSTIME %.6f s' % (self.statustime))
    e.write('\n# ENERGYTIME %.6f ps' % (self.energytime))
    e.write('\n# GEOMTIME %.6f ps' % (self.geomtime))
    e.write('\n# TOTALTIME %.6f ps' % (self.tottime))
    e.write('\n# TIMESTEP %.6f ps' % (self.timestep))
    e.write('\n# EQTIME %.6f ps' % (self.eqtime))
    e.write('\n# EQRATE %.6f ps' % (self.eqrate))
    e.write('\n#\n# -- ENERGY DATA --\n#')
    e.write('\n# energy terms [kcal/mol]\n#  time      e_total      ')
    e.write('e_kin      e_pot  e_nonbond   e_bonded e_boundary      ')
    e.write('e_vdw     e_elst     e_bond    e_angle     e_tors      ')
    e.write('e_oop\n')


class MonteCarlo(Simulation):
  """Metropolis Monte-Carlo class derived from Simulation base class.
  
  Metropolis Monte-Carlo generates a set of configurations for a system by
  propogating an initial configuration forward in time through random
  displacements of coordinates which are accepted or rejected based on the
  relative potential energy Boltzmann factor.
  
  Args:
    infile_name (str): Path of simulation input file. May be relative or
    absolute, though absolute is safer.
      
  Attributes:
    totconf (int): Total number of configurations.
    conf (int): Current configuration number.
    dispmag (float): Magnitude of average random displacement [Angstrom].
    dispinc (float): Rate constant for 'dispmag' adjustment.
    n_accept (int): Number of accepted MMC trials.
    n_reject (int): Number of rejected MMC trials.
    dispconf (int): Number of configurations between adjusting 'dispmag' value.
    energyconf (int): Number of configurations between energy printing.
    geomconf (int): Number of configurations between geometry printing.
    rand_disp (float**): Random displacement vector of coordinates.
    cprintchar (int): Total characters for conf number printing.
  """
  def __init__(self, infile_name):
    self.simtype = 'mc'
    self.totconf = 1000
    self.conf = 0
    self.dispmag = 0.1
    self.dispinc = math.log(2.0)
    self.n_accept = 0
    self.n_reject = 0
    self.dispconf = 100
    self.energyconf = 100
    self.geomconf = 100
    Simulation.__init__(self, infile_name)
    self.rand_disp = numpy.zeros((self.mol.n_atoms, const.NUMDIM))
    self.cprintchar = 7

  def Run(self):
    """Run Metropolis Monte-Carlo according to simulation parameters.
    
    For every configuration, compute the potential energy and compare to the
    previous step. If the relative Boltzmann factor is above a random number,
    accept or else reject. When desired, alter the magnitude of random
    displacement to seek 50% acceptance. Print molecular geometry and/or energy
    data to output files as desired. Run until total configurations reached.
    """
    self._OpenOutputFiles()
    self._ZeroVels()
    numpy.random.seed(self.random_seed)
    self.mol.GetEnergy('standard')
    self._CheckPrint(0, print_all=True)
    previous_energy = self.mol.e_total
    while self.conf < self.totconf:
      self._GetRandDisp()
      self._DispCoords(self.rand_disp)
      self.mol.GetEnergy('standard')
      delta_e = self.mol.e_total - previous_energy
      bf = math.exp(min(1.0, -delta_e / (const.KB * self.temperature)))
      if bf >= numpy.random.random():
        self._CheckPrint(1)
        self.conf += 1
        self.n_accept += 1
        previous_energy = self.mol.e_total
      else:
        self._DispCoords(-self.rand_disp)
        self.n_reject += 1
      self._CheckDisp()
    self._CheckPrint(0)
    self._CloseOutputFiles()

  def _GetRandDisp(self):
    """Generate random displacment vector for coordinates.
    
    Random trial displacements for MMC are selected from a Gaussian distribution
    of mu = 0.0 and sigma = 'dispmag' attribute for all 3N atomic coordinates.
    """
    self.rand_disp.fill(0.0)
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        randval = numpy.random.normal(0.0, self.dispmag)
        self.rand_disp[i][j] = numpy.random.normal(0.0, self.dispmag)

  def _ZeroVels(self):
    """Set all 3N atomic velocity components to zero."""
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        self.mol.atoms[i].vels[j] = 0.0

  def _DispCoords(self, disp_vector):
    """Displace all 3N atomic coordinates by specified vector.
    
    Args:
      disp_vector (float**): Nx3 atomic displacement array [Angstrom].
    """
    for i in range(self.mol.n_atoms):
      for j in range(const.NUMDIM):
        self.mol.atoms[i].coords[j] += disp_vector[i][j]
    self.mol.UpdateInternals()

  def _ChangeDisp(self):
    """Change root-mean-square magnitude of displacement vector.
    
    The MMC random displacement vector has mu = 0.0, and sigma = 'dispmag'
    chosen to best approach 50% acceptance ratio. Increase 'dispmag' when
    'p_accept' > 0.5 and vice versa.
    """
    p_accept = float(self.n_accept) / float(self.n_reject + self.n_accept)
    self.n_accept, self.n_reject = 0, 0
    self.dispmag *= math.exp(2.0 * self.dispinc * (p_accept - 0.5))

  def _CheckPrint(self, n_conf, print_all=False):
    """Check if printing of various mc data is need at current time.
    
    Args:
      n_conf (int): Simulation configurations between previous check.
      print_all (bool): Print regardless of configuration status.
    """
    if print_all or self.econf >= self.energyconf:
      self._PrintEnergy()
      self.econf = 0
    if print_all or self.gconf >= self.geomconf:
      self._PrintGeom()
      self.gconf = 0
    if print_all or (time.time() - self.stime) > self.statustime:
      self._PrintStatus()
      self.stime = time.time()
    self.econf += n_conf
    self.gconf += n_conf

  def _CheckDisp(self):
    """Check if changing magnitude of random displacment vector needed."""
    if self.dconf >= self.dispconf:
      self._ChangeDisp()
      self.dconf = 0
    self.dconf += 1

  def _PrintEnergyHeader(self):
    """Print header of energy output columns to file."""
    e = self.efile
    e.write('#\n# INPUTFILE %s' % (self.infile))
    e.write('\n#\n# -- INPUT DATA --\n#')
    e.write('\n# MOLFILE %s' % (self.mol.infile))
    e.write('\n# ENERGYOUT %s' % (self.energyout))
    e.write('\n# GEOMOUT %s' % (self.geomout))
    e.write('\n# RANDOMSEED %i' % (self.random_seed))
    e.write('\n# TEMPERATURE %.6f K' % (self.temperature))
    e.write('\n# BOUNDARY %.6f A' % (self.mol.boundary))
    e.write('\n# BOUNDARYSPRING %.6f kcal/(mol*A^2)' % (self.mol.k_box))
    e.write('\n# BOUNDARYTYPE %s' % (self.mol.boundary_type))
    e.write('\n# STATUSTIME %.6f s' % (self.statustime))
    e.write('\n# ENERGYCONF %i' % (self.energyconf))
    e.write('\n# GEOMCONF %i' % (self.geomconf))
    e.write('\n# TOTALCONF %i' % (self.totconf))
    e.write('\n#\n# -- ENERGY DATA --\n#')
    e.write('\n# energy terms [kcal/mol] vs. configuration\n')
    e.write('#  conf        e_pot  e_nonbond   ')
    e.write('e_bonded e_boundary      e_vdw     e_elst     e_bond    ')
    e.write('e_angle     e_tors      e_oop\n')
