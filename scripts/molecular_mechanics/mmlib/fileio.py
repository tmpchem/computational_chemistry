"""Functions for printing, reading, and writing mmlib data.

Provides the interface to print molecular energy, gradient, geometry,
topology, parameter, simulation, and optimization data to screen or file.
"""

import math
import numpy
import os
import sys

from mmlib import constants as const
from mmlib import geomcalc
from mmlib import molecule
from mmlib import param
from mmlib import topology

# Formatting constants for print string creation functions.

_GEOM_PRINT_HEADER = ' Molecular Geometry and Non-bonded Parameters '
_GEOM_PRINT_BANNER_CHARS = 65
_GEOM_PRINT_PARAMS = ['type', 'x', 'y', 'z', 'q', 'ro/2', 'eps']
_GEOM_PRINT_SPACES = [6, 5, 9, 9, 7, 6, 4]
_GEOM_PRINT_ABSENT = '\nNo Atoms Detected'

_BOND_PRINT_HEADER = ' Bond Length Data '
_BOND_PRINT_BANNER_CHARS = 57
_BOND_PRINT_PARAMS = ['k_b', 'r_eq', 'r_ij', 'types', 'energy', 'atoms']
_BOND_PRINT_SPACES = [10, 5, 5, 3, 4, 1]
_BOND_PRINT_ABSENT = '\n No Bonds Detected'

_ANGLE_PRINT_HEADER = ' Bond Angle Data '
_ANGLE_PRINT_BANNER_CHARS = 58
_ANGLE_PRINT_PARAMS = ['k_a', 'a_eq', 'a_ijk', 'types', 'energy', 'atoms']
_ANGLE_PRINT_SPACES = [9, 3, 4, 4, 5, 2]
_ANGLE_PRINT_ABSENT = '\n No Bond Angles Detected'

_TORSION_PRINT_HEADER = ' Torsion Angle Data '
_TORSION_PRINT_BANNER_CHARS = 67
_TORSION_PRINT_PARAMS = ['vn/2', 'gamma', 't_ijkl n p', 'types', 'energy',
                         'atoms']
_TORSION_PRINT_SPACES = [9, 2, 3, 5, 6, 3]
_TORSION_PRINT_ABSENT = '\n No Torsion Angles Detected'

_OUTOFPLANE_PRINT_HEADER = ' Out-of-plane Angle Data '
_OUTOFPLANE_PRINT_BANNER_CHARS = 55
_OUTOFPLANE_PRINT_PARAMS = ['vn/2', 'o_ijkl', 'types', 'energy', 'atoms']
_OUTOFPLANE_PRINT_SPACES = [9, 2, 5, 6, 3]
_OUTOFPLANE_PRINT_ABSENT = '\n No Out-of-plane Angles Detected'

_ENERGY_PRINT_HEADER = ' Energy Values '
_ENERGY_PRINT_BANNER_CHARS = 33
_ENERGY_PRINT_PARAMS = ['component', '[kcal/mol]']
_ENERGY_PRINT_SPACES = [3, 9]

_AVERAGE_PRINT_HEADER = ' Energy Component Properties [kcal/mol] '
_AVERAGE_PRINT_BANNER_CHARS = 68
_AVERAGE_PRINT_PARAMS = ['component', 'avg', 'std', 'min', 'max']
_AVERAGE_PRINT_SPACES = [3, 11, 9, 9, 9]

# Formatting column headers and class attributes for energy components
_ENERGY_PRINT_LABELS = [
    'Total', 'Kinetic', 'Potential', 'Non-bonded', 'Bonded', 'Boundary',
    'van der Waals', 'Electrostatic', 'Bonds', 'Angles', 'Torsions',
    'Out-of-planes']
_ENERGY_PRINT_ATTRIBUTES = [
    'e_total', 'e_kinetic', 'e_potential', 'e_nonbonded', 'e_bonded',
    'e_bound', 'e_vdw', 'e_elst', 'e_bonds', 'e_angles', 'e_torsions',
    'e_outofplanes']

# Syntax use messages for programs within the molecular mechanics repository
_PROGRAM_MESSAGES = {
  'mm.py': 'xyzq or prm file for molecular mechanics\n',
  'md.py': 'simluation file for molecular dynamics\n',
  'mc.py': 'simulation file for Metropolis Monte Carlo\n',
  'opt.py': 'optimization file for energy minimization\n',
  'ana.py': 'plot file for data analysis\n'}

def GetFileStringArray(infile_name):
  """Creates a 2d array of strings from input file name.

  Each newline character creates a separate line element in the array.
  Each line is split by whitespace into an array of strings.

  Args:
    infile_name (str): Path to desired text file for reading. May be relative or
        absolute. Exits on error if doesn't exist.
  
  Returns:
    (str**): Contents of text file as an array of arrays of strings
  """
  if not os.path.exists(infile_name):
    raise ValueError(
        'attempted to read from file which does not exist: %s' % (infile_name))

  with open(infile_name, 'r') as infile:
    lines = infile.readlines()

  return [line.split() for line in lines]


def _IsType(type, value):
  """Determines if a value can be cast as given type.

  Args:
    type (type): Type to attempt casting to.
    value (any): Value to attempt being casted.

  Returns:
    is_type (bool): Whether 'value' can be cast to 'type'.
  """
  try:
    type(value)
  except ValueError:
    return False
  return True


def _AreAllType(type, values):
  """Determines if an array of values can all be cast as given type.

  Args:
    type (type): Type to attempt casting to.
    values (any*): Array of values to attempt being casted.

  Returns:
    are_type (bool): Whether all 'values' can be cast to 'type'.
  """
  for value in values:
    if not _IsType(type, value):
      return False
  return True


def GetAtomFromXyzq(row):
  """Parses and validates a row of input from xyzq file into an Atom object.

  Args:
    row (str*): Array of strings from row of xyzq file.

  Returns:
    atom (mmlib.molecule.Atom): Atom object.

  Raises:
    IndexError: If insufficient number of fields.
    ValueError: If incorrect data type of any field.
  """
  if len(row) < 2+const.NUMDIM:
    raise IndexError('Insufficient columns to parse xyzq file row into Atom '
                     'object: %s' % ' '.join(row))

  type_ = row[0]
  coords = row[1:1+const.NUMDIM]
  charge = row[1+const.NUMDIM]

  if not type_ or not type_[0].isalpha():
    raise ValueError('Atom type must begin with letter: %s' % type_)
  if not _AreAllType(float, coords):
    raise ValueError(
        'Atomic coordinates must be numeric values: %s' % ' '.join(coords))
  if not _IsType(float, charge):
    raise ValueError('Atomic charge must be numeric value: %s' % charge)

  return molecule.Atom(type_, numpy.fromiter(coords, float), float(charge))


def GetAtomsFromXyzq(rows):
  """Parses molecular geometry data from xyzq file.
  
  First line contains (int) number of atoms. Second line is ignored comment.
  Each line afterward (3 to [n+2]) contains atom type, (float) 3 xyz Cartesian
  coordinates [Angstrom], and (float) charge [e].

  Args:
    rows (str**): 2d array of string from xyzq input file.

  Returns:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Raises:
    EOFError: If input is empty or of insufficient size.
    IndexError: If first line is empty.
  """
  if not rows:
    raise EOFError('XYZQ file is empty.')
  if not rows[0]:
    raise IndexError('First line of XYZQ file must be number of atoms.')

  n_atoms = int(rows[0][0])

  if n_atoms <= 0:
    raise ValueError('Number of atoms must be positive: %i' % n_atoms)
  if len(rows) < n_atoms+2:
    raise EOFError('XYZQ file does not contain enough lines for stated number '
                   'of atoms: %i' % n_atoms)

  atoms = []
  for row in rows[2:n_atoms+2]:
    atoms.append(GetAtomFromXyzq(row))
  return atoms


def GetAtomFromPrm(row):
  """Parses atom row into an Atom object.

  Args:
    row (str*): Array of strings from line of prm file.

  Returns:
    atom (mmlib.molecule.Atom): Atom object with attributes from row.

  Raises:
    IndexError: If insufficient number of fields.
    ValueError: If incorrect data type of any field.
  """
  if len(row) < 6+const.NUMDIM:
    raise IndexError('Insufficient columns to parse prm file row into Atom '
                     'object: %s' % ' '.join(row))

  type_ = row[2]
  coords = row[3:3+const.NUMDIM]
  charge, ro, eps = row[3+const.NUMDIM:6+const.NUMDIM]

  if not type_ or not type_[0].isalpha():
    raise ValueError('Atom type must begin with letter: %s' % type_)
  if not _AreAllType(float, coords):
    raise ValueError(
        'Atomic coordinates must be numeric values: %s' % ' '.join(coords))
  if not _IsType(float, charge):
    raise ValueError('Atomic charge must be numeric value: %s' % charge)
  if not _IsType(float, ro) or not float(ro) > 0.0:
    raise ValueError(
        'Atomic vdw radius must be positive numeric value: %s' % ro)
  if not _IsType(float, eps) or not float(eps) >= 0.0:
    raise ValueError(
        'Atomic epsilon must be non-negative numeric value: %s' % eps)

  return molecule.Atom(type_, numpy.fromiter(coords, float), float(charge), 
                       float(ro), float(eps))


def GetBondFromPrm(row):
  """Parses bond row into a Bond object.

  Args:
    row (str*): Array of strings from line of prm file.

  Returns:
    bond (mmlib.molecule.Bond): Bond object with attributes from row.

  Raises:
    IndexError: If insufficient number of fields.
    ValueError: If incorrect data type of any field.
  """
  if len(row) < 5:
    raise IndexError('Insufficient columns to parse prm file row into Bond '
                     'object: %s' % ' '.join(row))

  at1, at2, k_b, r_eq = row[1:5]

  if not _IsType(int, at1) or not int(at1) > 0:
    raise ValueError('Bond atom1 index must be positive integer: %s' % at1)
  if not _IsType(int, at2) or not int(at2) > 0:
    raise ValueError('Bond atom2 index must be positive integer: %s' % at2)
  if not _IsType(float, k_b) or not float(k_b) > 0.0:
    raise ValueError(
        'Bond spring constant must be positive numeric value: %s' % k_b)
  if not _IsType(float, r_eq) or not float(r_eq) > 0.0:
    raise ValueError(
        'Equilibrium bond length must be positive numeric value: %s' % r_eq)

  return molecule.Bond(int(at1)-1, int(at2)-1, float(k_b), float(r_eq))


def GetAngleFromPrm(row):
  """Parses angle row into an Angle object.

  Args:
    row (str*): Array of strings from line of prm file.

  Returns:
    angle (mmlib.molecule.Angle): Angle object with attributes from row.

  Raises:
    IndexError: If insufficient number of fields.
    ValueError: If incorrect data type of any field.
  """
  if len(row) < 6:
    raise IndexError('Insufficient columns to parse prm file row into Angle '
                     'object: %s' % ' '.join(row))

  at1, at2, at3, k_a, a_eq = row[1:6]

  if not _IsType(int, at1) or not int(at1) > 0:
    raise ValueError('Angle atom1 index must be positive integer: %s' % at1)
  if not _IsType(int, at2) or not int(at2) > 0:
    raise ValueError('Angle atom2 index must be positive integer: %s' % at2)
  if not _IsType(int, at3) or not int(at3) > 0:
    raise ValueError('Angle atom3 index must be positive integer: %s' % at3)
  if not _IsType(float, k_a) or not float(k_a) > 0.0:
    raise ValueError(
        'Angle spring constant must be positive numeric value: %s' % k_a)
  if not _IsType(float, a_eq) or not 0.0 <= float(a_eq) <= 180.0:
    raise ValueError(
        'Equilibrium bond angle must be between 0.0 and 180.0: %s' % a_eq)

  return molecule.Angle(
      int(at1)-1, int(at2)-1, int(at3)-1, float(k_a), float(a_eq))


def GetTorsionFromPrm(row):
  """Parses torsion row into an Torsion object.

  Args:
    row (str*): Array of strings from line of prm file.

  Returns:
    torsion (mmlib.molecule.Torsion): Torsion object.

  Raises:
    IndexError: If insufficient number of fields.
    ValueError: If incorrect data type of any field.
  """
  if len(row) < 9:
    raise IndexError('Insufficient columns to parse prm file row into Torsion '
                     'object: %s' % ' '.join(row))

  at1, at2, at3, at4, v_n, gamma, nfold, paths = row[1:9]

  if not _IsType(int, at1) or not int(at1) > 0:
    raise ValueError('Torsion atom1 index must be positive integer: %s' % at1)
  if not _IsType(int, at2) or not int(at2) > 0:
    raise ValueError('Torsion atom2 index must be positive integer: %s' % at2)
  if not _IsType(int, at3) or not int(at3) > 0:
    raise ValueError('Torsion atom3 index must be positive integer: %s' % at3)
  if not _IsType(int, at4) or not int(at4) > 0:
    raise ValueError('Torsion atom4 index must be positive integer: %s' % at4)
  if not _IsType(float, v_n) or not float(v_n) > 0.0:
    raise ValueError(
        'Torsion half-barrier must be positive numeric value: %s' % v_n)
  if not _IsType(float, gamma) or not -180.0 <= float(gamma) <= 180.0:
    raise ValueError(
        'Torsion offset angle must be between -180.0 and 180.0: %s' % gamma) 
  if not nfold.isdigit() or not int(nfold) > 0:
    raise ValueError(
        'Torsion barrier frequency must be positive integer: %s' % nfold)
  if not paths.isdigit() or not int(paths) > 0:
    raise ValueError(
        'Torsion path number must be positive integer: %s' % paths)

  return molecule.Torsion(
      int(at1)-1, int(at2)-1, int(at3)-1, int(at4)-1, float(v_n), float(gamma),
      int(nfold), int(paths))


def GetOutofplaneFromPrm(row):
  """Parses outofplane row into an Outofplane object.

  Args:
    row (str*): Array of strings from line of prm file.

  Returns:
    outofplane (mmlib.molecule.Outofplane): Outofplane object.
  """
  if len(row) < 6:
    raise IndexError('Insufficient columns to parse prm file row in Outofplane '
                     'object: %s' % ' '.join(row))

  at1, at2, at3, at4, v_n = row[1:6]

  if not _IsType(int, at1) or not int(at1) > 0:
    raise ValueError(
        'Outofplane atom1 index must be positive integer: %s' % at1)
  if not _IsType(int, at2) or not int(at2) > 0:
    raise ValueError(
        'Outofplane atom2 index must be positive integer: %s' % at2)
  if not _IsType(int, at3) or not int(at3) > 0:
    raise ValueError(
        'Outofplane atom3 index must be positive integer: %s' % at3)
  if not _IsType(int, at4) or not int(at4) > 0:
    raise ValueError(
        'Outofplane atom4 index must be positive integer: %s' % at4)
  if not _IsType(float, v_n) or not float(v_n) > 0.0:
    raise ValueError('Outofplane half-barrier must be numeric value: %s' % v_n)

  return molecule.Outofplane(
      int(at1)-1, int(at2)-1, int(at3)-1, int(at4)-1, float(v_n))


def GetAtomsFromPrm(rows):
  """Parses atom rows into an array of Atom objects.

  Args:
    rows (str**): 2d array of strings from lines of prm file.

  Returns:
    atoms (mmlib.molecule.Atom*): Array of Atom object with parameters from
        rows.
  """
  atoms = []
  for row in rows:
    if row and row[0].upper() == 'ATOM':
      atoms.append(GetAtomFromPrm(row))
  return atoms


def GetBondsFromPrm(rows):
  """Parses bond rows into an array of Bond objects.
  
  Args:
    rows (str**): 2d array of strings from lines of prm file.

  Returns:
    bonds (mmlib.molecule.Bond*): Array of Bond objects with parameters from
        rows.
  """
  bonds = []
  for row in rows:
    if row and row[0].upper() == 'BOND':
      bonds.append(GetBondFromPrm(row))
  return bonds


def GetAnglesFromPrm(rows):
  """Parses angle rows into an array of Angle objects.
  
  Args:
    rows (str**): 2d array of strings from lines of prm file.

  Returns:
    angles (mmlib.molecule.Angle): Array of Angle objects with parameters from
        rows.
  """
  angles = []
  for row in rows:
    if row and row[0].upper() == 'ANGLE':
      angles.append(GetAngleFromPrm(row))
  return angles


def GetTorsionsFromPrm(rows):
  """Parses torsion rows into an array of Torsion objects.
  
  Args:
    rows (str**): 2d array of strings from lines of prm file.

  Returns:
    torsions (mmlib.molecule.Torsion*): Array of Torsion objects with parameters
        from row.
  """
  torsions = []
  for row in rows:
    if row and row[0].upper() == 'TORSION':
      torsions.append(GetTorsionFromPrm(row))
  return torsions


def GetOutofplanesFromPrm(rows):
  """Parses outofplane row into Outofplane object.
  
  Args:
    rows (str**): 2d array of strings from lines of prm file.

  Returns:
    outofplanes (mmlib.molecule.Outofplane*): Array of Outofplane objects with
        parameters from row.
  """
  outofplanes = []
  for row in rows:
    if row and row[0].upper() == 'OUTOFPLANE':
      outofplanes.append(GetOutofplaneFromPrm(row))
  return outofplanes


def GetSimData(sim):
  """Parses contents of sim file into molecular simulation data.
  
  Many molecular simulation parameters (dynamics, monte carlo, etc.) can be
  determined by default, or overriden in a simulation file. All listed values
  below can be set through the given keyword arguments.
  
  Mandatory values include (str) molecule [file path], (float) temperature
  [Kelvin], and (float / int) total time/confs [ps / none] (md / mc).
  
  Args:
    sim (mmlib.simulate.Simulation): Simulation object to append data.
  """
  infile_array = GetFileStringArray(sim.infile)
  cwd = os.getcwd()
  os.chdir(sim.indir)
  for q in range(len(infile_array)):
    if len(infile_array[q]) < 2:
      continue
    kwarg = infile_array[q][0].lower()
    kwargval = infile_array[q][1]
    kwargarr = infile_array[q][1:]
    if kwarg == 'molecule':
      sim.mol = molecule.Molecule(os.path.realpath(kwargval))
    elif kwarg == 'temperature':
      sim.temperature = float(kwargval)
    elif kwarg == 'pressure':
      sim.pressure = float(kwargval)
    elif kwarg == 'boundaryspring':
      sim.mol.k_box = float(kwargval)
    elif kwarg == 'boundary':
      sim.mol.boundary = float(kwargval)
      sim.mol.GetVolume()
    elif kwarg == 'boundarytype':
      sim.mol.boundary_type = kwargval.lower()
      sim.mol.GetVolume()
    elif kwarg == 'origin':
      sim.mol.origin = list(map(float, kwargarr[:const.NUMDIM]))
    elif kwarg == 'totaltime':
      sim.tottime = float(kwargval)
    elif kwarg == 'totalconf':
      sim.totconf = int(kwargval)
    elif kwarg == 'timestep':
      sim.timestep = float(kwargval)
    elif kwarg == 'geomtime':
      sim.geomtime = float(kwargval)
    elif kwarg == 'geomconf':
      sim.geomconf = int(kwargval)
    elif kwarg == 'geomout':
      sim.geomout = os.path.realpath(kwargval)
    elif kwarg == 'energytime':
      sim.energytime = float(kwargval)
    elif kwarg == 'energyconf':
      sim.energyconf = int(kwargval)
    elif kwarg == 'energyout':
      sim.energyout = os.path.realpath(kwargval)
    elif kwarg == 'statustime':
      sim.statustime = float(kwargval)
    elif kwarg == 'eqtime':
      sim.eqtime = float(kwargval)
    elif kwarg == 'eqrate':
      sim.eqrate = float(kwargval)
    elif kwarg == 'randomseed':
      sim.random_seed = int(kwargval) % 2**32
  os.chdir(cwd)


def GetOptData(opt):
  """Parses contents of sim file into molecular simulation data.
 
  Many molecular energy minimization parameters can be determined by default, or
  overridden in an optimization file. All listed values below can be set through
  the give keyword arguments.

  The only mandatory value is (str) molecule [file path]. Setting (str) geomout
  and (str) energyout is also strongly recommended.

  Args:
    opt (mmlib.optimize.Optimization): Optimization object to append data.
  """
  infile_array = GetFileStringArray(opt.infile)
  cwd = os.getcwd()
  os.chdir(opt.indir)
  for q in range(len(infile_array)):
    if len(infile_array[q]) < 2:
      continue
    kwarg = infile_array[q][0].lower()
    kwargval = infile_array[q][1]
    kwargarr = infile_array[q][1:]
    if kwarg == 'molecule':
      opt.mol = molecule.Molecule(os.path.realpath(kwargval))
    elif kwarg == 'opttype':
      opt.opt_type = kwargval.lower()
    elif kwarg == 'optcriteria':
      opt.opt_str = kwargval.lower()
    elif kwarg == 'e_converge':
      opt.conv_delta_e = float(kwargval)
    elif kwarg == 'grms_converge':
      opt.conv_grad_rms = float(kwargval)
    elif kwarg == 'gmax_converge':
      opt.conv_grad_max = float(kwargval)
    elif kwarg == 'drms_converge':
      opt.conv_disp_rms = float(kwargval)
    elif kwarg == 'dmax_converge':
      opt.conv_disp_max = float(kwargval)
    elif kwarg == 'nmaxiter':
      opt.n_maxiter = float(kwargval)
    elif kwarg == 'geomout':
      opt.geomout = os.path.realpath(kwargval)
    elif kwarg == 'energyout':
      opt.energyout = os.path.realpath(kwargval)
  os.chdir(cwd)


def GetAnalysisData(ana):
  """Parses contents of plt file into ensemble analysis data.
 
  Many simulation data analysis keywords can be determined by default, or
  overridden in a plot file. All listed values below can be set through the
  given keyword arguments.

  Mandatory values include (str) 'input' [file path] and (str) 'simtype'.
  Setting 'plotout' is also strongly recommended.
  
  Args:
    ana (mmlib.analyze.Analysis): Analysis object to append data.
  """
  infile_array = GetFileStringArray(ana.infile)
  cwd = os.getcwd()
  os.chdir(ana.indir)
  for q in range(len(infile_array)):
    if len(infile_array[q]) < 2:
      continue
    kwarg = infile_array[q][0].lower()
    kwargval = infile_array[q][1]
    kwargarr = infile_array[q][1:]
    if kwarg == 'input':
      ana.simfile = os.path.realpath(kwargval)
      ana.simdir = os.path.dirname(ana.simfile)
    elif kwarg == 'simtype':
      ana.simtype = kwargval.lower()
    elif kwarg == 'plotout':
      ana.plotout = os.path.realpath(kwargval)
    elif kwarg == 'percentstart':
      ana.percent_start = float(kwargval)
    elif kwarg == 'percentstop':
      ana.percent_stop = float(kwargval)
  os.chdir(cwd)


def GetProperties(prop_file):
  """Parses molecular property sequences from simulation data file.
  
  Input file contains a commented (#) header with column identifier keys and
  lines of snapshot data at various configurations of a molecular simulation,
  identified either by time [ps] or configuration number.
  
  First find the line the key labels are on, then find how many and which lines
  contain data. Then find which column corresponds to each key, and populate the
  data arrays into the dictionary entry of each key.
  
  Args:
    prop_file (str): Path to input property file.
  
  Returns:
    prop (float**): Dictionary of property keys with array values from each
        configuration of molecule during trajectory.
  """
  prop_array = GetFileStringArray(prop_file)
  n_lines = len(prop_array)
  prop_keys = const.PROPERTYKEYS
  key1 = prop_keys[2]

  key_line = 0
  for i in range(n_lines):
    if key1 in prop_array[i]:
      key_line = i
      break
  key1_col = prop_array[key_line].index(key1)
  n_keys = len(prop_array[key_line]) - 1

  n_confs = 0
  excluded_lines = []
  for i in range(len(prop_array)):
    if '#' in prop_array[i][0] or not len(prop_array[i]) == n_keys:
      excluded_lines.append(i)
    else:   
      n_confs += 1

  prop = {}
  for j in range(n_keys):
    key = prop_array[key_line][j+1]
    prop[key] = numpy.zeros(n_confs)
    confnum = 0
    for i in range(n_lines):
      if not i in excluded_lines:
        prop[key][confnum] = float(prop_array[i][j])
        confnum += 1

  return prop


def GetTrajectory(traj_file):
  """Parses molecular xyz coordinate sequences from xyz file.
  
  Input file contains 'n_confs' molecular xyz-coordinate snapshots with
  'n_atoms' atoms each. Start by reading in 'n_atoms', and continue until file
  ends.
  
  Args:
    traj_file (str): Path to input trajectory file.
  
  Returns:
    traj (float***): Array of xyz-coordinates at each configuration of molecule
        during trajectory.
  """
  traj_array = GetFileStringArray(traj_file)
  n_lines = len(traj_array)
  n_atoms = int(traj_array[0][0])
  n_confs = int(math.floor(n_lines / (n_atoms+2)))
  traj = numpy.zeros((n_confs, n_atoms, const.NUMDIM))
  for p in range(n_confs):
    geom_start = p * (n_atoms+2)
    for i in range(n_atoms):
      atom_start = geom_start + i + 2
      for j in range(const.NUMDIM):
        traj[p][i][j] = float(traj_array[atom_start][j+1])
  return traj


def GetPrintCoordsXyzString(atoms, comment, total_chars=12, decimal_chars=6):
  """Creates string of xyz-formatted coordinates for a set of atoms.
  
  Print to screen all (float) 3N atomic cartesian coordinates [Angstrom] from
  mol in xyz file format with (str) 'comment' for the comment line.
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    comment (str): Comment string for xyz file comment line.
    total_chars (int): Number of characters in float print string.
    decimal_chars (int): Number of post-decimal digits in float print string.

  Returns:
    string (str): String of xyz-formatted coordinate snapshot for printing.
  """
  string = ['%i\n%s\n' % (len(atoms), comment)]
  for atom in atoms:
    string.append('%-2s' % atom.element)
    for j in range(const.NUMDIM):
      string.append(' %*.*f' % (total_chars, decimal_chars, atom.coords[j]))
    string.append('\n')
  return ''.join(string)


def GetPrintGradientString(atoms, gradient, comment, total_chars=12,
                           decimal_chars=6):
  """Creates specified atomic gradient type string for a set of atoms.
  
  Generate print string for all (float) 3N atomic cartesian gradient components
  [kcal/(mol*Angstrom)] from mol in xyz file format. Gradient is partial
  derivative of 'grad_type' energy [kcal/mol] with respect to each (float)
  cartesian coordinate [Angstrom].
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.
    gradient (numpy.float**): Nx3 Molecular energy gradient [kcal/(mol*A)].
    comment (str): Comment on gradient type / source / etc. to print to screen.
    total_chars (int): Number of characters in float print string.
    decimal_chars (int): Number of post-decimal digits in float print string.

  Returns:
    string (str): String of xyz-formatted gradient snapshot for printing.
  """
  string = ['\n %s\n' % (comment)]
  for i in range(len(gradient)):
    string.append('%-2s' % atoms[i].type_)
    for j in range(const.NUMDIM):
      string.append(' %*.*f' % (total_chars, decimal_chars, grad[i][j]))
    string.append('\n')
  return ''.join(string)


def _GetBannerString(title, length, lines1, lines2):
  """Creates a string in the center of a banner with dashes on each side.

  Args:
    title (str): Banner header title.
    length (int): Total number of characters in banner.
    lines1 (int): Number of leading newlines.
    lines2 (int): Number of trailing newlines.

  Returns:
    string (str): Resulting print header string.
  """
  dashes1 = (length - len(title))//2 - 1
  dashes2 = length - len(title) - dashes1 - 2
  return lines1*'\n' + dashes1*'-' + title + dashes2*'-' + lines2*'\n'


def _GetPaddedString(fields, spacings):
  """Creates a string of an array of strings padded by spaces.

  Args:
    fields (str*): Array of string fields to be printed.
    spacings (int*): Array of number of spaces to be printed prior to each
        'strings' element.

  Returns:
    string (str): Resulting print field string.
  """
  string = []
  for spacing, field in zip(spacings, fields):
    string.append(spacing*' ' + field)
  return ''.join(string) + '\n'


def _GetHeaderString(header, n_banner, params, spacings):
  """Creates a banner header string for a section of parameter outputs.
  
  Args:
    header (str): Banner header title.
    n_banner (int): Total number of characters in banner.
    params (str*): Array of strings to be printed.
    spacings (int*): Array of number of spaces to be printed prior to each
        'params' element.

  Returns:
    string (str): Resulting print banner string.
  """
  return ''.join([
    _GetBannerString(header, n_banner, 1, 1),
    _GetPaddedString(params, spacings),
    _GetBannerString('', n_banner, 0, 0)])


def GetPrintGeomString(atoms):
  """Creates a string of molecular geometry and non-bonded parameters.
  
  Create a header banner for the section, and append the following fields for
  each Atom object:
    * (int) atomic index
    * (str) atom type
    * (float) 3 xyz cartesian coordinates [Angstrom]
    * (float) atomic partial charge [e]
    * (float) atomic van der waals radius [Angstrom]
    * (float) atomic van der waals epsilon [kcal/mol]
  
  Args:
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    string (str): String of formatted atom coordinates and parameters to print.
  """
  if not atoms:
    return _ATOM_PRINT_ABSENT

  string = [_GetHeaderString(_GEOM_PRINT_HEADER, _GEOM_PRINT_BANNER_CHARS,
                             _GEOM_PRINT_PARAMS, _GEOM_PRINT_SPACES), '\n']

  for i, atom in enumerate(atoms):
    string.append('%4i | %-2s' % (i+1, atom.type_))
    for j in range(const.NUMDIM):
      string.append('%10.4f' % atom.coords[j])
    string.append(' %7.4f %7.4f %7.4f\n' % (atom.charge, atom.ro, atom.eps))
  return ''.join(string)


def GetPrintBondsString(bonds, atoms):
  """Creates a string of molecular bond parameters for output.

  Create a header banner for the section, and append the following fields for
  each Bond object:
    * (int) bond index
    * (float) bond spring constant [kcal/(mol*A^2)]
    * (float) equilibrium bond length [Angstrom]
    * (float) bond length [Angstrom]
    * (str) 2 bond atom types
    * (float) bond energy [kcal/mol]
    * (int) 2 bond atom atomic indices
  
  Args:
    bonds (mmlib.molecule.Bond*): Array of molecule's Bond objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    string (str): String of formatted bond parameters to print.
  """
  if not bonds:
    return _BOND_PRINT_ABSENT

  string = [_GetHeaderString(_BOND_PRINT_HEADER, _BOND_PRINT_BANNER_CHARS,
                             _BOND_PRINT_PARAMS, _BOND_PRINT_SPACES)]

  for p, bond in enumerate(bonds):
    t1, t2 = atoms[bond.at1].type_, atoms[bond.at2].type_
    string.append('%4i | %7.2f %8.4f %8.4f (%2s-%2s) %8.4f (%i-%i)' % (
        p+1, bond.k_b, bond.r_eq, bond.r_ij, t1, t2, bond.energy, bond.at1+1,
        bond.at2+1))
  return '\n'.join(string)


def GetPrintAnglesString(angles, atoms):
  """Creates a string of molecular angle parameters for output.

  Create a header banner for the section, and append the following fields for
  each Angle object:
    * (int) angle index
    * (float) angle spring constant [kcal/(mol*A^2)]
    * (float) equilibrium bond angle [degrees]
    * (float) bond angle [degrees]
    * (str) 3 angle atom types
    * (float) angle energy [kcal/mol]
    * (int) 3 angle atom atomic indices
  
  Args:
    angles (mmlib.molecule.Angle*): Array of molecule's Angle objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    string (str): String of formatted angle parameters to print.
  """
  if not angles:
    return _ANGLE_PRINT_ABSENT

  string = [_GetHeaderString(_ANGLE_PRINT_HEADER, _ANGLE_PRINT_BANNER_CHARS,
                             _ANGLE_PRINT_PARAMS, _ANGLE_PRINT_SPACES)]

  for p, angle in enumerate(angles):
    t1, t2, t3 = (
        atoms[angle.at1].type_, atoms[angle.at2].type_, atoms[angle.at3].type_)
    string.append('%4i | %6.2f %7.3f %7.3f (%2s-%2s-%2s) %7.4f (%i-%i-%i)' % (
        p+1, angle.k_a, angle.a_eq, angle.a_ijk, t1, t2, t3, angle.energy,
        angle.at1+1, angle.at2+1, angle.at3+1))
  return '\n'.join(string)


def GetPrintTorsionsString(torsions, atoms):
  """Creates a string of molecular torsion parameters for output.

  Create a header banner for the section, and append the following fields for
  each Torsion object:
    * (int) torsion index
    * (float) torsion half-barrier [kcal/mol]
    * (float) torsion offset angle [degrees]
    * (float) torsion angle [degrees]
    * (int) n-fold periodicity of barrier
    * (int) number of unique paths around central atom pair
    * (str) 4 torsion atom types
    * (float) torsion energy [kcal/mol]
    * (int) 4 torsion atom atomic indices

  Args:
    torsions (mmlib.molecule.Torsion*): Array of molecule's Torsion objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    string (str): String of formatted torsion parameters to print.
  """
  if not torsions:
    return _TORSION_PRINT_ABSENT

  string = [_GetHeaderString(_TORSION_PRINT_HEADER, _TORSION_PRINT_BANNER_CHARS,
                             _TORSION_PRINT_PARAMS, _TORSION_PRINT_SPACES)]

  for p, torsion in enumerate(torsions):
    t1, t2, t3, t4 = (
        atoms[torsion.at1].type_, atoms[torsion.at2].type_,
        atoms[torsion.at3].type_, atoms[torsion.at4].type_)
    string.append(
        '%4i | %6.2f %6.1f %8.3f %i %i (%2s-%2s-%2s-%2s) %7.4f (%i-%i-%i-%i)'
        % (p+1, torsion.v_n, torsion.gam, torsion.t_ijkl, torsion.n,
        torsion.paths, t1, t2, t3, t4, torsion.energy, torsion.at1+1,
        torsion.at2+1, torsion.at3+1, torsion.at4+1))
  return '\n'.join(string)


def GetPrintOutofplanesString(outofplanes, atoms):
  """Creates a string of molecular outofplane parameters for output.

  Create a header banner for the section, and append the following fields for
  each Outofplane object:
    * (int) outofplane index
    * (float) outofplane half-barrier [kcal/mol]
    * (float) outofplane angle [degrees]
    * (str) 4 outofplane atom types
    * (float) outofplane energy [kcal/mol]
    * (int) 4 outofplane atom atomic indices

  Args:
    outofplanes (mmlib.molecule.Outofplane*):Array of Outofplane objects.
    atoms (mmlib.molecule.Atom*): Array of molecule's Atom objects.

  Returns:
    string (str): String of formatted outofplane parameters to print.
  """
  if not outofplanes:
    return _OUTOFPLANE_PRINT_ABSENT

  string = [_GetHeaderString(
      _OUTOFPLANE_PRINT_HEADER, _OUTOFPLANE_PRINT_BANNER_CHARS,
      _OUTOFPLANE_PRINT_PARAMS, _OUTOFPLANE_PRINT_SPACES)]
  
  for p, outofplane in enumerate(outofplanes):
    t1, t2, t3, t4 = (
        atoms[outofplane.at1].type_, atoms[outofplane.at2].type_,
        atoms[outofplane.at3].type_, atoms[outofplane.at4].type_)
    string.append('%4i | %6.2f %7.3f (%2s-%2s-%2s-%2s) %7.4f (%i-%i-%i-%i)' % (
        p+1, outofplane.v_n, outofplane.o_ijkl, t1, t2, t3, t4,
        outofplane.energy, outofplane.at1+1, outofplane.at2+1, outofplane.at3+1,
        outofplane.at4+1))
  return '\n'.join(string)


def GetPrintEnergyString(mol):
  """Creates a string of molecular energy components for output.

  Create a header banner for the section, and append each energy component label
  and value [kcal/mol].

  Args:
    mol (mmlib.molecule.Molecule): Molecule object with energy component data.

  Returns:
    string (str): String of formatted energy component values to print.
  """
  string = [_GetHeaderString(_ENERGY_PRINT_HEADER, _ENERGY_PRINT_BANNER_CHARS,
                             _ENERGY_PRINT_PARAMS, _ENERGY_PRINT_SPACES)]
  
  for label, attribute in zip(_ENERGY_PRINT_LABELS, _ENERGY_PRINT_ATTRIBUTES):
    string.append('   %-13s | %10.4f' % (label, getattr(mol, attribute)))
  return '\n'.join(string)


def PrintAverages(ana):
  """Print list of expectation values in a table to screen.
  
  For each energy term in mmlib.analyze.Analysis.pdict print (float) average,
  standard deviation, minimum, and maximum value [kcal/mol] to screen.

  Args:
    ana (mmlib.analyze.Analyze): Analyze object with energy component
        expectation values.
  """
  string = [_GetHeaderString(_AVERAGE_PRINT_HEADER, _AVERAGE_PRINT_BANNER_CHARS,
                             _AVERAGE_PRINT_PARAMS, _AVERAGE_PRINT_SPACES)]
  
  pdict = const.PROPERTYDICTIONARY
  keys = sorted(list(pdict.keys()), key = lambda x: pdict[x][3])
  keys = [key for key in keys if key in ana.prop]
  labels = [pdict[key][0] for key in keys]
  for key, label in zip(keys, labels):
    string.append('   %-13s | %11.4e %11.4e %11.4e %11.4e' % (
        label, ana.eavg[key], ana.estd[key], ana.emin[key], ana.emax[key]))
  return '\n'.join(string)

def ValidateInput(program_path, argv):
  """Checks for proper input argument syntax and return parsed result.
  
  Check that the command line input contains two strings or throw an error and
  print usage guidance. If correct, return the name of the input file given as
  the second command line input string.

  Args:
    program_path (str): Path to main file calling this function.
    argv (str*): Array of strings from command line input.

  Returns:
    infile_name (str): Name of input file given from command line.

  Raises:
    ValueError: If program_name not recognized.
    FileNotFoundError: If input file is not in file system path.
  """
  if (len(argv) < 2):
    program_name = program_path.split('/')[-1]
    if program_name in _PROGRAM_MESSAGES:
      print('\nUsage: python %s INPUT_FILE\n' % program_name)
      print('INPUT_FILE: %s' % _PROGRAM_MESSAGES[program_name])
      sys.exit()
    else:
      raise ValueError('Program name not recognized: %s' % program_name)
  
  input_file = argv[1]
  if os.path.isfile(input_file):
    return input_file
  else:
    raise FileNotFoundError('Specified input is not a file: %s' % input_file)
