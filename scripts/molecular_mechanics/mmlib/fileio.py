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

_GEOM_PRINT_HEADER = ' Molecular Geometry and Non-bonded Parameters '
_GEOM_PRINT_BANNER_CHARS = 65
_GEOM_PRINT_PARAMS = ['type', 'x', 'y', 'z', 'q', 'ro/2', 'eps']
_GEOM_PRINT_SPACES = [6, 5, 9, 9, 7, 6, 4]

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

_ENERGY_PRINT_LABELS = [
    'Total', 'Kinetic', 'Potential', 'Non-bonded', 'Bonded', 'Boundary',
    'van der Waals', 'Electrostatic', 'Bonds', 'Angles', 'Torsions',
    'Out-of-planes']
_ENERGY_PRINT_ATTRIBUTES = [
    'e_total', 'e_kinetic', 'e_potential', 'e_nonbonded', 'e_bonded',
    'e_bound', 'e_vdw', 'e_elst', 'e_bonds', 'e_angles', 'e_torsions',
    'e_outofplanes']

_AVERAGE_PRINT_HEADER = ' Energy Component Properties [kcal/mol] '
_AVERAGE_PRINT_BANNER_CHARS = 68
_AVERAGE_PRINT_PARAMS = ['component', 'avg', 'std', 'min', 'max']
_AVERAGE_PRINT_SPACES = [3, 11, 9, 9, 9]

_PROGRAM_MESSAGES = {
  'mm.py': 'xyzq or prm file for molecular mechanics\n',
  'md.py': 'simluation file for molecular dynamics\n',
  'mc.py': 'simulation file for Metropolis Monte Carlo\n',
  'opt.py': 'optimization file for energy minimization\n',
  'ana.py': 'plot file for data analysis\n'}

def _GetFileStringArray(infile_name):
  """Create a 2-d array of strings from input file name.

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
  infile = open(infile_name, 'r')
  infile_data = infile.readlines()
  infile.close()

  infile_array = []
  for line in infile_data:
    infile_array.append(line.split())
  return infile_array


def GetElement(at_type):
  """Infer atomic element from atom type.

  If atom type is a single character, or second character is uppercase,
  return uppercase first letter. Otherwise, return capitalized first two
  characters.

  Args:
    at_type (str): Atom type.

  Returns:
    at_element (str): Atomic element.
  """
  if len(at_type) == 1 or not at_type[1].islower():
    return at_type[0].upper()
  else:
    return at_type[0:2].capitalize()


def GetGeom(mol):
  """Read in molecular geometry data from molecule xyzq file.
  
  Parse 2-d array of strings from xyzq file into atomic data. First line
  contains (int) number of atoms. Second line is ignored comment. Each line
  after (3 to [n+2]) contains atom type, (float) 3 xyz cartesian coordinates
  [Angstrom], and (float) charge [e].
  
  Args:
    mol (mmlib.molecule.Molecule): molecule with an associated xyzq input file.
  """
  infile_array = _GetFileStringArray(mol.infile)
  mol.n_atoms = int(infile_array[0][0])
  for i in range(mol.n_atoms):
    at_type = infile_array[i+2][0]
    at_coords = numpy.array(
        list(map(float, infile_array[i+2][1:1+const.NUMDIM])))
    at_charge = float(infile_array[i+2][4])

    at_element = GetElement(at_type)
    at_mass = param.GetAtMass(at_element)
    at_ro, at_eps = param.GetVdwParam(at_type)
    atom = molecule.Atom(at_type, at_coords, at_charge, at_ro, at_eps, at_mass)
    atom.SetCovRad(param.GetCovRad(at_element))
    mol.atoms.append(atom)
    mol.mass += at_mass


def _GetAtom(mol, record):
  """Parse atom record into an atom object and append to molecule.

  Appends mmlib.molecule.Atom object to mmlib.molecule.Molecule object. Contents
  of atom object include (float*) xyz cartesian coordinates [Angstrom], (float)
  partial charge [e], (float) van der Waals radius [Angstrom], (float) van der
  Waals epsilon [kcal/mol], (str) atom type, (str) atomic element, (float)
  covalent radius, [Angstrom], and (float) mass [amu].

  Args:
    mol (mmlib.molecule.Molecule): Molecule to append atom.
    record (str*): Array of strings from line of prm file.
  """
  at_type = record[2]
  at_coords = numpy.array(list(map(float, record[3:3+const.NUMDIM])))
  at_charge, at_ro, at_eps = list(map(float, record[6:9]))

  at_element = GetElement(at_type)
  at_mass = param.GetAtMass(at_element)
  atom = molecule.Atom(at_type, at_coords, at_charge, at_ro, at_eps, at_mass)
  atom.SetCovRad(param.GetCovRad(at_element))
  mol.atoms.append(atom)
  mol.mass += at_mass


def _GetBond(mol, record):
  """Parse bond record into a bond object and append to molecule.
  
  Appends mmlib.molecule.Bond object to mmlib.molecule.Molecule object. Contents
  of bond object include (int) 2 atomic indices, (float) spring constant
  [kcal/(mol*A^2)], (float) equilibrium bond length [Angstrom], (float) bond
  length [Angstrom].
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule to append bond.
    record (str*): Array of strings from line of prm file.
  """
  at1, at2 = [x-1 for x in list(map(int, record[1:3]))]
  k_b, r_eq = list(map(float, record[3:5]))
  c1, c2 = [mol.atoms[i].coords for i in [at1, at2]]

  r_ij = geomcalc.GetRij(c1, c2)
  bond = molecule.Bond(at1, at2, r_ij, r_eq, k_b)
  mol.bonds.append(bond)
  mol.bond_graph[at1][at2] = r_ij
  mol.bond_graph[at2][at1] = r_ij


def _GetAngle(mol, record):
  """Parse angle record into an angle object and append to molecule.
  
  Appends mmlib.molecule.Angle object to mmlib.molecule.Molecule object.
  Contents of angle object include (int) 3 atomic indices, (float) spring
  constant [kcal/(mol*radian^2)], (float) equilibrium bond angle [degrees], and
  (float) bond angle [degrees].
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule to append angle.
    record (str*): Array of strings from line of prm file.
  """
  at1, at2, at3 = [x-1 for x in list(map(int, record[1:4]))]
  k_a, a_eq = list(map(float, record[4:6]))
  c1, c2, c3 = [mol.atoms[i].coords for i in [at1, at2, at3]]

  a_ijk = geomcalc.GetAijk(c1, c2, c3)
  angle = molecule.Angle(at1, at2, at3, a_ijk, a_eq, k_a)
  mol.angles.append(angle)


def _GetTorsion(mol, record):
  """Parse torsion record into a torsion object and append to molecule.
  
  Appends mmlib.molecule.Torsion object to mmlib.molecule.Molecule object.
  Contents of torsion object include (int) 4 atomic indices, (float)
  half-barrier height [kcal/mol], (float) barrier offset [degrees], (int)
  barrier frequency, (int) barrier paths, and (float) torsion angle [degrees].
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule to append torsion.
    record (str*): Array of strings from line of prm file.
  """
  at1, at2, at3, at4 = [x-1 for x in list(map(int, record[1:5]))]
  v_n, gamma = list(map(float, record[5:7]))
  nfold, paths = list(map(int, record[7:9]))
  c1, c2, c3, c4 = [mol.atoms[i].coords for i in [at1, at2, at3, at4]]

  t_ijkl = geomcalc.GetTijkl(c1, c2, c3, c4)
  tors = molecule.Torsion(at1, at2, at3, at4, t_ijkl, v_n, gamma, nfold, paths)
  mol.torsions.append(tors)


def _GetOutofplane(mol, record):
  """Parse outofplane record into object and append to molecule.
  
  Appends mmlib.molecule.Outofplane object to mmlib.molecule.Molecule object.
  Contents of outofplane object include (int) 4 atomic indices, (float) barrier
  height [kcal/mol], and (float) outofplane angle [degrees].
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule to append outofplane.
    record (str*): Array of string from line of prm file.
  """
  at1, at2, at3, at4 = [x-1 for x in list(map(int, record[1:5]))]
  v_n = float(record[5])
  c1, c2, c3, c4 = [mol.atoms[i].coords for i in [at1, at2, at3, at4]]

  o_ijkl = geomcalc.GetOijkl(c1, c2, c3, c4) 
  outofplane = molecule.Outofplane(at1, at2, at3, at4, o_ijkl, v_n)
  mol.outofplanes.append(outofplane)


def GetPrm(mol):
  """Parse contents of prm file into molecular topology / geometry data.
  
  Reads in and organizes contents of a prm parameter file into given
  mmlib.molecule.Molecule object. Contents of molecule object include Atom
  objects, Bond objects, Angle objects, Torsion objects, and Outofplane objects
  and associated parameters and values.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule to append data
  """
  infile_array = _GetFileStringArray(mol.infile)
  for i in range(len(infile_array)):
    record = infile_array[i]
    rec_type = record[0].lower()
    if rec_type == 'atom':
      _GetAtom(mol, record)
  mol.n_atoms = len(mol.atoms)

  parse_functions = {
    'bond': _GetBond,
    'angle': _GetAngle,
    'torsion': _GetTorsion,
    'outofplane': _GetOutofplane}

  mol.bond_graph = {i:{} for i in range(mol.n_atoms)}
  for i in range(len(infile_array)):
    record = infile_array[i]
    record_type = record[0].lower()
    if record_type in parse_functions:
      parse_function = parse_functions[record_type]
      parse_function(mol, record)

  mol.n_bonds = len(mol.bonds)
  mol.n_angles = len(mol.angles)
  mol.n_torsions = len(mol.torsions)
  mol.n_outofplanes = len(mol.outofplanes)
  topology.GetNonints(mol)


def GetSimData(sim):
  """Parse contents of sim file into molecular simulation data.
  
  Many molecular simulation parameters (dynamics, monte carlo, etc.) can be
  determined by default, or overriden in a simulation file. All listed values
  below can be set through the given keyword arguments.
  
  Mandatory values include (str) molecule [file path], (float) temperature
  [Kelvin], and (float / int) total time/confs [ps / none] (md / mc).
  
  Args:
    sim (mmlib.simulate.Simulation): Simulation object to append data.
  """
  infile_array = _GetFileStringArray(sim.infile)
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
      sim.temp = float(kwargval)
    elif kwarg == 'pressure':
      sim.press = float(kwargval)
    elif kwarg == 'boundaryspring':
      sim.mol.k_box = float(kwargval)
    elif kwarg == 'boundary':
      sim.mol.bound = float(kwargval)
      sim.mol.GetVolume()
    elif kwarg == 'boundarytype':
      sim.mol.boundtype = kwargval.lower()
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
  """Parse contents of sim file into molecular simulation data.
  
  Many molecular energy minimization parameters can be determined by default, or
  overridden in an optimization file. All listed values below can be set through
  the give keyword arguments.

  The only mandatory value is (str) molecule [file path]. Setting (str) geomout
  and (str) energyout is also strongly recommended.

  Args:
    opt (mmlib.optimize.Optimization): Optimization object to append data.
  """
  infile_array = _GetFileStringArray(opt.infile)
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
  """Parse contents of plt file into ensemble analysis data.
  
  Many simulation data analysis keywords can be determined by default, or
  overridden in a plot file. All listed values below can be set through the
  given keyword arguments.

  Mandatory values include (str) 'input' [file path] and (str) 'simtype'.
  Setting 'plotout' is also strongly recommended.
  
  Args:
    ana (mmlib.analyze.Analysis): Analysis object to append data.
  """
  infile_array = _GetFileStringArray(ana.infile)
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
  """Read in molecular property sequences from simulation data file.
  
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
  prop_array = _GetFileStringArray(prop_file)
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
  """Read in molecular xyz coordinate sequences from xyz file.
  
  Input file contains 'n_confs' molecular xyz-coordinate snapshots with
  'n_atoms' atoms each. Start by reading in 'n_atoms', and continue until file
  ends.
  
  Args:
    traj_file (str): Path to input trajectory file.
  
  Returns:
    traj (float***): Array of xyz-coordinates at each configuration of molecule
        during trajectory.
  """
  traj_array = _GetFileStringArray(traj_file)
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


def PrintCoords(mol, comment, ofile=None):
  """Print atomic coordinates for a set of atoms.
  
  Print to screen all (float) 3N atomic cartesian coordinates [Angstrom] from
  mol in xyz file format with (str) 'comment' for the comment line.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with (float) 3N atomic
        cartesian coordinates [Angstrom].
    comment (str): Comment string for xyz file comment line.
  """
  print('%i\n%s\n' % (mol.n_atoms, comment), end='')
  for i in range(mol.n_atoms):
    print('%-2s' % (mol.atoms[i].type_), end='')
    for j in range(const.NUMDIM):
      print(' %12.6f' % (mol.atoms[i].coords[j]), end='')
    print('\n', end='')
  print('\n', end='')


def PrintCoordsFile(ofile, mol, comment, totchar, decchar):
  """Print atomic coordinates to an open output file.
  
  Print to ofile all (float) 3N atomic cartesian coordinates [Angstrom] from mol
  in xyz file format with (str) 'comment' in comment line.
  
  Args:
    ofile (file): Open file handle for printing coordinate data.
    mol (mmlib.molecule.Molecule): Molecule object with (float) 3N atomic
        cartesian coordinates [Angstrom].
    comment (str): Comment string for xyz file comment line.
    totchar (str): Total number of characters in coordinate printing.
    decchar (str): Post-decimal characters in coordinate printing.
  """
  ofile.write('%i\n%s\n' % (mol.n_atoms, comment))
  for i in range(mol.n_atoms):
    ofile.write('%-2s' % (mol.atoms[i].element))
    for j in range(const.NUMDIM):
      ofile.write(' %*.*f' % (totchar, decchar, mol.atoms[i].coords[j]))
    ofile.write('\n')


def PrintGradient(grad, comment):
  """Print specified atomic gradient type for a set of atoms.
  
  Print to screen all (float) 3N atomic cartesian gradient components
  [kcal/(mol*Angstrom)] from mol in xyz file format. Gradient is partial
  derivative of 'grad_type' energy [kcal/mol] with respect to each (float)
  cartesian coordinate [Angstrom].
  
  Args:
    grad (numpy.float**): Molecular energy gradient (or component)
        [kcal/(mol*A)].
    comment (str): Comment on gradient type / source / etc. to print to screen.
  """
  print('\n %s\n' % (comment))
  for i in range(mol.n_atoms):
    print('%-2s' % (mol.atoms[i].type_), end='')
    for j in range(const.NUMDIM):
      print(' %12.6f' % (grad[i][j]), end='')
    print('\n', end='')
  print('\n', end='')


def _PrintBanner(title, length, lines1, lines2):
  """Print a string in the center of a banner with dashes on each side
  
  Print leading newlines, one space, dashes to center string, header string,
  dashes to end of line, and trailing newlines to screen.
  
  Args:
    title (str): Banner header title.
    length (int): Total number of characters in banner.
    lines1 (int): Number of leading newlines.
    lines2 (int): Number of trailing newlines.
  """
  dashes1 = (length - len(title))//2 - 1
  dashes2 = length - len(title) - dashes1 - 2
  return lines1*'\n' + dashes1*'-' + title + dashes2*'-' + lines2*'\n'


def _PrintPadded(fields, spacings):
  """Print an array of strings padded by spaces.
  
  Print leading number of spaces prior to each string, each elements of the
  arrays 'strings' and 'spacings', respectively.
  
  Args:
    fields (str*): Array of string fields to be printed.
    spacings (int*): Array of number of spaces to be printed prior to each
        'strings' element.
  """
  string = []
  for spacing, field in zip(spacings, fields):
    string.append(spacing*' ' + field)
  return ''.join(string) + '\n'


def _PrintHeader(header, n_banner, params, spacings):
  """Print banner header for a section of output.
  
  Print a header string in the center of a dash banner, followed by a set of
  column headers, and a trailing line of dashes.
  
  Args:
    header (str): Banner header title.
    n_banner (int): Total number of characters in banner.
    params (str*): Array of strings to be printed.
    spacings (int*): Array of number of spaces to be printed prior to each
        'params' element.
  """
  string = ''.join([
    _PrintBanner(header, n_banner, 1, 1),
    _PrintPadded(params, spacings),
    _PrintBanner('', n_banner, 0, 0)])
  print(string)


def PrintGeom(mol):
  """Print geometry and non-bonded parameters for a molecule to screen.
  
  Print a header banner for the section, and for each Atom object in a Molecule
  object print the (str) atom type, (float) 3 xyz cartesian coordinates
  [Angstrom], (float) partial charge [e], (float) van der waals radius
  [Angstrom], and (float) van der waals epsilon [kcal/mol].
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule with Atom objects with data for
        printing.
  """
  _PrintHeader(_GEOM_PRINT_HEADER, _GEOM_PRINT_BANNER_CHARS, _GEOM_PRINT_PARAMS,
               _GEOM_PRINT_SPACES)
  
  for i in range(mol.n_atoms):
    print('%4i | %-2s' % (i+1, mol.atoms[i].type_), end='')
    for j in range(const.NUMDIM):
      print('%10.4f' % (mol.atoms[i].coords[j]), end='')
    print(' %7.4f %7.4f %7.4f' % (mol.atoms[i].charge,
      mol.atoms[i].ro, mol.atoms[i].eps))


def PrintGeomFile(outfile, mol):
  """Write out molecular geometry and non-bonded parameters to prm file.
  
  For each Atom object in Molecule object, write to outfile an atom record,
  containing (int) atomic index, (str) atom type, (float) 3 cartesian
  coordinates [Angstrom], (float) partial charge [e], (float) van der waals
  radius [Angstrom], and (float) van der waals epsilon [kcal/mol].
  
  Args:
    outfile (_io.TextIOWrapper): Output stream to open prm file.
    mol (mmlib.molecule.Molecule): Molecule with geometry and non-bonded
        parameter data for printing.
  """
  outfile.write('# %s Atoms (at, type, x, y, z, q, ro, eps)\n' % (mol.n_atoms))
  for i in range(mol.n_atoms):
    outfile.write('ATOM %4i %-2s' % (i+1, mol.atoms[i].type_))
    for j in range(const.NUMDIM):
      outfile.write(' %11.6f' % (mol.atoms[i].coords[j]))
    outfile.write(' %8.5f %7.4f %7.4f\n' % (
        mol.atoms[i].charge, mol.atoms[i].ro, mol.atoms[i].eps))


def PrintBonds(mol):
  """Print bond topology and parameters for a molecule to screen.
  
  Print a header banner for the section, and for each Bond object in a Molecule
  object print (float) spring constant [kcal/(mol*A^2)], (float) equilibrium
  bond length [Angstrom], (float) bond length [Angstrom], (str) 2 atom types,
  (float) bond energy [kcal/mol], and (int) 2 atomic indices.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule with Bond objects with data for
        printing.
  """
  if mol.n_bonds:
    _PrintHeader(_BOND_PRINT_HEADER, _BOND_PRINT_BANNER_CHARS,
                 _BOND_PRINT_PARAMS, _BOND_PRINT_SPACES)
  else:
    print(_BOND_PRINT_ABSENT)
  
  a = mol.atoms
  for p in range(mol.n_bonds):
    b = mol.bonds[p]
    t1, t2 = a[b.at1].type_, a[b.at2].type_
    pstr = '%4i | %7.2f %8.4f %8.4f (%2s-%2s) %8.4f (%i-%i)' % (
        p+1, b.k_b, b.r_eq, b.r_ij, t1, t2, b.energy, b.at1+1, b.at2+1)
    print(pstr)


def PrintBondsFile(outfile, mol):
  """Write out molecular bond data and parameters to prm file.
  
  For each Bond object in Molecule object, write to outfile a bond record,
  containing (int) 2 atomic indices, (float) spring constant [kcal/(mol*A^2)],
  and (float) equilibrium bond length [Angstrom].
  
  Args:
    outfile (_io.TextIOWrapper): Output stream to open prm file.
    mol (mmlib.molecule.Molecule): Molecule with bond data and parameters for
        printing.
  """
  outfile.write('# %i Bonds (At1, At2, K_b, R_eq)\n' % (mol.n_bonds))
  for p in range(mol.n_bonds):
    b = mol.bonds[p]
    outfile.write('BOND %4i %4i %7.2f %7.4f\n' % (
        b.at1+1, b.at2+1, b.k_b, b.r_eq))


def PrintAngles(mol):
  """Print angle topology and parameters for a molecule to screen.
  
  Print a header banner for the section, and for each Angle object in a Molecule
  object print (float) spring constant [kcal/(mol*rad^2)], (float) equilibrium
  bond angle [degrees], (float) bond angle [degrees], (str) 3 atom types,
  (float) angle energy [kcal/mol], and (int) 3 atomic indices.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule with Angle objects with data for
        printing.
  """
  if mol.n_angles:
    _PrintHeader(_ANGLE_PRINT_HEADER, _ANGLE_PRINT_BANNER_CHARS,
                 _ANGLE_PRINT_PARAMS, _ANGLE_PRINT_SPACES)
  else:
    print(_ANGLE_PRINT_ABSENT)
  
  at = mol.atoms
  for p in range(mol.n_angles):
    a = mol.angles[p]
    t1, t2, t3 = at[a.at1].type_, at[a.at2].type_, at[a.at3].type_
    pstr = '%4i | %6.2f %7.3f %7.3f (%2s-%2s-%2s) %7.4f (%i-%i-%i)' % (p+1,
        a.k_a, a.a_eq, a.a_ijk, t1, t2, t3, a.energy, a.at1+1, a.at2+1, a.at3+1)
    print(pstr)


def PrintAnglesFile(outfile, mol):
  """Write out molecular angle data and parameters to prm file.
  
  For each Angle object in Molecule object, write to outfile an angle record,
  containing (int) 3 atomic indices, (float) spring constant [kcal/(mol*rad^2)],
  and (float) equilibrium bond angle [degrees].
  
  Args:
    outfile (_io.TextIOWrapper): Output stream to open prm file.
    mol (mmlib.molecule.Molecule): Molecule with angle data and parameters for
        printing.
  """
  outfile.write('# %i Angles (At1, At2, At3, K_a, A_eq)\n' % (mol.n_angles))
  for p in range(mol.n_angles):
    a = mol.angles[p]
    outfile.write('ANGLE %4i %4i %4i %7.4f %8.4f\n' % (
        a.at1+1, a.at2+1, a.at3+1, a.k_a, a.a_eq))


def PrintTorsions(mol):
  """Print torsion topology and parameters for a molecule to screen.
  
  Print a header banner for the section, and for each Torsion object in a
  Molecule object print (float) barrier height [kcal/mol], (float) barrier
  offset [degrees], (int) barrier frequency, (int) barrier paths, (str) 4 atom
  types, (float) torsion energy [kcal/mol], and (int) 4 atomic indices.

  Args:
    mol (mmlib.molecule.Molecule): Molecule with Torsion objects with data for
        printing.
  """
  if mol.n_torsions:
    _PrintHeader(_TORSION_PRINT_HEADER, _TORSION_PRINT_BANNER_CHARS,
                 _TORSION_PRINT_PARAMS, _TORSION_PRINT_SPACES)
  else:
    print(_TORSION_PRINT_ABSENT)
  
  a = mol.atoms
  for p in range(mol.n_torsions):
    t = mol.torsions[p]
    t1, t2 = a[t.at1].type_, a[t.at2].type_
    t3, t4 = a[t.at3].type_, a[t.at4].type_
    print('%4i | %6.2f %6.1f %8.3f %i %i (%2s-%2s-%2s-%2s) %7.4f (%i-%i-%i-%i)'
          % (p+1, t.v_n, t.gam, t.t_ijkl, t.n, t.paths, t1, t2, t3, t4,
             t.energy, t.at1+1, t.at2+1, t.at3+1, t.at4+1))


def PrintTorsionsFile(outfile, mol):
  """Write out molecular torsion data and parameters to prm file.
  
  For each Torsion object in Molecule object, write to outfile a torsion record,
  ontaining (int) 4 atomic indices, (float) half-barrier height [kcal/mol],
  (float) barrier offset [degrees], (int) barrier frequency, and (int) barrier
  paths.
  
  Args:
    outfile (_io.TextIOWrapper): Output stream to open prm file.
    mol (mmlib.molecule.Molecule): Molecule with torsion data and parameters for
        printing.
  """
  outfile.write('# %i Torsions (At1, At2, At3, At4, '
                'V_n, Gamma, N_f, paths)\n' % (mol.n_torsions))
  for p in range(mol.n_torsions):
    t = mol.torsions[p]
    outfile.write('TORSION %4i %4i %4i %4i %6.2f %6.1f %i %i\n' % (
        t.at1+1, t.at2+1, t.at3+1, t.at4+1, t.v_n, t.gam, t.n, t.paths))


def PrintOutofplanes(mol):
  """Print outofplane topology and parameters for a molecule to screen.
  
  Print a header banner for the section, and for each Outofplane object in a
  Molecule object print (float) half-barrier height [kcal/mol], (float)
  outofplane angle [degrees], (str) 4 atom types, (float) torsion energy
  [kcal/mol], and (int) 4 atomic indices.

  Args:
    mol (mmlib.molecule.Molecule): Molecule with Outofplane objects with data
        for printing.
  """
  if mol.n_outofplanes:
    _PrintHeader(_OUTOFPLANE_PRINT_HEADER, _OUTOFPLANE_PRINT_BANNER_CHARS,
                 _OUTOFPLANE_PRINT_PARAMS, _OUTOFPLANE_PRINT_SPACES)
  else:
    print(_OUTOFPLANE_PRINT_ABSENT)
  
  a = mol.atoms
  for p in range(mol.n_outofplanes):
    o = mol.outofplanes[p]
    t1, t2 = a[o.at1].type_, a[o.at2].type_
    t3, t4 = a[o.at3].type_, a[o.at4].type_
    print('%4i | %6.2f %7.3f (%2s-%2s-%2s-%2s) %7.4f (%i-%i-%i-%i)'
          % (p+1, o.v_n, o.o_ijkl, t1, t2, t3, t4, o.energy, o.at1+1, o.at2+1,
             o.at3+1, o.at4+1))


def PrintOutofplanesFile(outfile, mol):
  """Write out molecular outofplane data and parameters to prm file.
  
  For each Outofplane object in Molecule object, write to outfile an outofplane
  record, containing (int) 4 atomic indices, and (float) half-barrier height
  [kcal/mol].
  
  Args:
    outfile (_io.TextIOWrapper): Output stream to open prm file.
    mol (mmlib.molecule.Molecule): Molecule with outofplane data and parameters
        for printing.
  """
  outfile.write('# %i Outofplanes (At1, At2, At3, At4, V_n, Gamma, N_f)\n'
                % mol.n_outofplanes)
  for p in range(mol.n_outofplanes):
    o = mol.outofplanes[p]
    outfile.write('OUTOFPLANE %4i %4i %4i %4i %6.2f\n'
                  % (o.at1+1, o.at2+1, o.at3+1, o.at4+1, o.v_n))


def PrintEnergy(mol):
  """Print list of energy values in a table to screen.
  
  For each energy term in 'labels' print one (float) value [kcal/mol] to screen
  per line.
  
  Args:
    mol (mmlib.molecule.Molecule): Molecule object with energy component data.
  """
  _PrintHeader(_ENERGY_PRINT_HEADER, _ENERGY_PRINT_BANNER_CHARS,
               _ENERGY_PRINT_PARAMS, _ENERGY_PRINT_SPACES)
  
  for label, attribute in zip(_ENERGY_PRINT_LABELS, _ENERGY_PRINT_ATTRIBUTES):
    print('   %-13s | %10.4f' % (label, getattr(mol, attribute)))


def PrintAverages(ana):
  """Print list of expectation values in a table to screen.
  
  For each energy term in mmlib.analyze.Analysis.pdict print (float) average,
  standard deviation, minimum, and maximum value [kcal/mol] to screen.

  Args:
    ana (mmlib.analyze.Analyze): Analyze object with energy component
        expectation values.
  """
  _PrintHeader(_AVERAGE_PRINT_HEADER, _AVERAGE_PRINT_BANNER_CHARS,
               _AVERAGE_PRINT_PARAMS, _AVERAGE_PRINT_SPACES)
  
  pdict = const.PROPERTYDICTIONARY
  keys = sorted(list(pdict.keys()), key = lambda x: pdict[x][3])
  keys = [key for key in keys if key in ana.prop]
  labels = [pdict[key][0] for key in keys]
  for key, label in zip(keys, labels):
    print('   %-13s | %11.4e %11.4e %11.4e %11.4e' % (
        label, ana.eavg[key], ana.estd[key], ana.emin[key], ana.emax[key]))


def ValidateInput(program_name):
  """Checks for proper input argument syntax and return parsed result.
  
  Check that the command line input contains two strings or throw an error and
  print usage guidance. If correct, return the name of the input file given as
  the second command line input string.

  Args:
    program_name (str): Name of main file calling this function.

  Returns:
    infile_name (str): Name of input file given from command line.

  Raises:
    ValueError: If program_name not recognized.
  """
  if (len(sys.argv) < 2):
    if program_name in _PROGRAM_MESSAGES:
      print('\nUsage: python %s INPUT_FILE\n' % program_name)
      print('INPUT_FILE: %s' % _PROGRAM_MESSAGES[program_name])
      sys.exit()
    else:
      raise ValueError('Program name not recognized: %s' % program_name)
  
  input_file = sys.argv[1]
  if os.path.isfile(input_file):
    return input_file
  else:
    raise ValueError('Specified input is not a file: %s' % input_file)
