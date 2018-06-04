"""Classes and functions for handling molecular system data.

Includes classes for representing geometry and parameter data for atoms, bonds,
bond angles, torsion angles, outofplane angles, and entire molecular systems.
"""

import math
import numpy
import os

from mmlib import constants as const
from mmlib import energy
from mmlib import fileio
from mmlib import geomcalc
from mmlib import gradient
from mmlib import param
from mmlib import topology

class Atom:
  """Atom class for atomic geometry and parameter data.
  
  Initialize attributes to corresponding specified argument values, look up in
  parameter tables, or set to zero.
  
  Args / Attributes:
    type_ (str): AMBER94 mm atom type.
    coords (float*): 3 cartesian coordinates [Angstrom].
    charge (float): atomic partial charge [e].
    ro (float): vdw radius [Angstrom].
    eps (float): vdw epsilon [kcal/mol].
    mass (float): atomic mass [g/mol].

  Attributes:
    covrad (float): covalent radius [Angstrom].
    sreps (float): square root of vdw epsilon [(kcal/mol)^0.5].
    vels (float*): 3 cartesian velocity components [Angstrom/ps].
    accs (float*): 3 cartesian acceleration components [A/(ps^2)].
    pvels (float*): 3 previous 'vels' [Angstrom/ps].
    paccs (float*): 3 previous 'accs' [Angstrom/(ps^2)].
  """
  def __init__(self, type_, coords, charge, ro=None, eps=None):
    self.SetType(type_)
    self.SetCoords(coords)
    self.SetCharge(charge)

    if ro == None or eps == None:
      ro, eps = param.GetVdwParam(self.type_)
    self.SetRo(ro)
    self.SetEps(eps)

    self.SetElement(param.GetElement(type_))
    self.SetMass(param.GetMass(self.element))
    self.SetCovRad(param.GetCovRad(self.element))

    self.SetVels(numpy.zeros(const.NUMDIM))
    self.SetAccs(numpy.zeros(const.NUMDIM))
    self.SetPVels(numpy.zeros(const.NUMDIM))
    self.SetPAccs(numpy.zeros(const.NUMDIM))

  def SetType(self, type_):
    """Set new (str) atom type."""
    self.type_ = type_

  def SetCoords(self, coords):
    """Set new (float*) coodinates [Angstrom]."""
    self.coords = coords

  def SetCoord(self, index, coord):
    """Set new (float) ith coordinate [Angstrom]."""
    self.coords[index] = coords

  def SetCharge(self, charge):
    """Set new (float) partial charge [e]."""
    self.charge = charge

  def SetRo(self, ro):
    """Set new (float) vdw radius [Angstrom]."""
    self.ro = ro

  def SetEps(self, eps):
    """Set new (float) vdw epsilon [kcal/mol]."""
    self.eps = eps
    self.sreps = math.sqrt(eps)

  def SetElement(self, element):
    """Set new (str) atomic element."""
    self.element = element

  def SetMass(self, mass):
    """Set new (float) atomic mass [g/mol]."""
    self.mass = mass

  def SetCovRad(self, covrad):
    """Set new (float) covalent radius [Angstrom]."""
    self.covrad = covrad

  def SetVels(self, vels):
    """Set new (float*) velocities [Angstrom/ps]."""
    self.vels = vels

  def SetAccs(self, accs):
    """Set new (float*) accelerations [Angstrom/(ps^2)]."""
    self.accs = accs

  def SetPVels(self, pvels):
    """Set new (float*) velocities [Angstrom/ps]."""
    self.pvels = pvels

  def SetPAccs(self, paccs):
    """Set new (float*) accelerations [Angstrom/(ps^2)]."""
    self.paccs = paccs


class Bond:
  """Bond class for bond geometry and parameter data.
  
  Initialize attributes to specified argument values. Change by calling
  appropriate 'Set[Param]' function.
  
  Args / Attributes:
    at1 (int): Atom1 atomic index in Molecule.
    at2 (int): Atom2 atomic index in Molecule.
    r_ij (float): Distance [Angstrom] between at1 and at2.
    r_eq (float): Equlibrium at1-at2 bond length [Angstrom].
    k_b (float): Bond spring constant [kcal/(mol*A^2)].

  Attributes:
    energy (float): Energy of bond [kcal/mol].
    grad_mag (float): Energy gradient magnitude of bond [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, r_ij, r_eq, k_b):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetRij(r_ij)
    self.SetReq(r_eq)
    self.SetKb(k_b)

  def SetAt1(self, at1):
    """Set new (int) atomic index 1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) atomice index 2."""
    self.at2 = at2

  def SetRij(self, r_ij):
    """Set new (float)  [Angstrom]."""
    self.r_ij = r_ij

  def SetReq(self, r_eq):
    """Set new (float)  [Angstrom]."""
    self.r_eq = r_eq

  def SetKb(self, k_b):
    """Set new (float)  [kcal/(mol*A^2)]."""
    self.k_b = k_b

  def GetEnergy(self):
    """Calculate bond energy (float) [kcal/mol]."""
    self.energy = energy.GetEBond(self.r_ij, self.r_eq, self.k_b)

  def GetGradientMagnitude(self):
    """Calculate bond gradient magnitude (float) [kcal/(mol*A)]."""
    self.grad_mag = gradient.GetGMagBond(self.r_ij, self.r_eq, self.k_b)


class Angle:
  """Angle class for angle geometry and parameter data.
  
  Initialize attributes to specified argument values. Change by calling
  appropriate 'Set[Param]' function.
  
  Args / Attributes:
    at1 (int): Atom1 atomic index in Molecule.
    at2 (int): Atom2 atomic index in Molecule.
    at3 (int): Atom3 atomic index in Molecule.
    a_ijk (float): Current at1-at2-at3 bond angle [degrees].
    a_eq (float): Equlibrium at1-at2-at3 bond angle [degrees].
    k_a (float): Angle spring constant [kcal/(mol*rad^2)].

  Attributes:
    energy (float): Energy of angle [kcal/mol].
    grad_mag (float): Energy gradient magnitude of angle [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, at3, a_ijk, a_eq, k_a):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetAt3(at3)
    self.SetAijk(a_ijk)
    self.SetAeq(a_eq)
    self.SetKa(k_a)

  def SetAt1(self, at1):
    """Set new (int) atomic index 1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) atomic index 2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) atomic index 3."""
    self.at3 = at3

  def SetAijk(self, a_ijk):
    """Set new (float) a_ijk [degrees]."""
    self.a_ijk = a_ijk

  def SetAeq(self, a_eq):
    """Set new (float) a_eq [degrees]."""
    self.a_eq = a_eq

  def SetKa(self, k_a):
    """Set new (float) k_a [kcal/(mol*rad^2)]."""
    self.k_a = k_a

  def GetEnergy(self):
    """Get energy (float) [kcal/mol]."""
    self.energy = energy.GetEAngle(self.a_ijk, self.a_eq, self.k_a)

  def GetGradientMagnitude(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.grad_mag = gradient.GetGMagAngle(self.a_ijk, self.a_eq, self.k_a)


class Torsion:
  """Torsion class for torsion geometry and parameter data.
  
  Initialize attributes to specified argument values. Change
  by calling appropriate 'set_[param]' function.
  
  Args / Attributes:
    at1 (int): Atom1 atomic index in Molecule.
    at2 (int): Atom2 atomic index in Molecule.
    at3 (int): Atom3 atomic index in Molecule.
    at4 (int): Atom4 atomic index in Molecule.
    t_ijkl (float): Current at1-at2-at3-at4 torsion angle [degrees].
    v_n (float): Rotation half-barrier height [kcal/mol].
    gamma (float): Barrier offset [degrees].
    nfold (int): Frequency of barrier.
    paths (int): Unique paths through torsion.

  Attributes:
    energy (float): Energy of torsion [kcal/mol].
    grad (float): Energy gradient magnitude of torsion [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, at3, at4, t_ijkl, v_n, gamma, nfold, paths):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetAt3(at3)
    self.SetAt4(at4)
    self.SetTijkl(t_ijkl)
    self.SetVn(v_n)
    self.SetGamma(gamma)
    self.SetNfold(nfold)
    self.SetPaths(paths)

  def SetAt1(self, at1):
    """Set new (int) atomic index 1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) atomic index 2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) atomic index 3."""
    self.at3 = at3

  def SetAt4(self, at4):
    """Set new (int) atomic index 4."""
    self.at4 = at4

  def SetTijkl(self, t_ijkl):
    """Set new (float) t_ijkl [degrees]."""
    self.t_ijkl = t_ijkl

  def SetVn(self, v_n):
    """Set new (float) v_n [kcal/mol]."""
    self.v_n = v_n

  def SetGamma(self, gamma):
    """Set new (float) gamma [degrees]."""
    self.gam = gamma

  def SetNfold(self, nfold):
    """Set new (int) nfold."""
    self.n = nfold

  def SetPaths(self, paths):
    """Set new (int) paths."""
    self.paths = paths

  def GetEnergy(self):
    """Get energy (float) [kcal/mol]."""
    self.energy = energy.GetETorsion(self.t_ijkl, self.v_n, self.gam, self.n,
                                     self.paths)

  def GetGradientMagnitude(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.grad_mag = gradient.GetGMagTorsion(self.t_ijkl, self.v_n, self.gam,
                                            self.n, self.paths)


class Outofplane:
  """Outofplane class for outofplane geometry and parameter data.
  
  Initialize attributes to specified argument values. Change by calling
  appropriate 'set_[param]' function.
  
  Args / Attributes:
    at1 (int): Atom1 atomic index in Molecule.
    at2 (int): Atom2 atomic index in Molecule.
    at3 (int): Atom3 atomic index in Molecule.
    at4 (int): Atom4 atomic index in Molecule.
    o_ijkl (float): Current at1-at2-at3-at4 outofplane angle [degrees].
    v_n (float): Half-barrier height [kcal/mol].

  Attributes:
    energy (float): Energy of outofplane [kcal/mol].
    grad (float): Energy gradient magnitude of outofplane [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, at3, at4, o_ijkl, v_n):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetAt3(at3)
    self.SetAt4(at4)
    self.SetOijkl(o_ijkl)
    self.SetVn(v_n)

  def SetAt1(self, at1):
    """Set new (int) atomic index 1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) atomic index 2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) atomic index 3."""
    self.at3 = at3

  def SetAt4(self, at4):
    """Set new (int) atomic index 4."""
    self.at4 = at4

  def SetOijkl(self, o_ijkl):
    """Set new (float) o_ijkl [degrees]."""
    self.o_ijkl = o_ijkl

  def SetVn(self, v_n):
    """Set new (float) v_n [kcal/mol]."""
    self.v_n = v_n

  def GetEnergy(self):
    """Get energy (float) [kcal/mol]."""
    self.energy = energy.GetEOutofplane(self.o_ijkl, self.v_n)

  def GetGradientMagnitude(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.grad_mag = gradient.GetGMagOutofplane(self.o_ijkl, self.v_n)


class Molecule:
  """
  Molecule class for molecular geometry / topology / energy data.
  
  Contains 'n_atoms' Atom objects, 'n_bonds' Bond objects, 'n_angles' Angle'
  objects, 'n_torsions' Torsion objects, and 'n_outofplanes' Outofplane objects,
  with all their associated data.
  
  Also contains (float) energy and (float**) gradient total values and their
  components: kinetic, potential, bonded, non-bonded, boundary, van der waals,
  electrostatic, bonds, angles, torsions, and outofplanes.
  
  Also contains bulk properties like volume, pressure, temperature, boundaries,
  dielectric, origin and virial.
  
  Args:
    infile_name (str): xyzq or prm input file with molecule data.
  
  Attributes:
    infile (str): Absolute path to 'infile_name'
    indir (str): Absolute path to infile directory.
    filetype (str): Input file format: 'xyzq' or 'prm'.
    name (str): Name of molecule from input file name.

    atoms (mmlib.molecule.Atom*): Array of Atom objects.
    bonds (mmlib.molecule.Bond*): Array of Bond objects.
    angles (mmlib.molecule.Angle*): Array of Atom objects.
    torsions (mmlib.molecule.Torsion*): Array of Torsion objects.
    outofplanes (mmlib.molecule.Outofplane*): Array of Outofplane objects.

    n_atoms (int): Number of atoms.
    n_bonds (int): Number of bonds.
    n_angles (int): Number of angles.
    n_torsions (int): Number of torsions.
    n_outofplanes (int): Number of outofplanes.

    nonints (set(int, int)): Array of covalently bonded atomic indices.
    bond_graph (dict(int:dict(int: float))): Nested dictionary keyed by atom
        pair indices with bond length as value.

    dielectric (float): Dielectric constant. Default = 1.0 (free space).
    mass (float): Sum of all atomic masses.
    k_box (float): Spring constant [kcal/(mol*A^2)] of boundary potential.
    boundary (float): (spherical / cubic) dimensions of system [Angstrom].
    boundary_type (str): Type of boundary shape, 'cube', 'sphere', or 'none'.
    origin (float*): Cartesian center of boundary object.
    volume (float): Volume of system [A^3].
    temperature (float): Instantaneous kinetic temperature [K].
    pressure (float): Instantaneous kinetic pressure [Pa].
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
    self.name = os.path.splitext(os.path.basename(self.infile))[0]

    self.atoms = []
    self.bonds = []
    self.angles = []
    self.torsions = []
    self.outofplanes = []

    self.n_atoms = 0
    self.n_bonds = 0
    self.n_angles = 0
    self.n_torsions = 0
    self.n_outofplanes = 0

    self.nonints = set()
    self.bond_graph = dict()

    self.dielectric = 1.0
    self.mass = 0.0
    self.k_box = 250.0
    self.boundary = 1.0E10
    self.boundary_type = 'sphere'
    self.origin = numpy.zeros(const.NUMDIM)
    self.volume = float('inf')
    self.temperature = 0.0
    self.pressure = 0.0
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
      self.ReadInXYZQ()
      self.GetTopology()
    elif (self.filetype == 'prm'):
      self.ReadInPrm()
        
    self.g_bonds = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_angles = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_torsions = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_outofplanes = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_vdw = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_elst = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_bound = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_bonded = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_nonbonded = numpy.zeros((self.n_atoms, const.NUMDIM))
    self.g_total = numpy.zeros((self.n_atoms, const.NUMDIM))

  def ReadInXYZQ(self):
    """Read in xyzq data from .xyzq input file."""
    input_rows = fileio.GetFileStringArray(self.infile)
    self.atoms = fileio.GetAtomsFromXyzq(input_rows)
    self.n_atoms = len(self.atoms)

  def ReadInPrm(self):
    """Read in prm data from .prm input file."""
    input_rows = fileio.GetFileStringArray(self.infile)

    self.atoms = fileio.GetAtomsFromPrm(input_rows)
    self.bonds = fileio.GetBondsFromPrm(input_rows, self.atoms)
    self.angles = fileio.GetAnglesFromPrm(input_rows, self.atoms)
    self.torsions = fileio.GetTorsionsFromPrm(input_rows, self.atoms)
    self.outofplanes = fileio.GetOutofplanesFromPrm(input_rows, self.atoms)

    self.n_atoms = len(self.atoms)
    self.n_bonds = len(self.bonds)
    self.n_angles = len(self.angles)
    self.n_torsions = len(self.torsions)
    self.n_outofplanes = len(self.outofplanes)

    self.bond_graph = topology.GetBondGraphFromBonds(self.bonds, self.n_atoms)
    self.nonints = topology.GetNonints(self.bonds, self.angles, self.torsions)

  def GetTopology(self):
    """Determine bonded topology of molecules from coordinates."""
    self.bond_graph = topology.GetBondGraph(self.atoms)
    self.bonds = topology.GetBonds(self.atoms, self.bond_graph)
    self.angles = topology.GetAngles(self.atoms, self.bond_graph)
    self.torsions = topology.GetTorsions(self.atoms, self.bond_graph)
    self.outofplanes = topology.GetOutofplanes(self.atoms, self.bond_graph)
    self.nonints = topology.GetNonints(self.bonds, self.angles, self.torsions)

    self.n_bonds = len(self.bonds)
    self.n_angles = len(self.angles)
    self.n_torsions = len(self.torsions)
    self.n_outofplanes = len(self.outofplanes)

  def GetEnergy(self, kintype=None):
    """Calculate (float) energy [kcal/mol] and all energy components."""
    self.e_bonds = energy.GetEBonds(self.bonds)
    self.e_angles = energy.GetEAngles(self.angles)
    self.e_torsions = energy.GetETorsions(self.torsions)
    self.e_outofplanes = energy.GetEOutofplanes(self.outofplanes)
    self.e_vdw, self.e_elst = energy.GetENonbonded(self.atoms, self.nonints,
                                                   self.dielectric)
    self.e_bound = energy.GetEBound(self.atoms, self.k_box, self.boundary,
                                    self.origin, self.boundary_type)
    self.e_kinetic = energy.GetEKinetic(self.atoms, kintype)

    self.e_bonded = (
        self.e_bonds +
        self.e_angles +
        self.e_torsions + 
        self.e_outofplanes)

    self.e_nonbonded = (
        self.e_vdw +
        self.e_elst)

    self.e_potential = (
        self.e_bonded +
        self.e_nonbonded +
        self.e_bound)

    self.e_total = (
        self.e_potential +
        self.e_kinetic)

  def GetGradient(self, grad_type='analytic'):
    """Calculate analytical or numerical gradient of energy.
    
    Args:
      grad_type (str): Type of gradient:
        'analytic': (default) exact, based on analytic derivatives.
        'numerical': approximate, based on numerical derivatives.
    """
    if (grad_type == 'analytic'):
      self.GetAnalyticGradient()
    elif (grad_type == 'numerical'):
      self.GetNumericalGradient()
    else:
      raise ValueError('Unexpected gradient type: %s\n'
                       "Use 'analytic' or 'numerical'." % grad_type)

    self.g_bonded.fill(0.0)
    self.g_bonded += self.g_bonds
    self.g_bonded += self.g_angles
    self.g_bonded += self.g_torsions
    self.g_bonded += self.g_outofplanes

    self.g_nonbonded.fill(0.0)
    self.g_nonbonded += self.g_vdw
    self.g_nonbonded += self.g_elst

    self.g_total.fill(0.0)
    self.g_total += self.g_bonded
    self.g_total += self.g_nonbonded
    self.g_total += self.g_bound

  def GetAnalyticGradient(self):
    """Calculate analytic (float**) gradient [kcal/(mol*A)] of energy."""
    gradient.GetGBonds(self.g_bonds, self.bonds, self.atoms)
    gradient.GetGAngles(self.g_angles, self.angles, self.atoms, self.bond_graph)
    gradient.GetGTorsions(
        self.g_torsions, self.torsions, self.atoms, self.bond_graph)
    gradient.GetGOutofplanes(
        self.g_outofplanes, self.outofplanes, self.atoms, self.bond_graph)

  def GetNumericalGradient(self):
    """Calculate numerical (float**) gradient [kcal/(mol*A)] of energy."""
    gradient.GetGNumerical(self)

  def UpdateInternals(self):
    """Update current values of internal degrees of freedom."""
    topology.UpdateBonds(self.bonds, self.atoms, self.bond_graph)
    topology.UpdateAngles(self.angles, self.atoms, self.bond_graph)
    topology.UpdateTorsions(self.torsions, self.atoms, self.bond_graph)
    topology.UpdateOutofplanes(self.outofplanes, self.atoms, self.bond_graph)

  def GetTemperature(self):
    """Calculate instantaneous kinetic temperature [K] of system."""
    self.temperature = energy.GetTemperature(self.e_kinetic, self.n_atoms)

  def GetPressure(self):
    """Calculate instantaneous kinetic pressure [Pa] of system."""
    self.virial = gradient.GetVirial(self.g_total, self.atoms)
    self.pressure = gradient.GetPressure(
        self.n_atoms, self.temperature, self.virial, self.volume)

  def GetVolume(self):
    """Caclculate approximate volume [A^3] of system."""
    self.volume = geomcalc.GetVolume(self.boundary, self.boundary_type)

  def PrintData(self):
    """Print energy / geometry / topology data of molecule to screen."""
    self.PrintEnergy()
    self.PrintGeom()
    self.PrintBonds()
    self.PrintAngles()
    self.PrintTorsions()
    self.PrintOutofplanes()

  def PrintEnergy(self):
    """Print energy and component data of molecule to screen."""
    print(fileio.GetPrintEnergyString(self))

  def PrintGeom(self):
    """Print geometry data of molecule to screen."""
    print(fileio.GetPrintGeomString(self.atoms))

  def PrintBonds(self):
    """Print bond data of molecule to screen."""
    print(fileio.GetPrintBondsString(self.bonds, self.atoms))

  def PrintAngles(self):
    """Print angle data of molecule to screen."""
    print(fileio.GetPrintAnglesString(self.angles, self.atoms))

  def PrintTorsions(self):
    """Print torsion data of molecule to screen."""
    print(fileio.GetPrintTorsionsString(self.torsions, self.atoms))

  def PrintOutofplanes(self):
    """Print outofplane data of molecule to screen."""
    print(fileio.GetPrintOutofplanesString(self.outofplanes, self.atoms))

  def PrintGradient(self):
    """Print gradient data to screen."""
    comment = mol.grad_type + ' total gradient'
    print(fileio.GetPrintGradientString(mol.atoms, self.g_total, comment))
