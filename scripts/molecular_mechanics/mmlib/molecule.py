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
  def __init__(self, type_, coords, charge, ro, eps, mass):
    self.type_ = type_
    self.element = fileio.GetElement(type_)

    self.SetCoords(coords)
    self.SetCharge(charge)
    self.SetRo(ro)
    self.SetEps(eps)
    self.SetMass(mass)
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

  def SetMass(self, mass):
    """Set new (float) atomic mass [g/mol]."""
    self.mass = mass

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
    grad (float): Energy gradient magnitude of bond [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, r_ij, r_eq, k_b):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetRij(r_ij)
    self.SetReq(r_eq)
    self.SetKb(k_b)

    self.GetEnergy()
    self.GetGradient()

  def SetAt1(self, at1):
    """Set new (int) at1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) at2."""
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

  def GetGradient(self):
    """Calculate bond gradient (float) [kcal/(mol*A)]."""
    self.grad = gradient.GetGBond(self.r_ij, self.r_eq, self.k_b)


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
    grad (float): Energy gradient magnitude of angle [kcal/(mol*A)].
  """
  def __init__(self, at1, at2, at3, a_ijk, a_eq, k_a):
    self.SetAt1(at1)
    self.SetAt2(at2)
    self.SetAt3(at3)
    self.SetAijk(a_ijk)
    self.SetAeq(a_eq)
    self.SetKa(k_a)

    self.GetEnergy()
    self.GetGradient()

  def SetAt1(self, at1):
    """Set new (int) at1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) at2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) at3."""
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

  def GetGradient(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.gradient = gradient.GetGAngle(self.a_ijk, self.a_eq, self.k_a)


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

    self.GetEnergy()
    self.GetGradient()

  def SetAt1(self, at1):
    """Set new (int) at1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) at2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) at3."""
    self.at3 = at3

  def SetAt4(self, at4):
    """Set new (int) at4."""
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

  def GetGradient(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.gradient = gradient.GetGTorsion(self.t_ijkl, self.v_n, self.gam,
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

    self.GetEnergy()
    self.GetGradient()

  def SetAt1(self, at1):
    """Set new (int) at1."""
    self.at1 = at1

  def SetAt2(self, at2):
    """Set new (int) at2."""
    self.at2 = at2

  def SetAt3(self, at3):
    """Set new (int) at3."""
    self.at3 = at3

  def SetAt4(self, at4):
    """Set new (int) at4."""
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

  def GetGradient(self):
    """Get gradient (float) [kcal/(mol*A)]."""
    self.gradient = gradient.GetGOutofplane(self.o_ijkl, self.v_n)


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
    boundtype (str): Type of boundary shape, 'cube', 'sphere', or 'none'.
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
    self.name = os.path.splitext(os.path.basename(self.infile))[0]
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
    self.k_box = 250.0
    self.bound = 1.0E10
    self.boundtype = 'sphere'
    self.origin = numpy.zeros(const.NUMDIM)
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
    fileio.GetGeom(self)

  def ReadInPrm(self):
    """Read in prm data from .prm input file."""
    fileio.GetPrm(self)

  def GetTopology(self):
    """Determine bonded topology of molecules from coordinates."""
    topology.GetBondGraph(self)
    topology.GetBonds(self)
    topology.GetAngles(self)
    topology.GetTorsions(self)
    topology.GetOutofplanes(self)
    topology.GetNonints(self)

  def GetEnergy(self, kintype):
    """Calculate (float) energy [kcal/mol] and all energy components."""
    energy.GetEBonds(self)
    energy.GetEAngles(self)
    energy.GetETorsions(self)
    energy.GetEOutofplanes(self)
    energy.GetENonbonded(self)
    energy.GetEBound(self)
    energy.GetEKinetic(self, kintype)
    energy.GetETotals(self)

  def GetGradient(self, grad_type):
    """Calculate analytical or numerical gradient of energy.
    
    Args:
      grad_type (str): Type of gradient:
        'analytic': exact, based on analytic derivatives.
        'numerical': approximate, based on numerical derivatives.
    """
    self.grad_type = grad_type
    if (grad_type == 'analytic'):
      self.GetAnalyticGradient()
    elif (grad_type == 'numerical'):
      self.GetNumericalGradient()
    else:
      raise ValueError('Unexpected gradient type: %s\n'
                       'Use \'analytic\' or \'numerical\'.' % (grad_type))

  def GetAnalyticGradient(self):
    """Calculate analytic (float**) gradient [kcal/(mol*A)] of energy."""
    gradient.GetGBonds(self)
    gradient.GetGAngles(self)
    gradient.GetGTorsions(self)
    gradient.GetGOutofplanes(self)
    gradient.GetGNonbonded(self)
    gradient.GetGBound(self)
    gradient.GetGTotals(self)

  def GetNumericalGradient(self):
    """Calculate numerical (float**) gradient [kcal/(mol*A)] of energy."""
    gradient.GetGNumerical(self)
    gradient.GetGTotals(self)

  def UpdateInternals(self):
    """Update current values of internal degrees of freedom."""
    topology.UpdateBonds(self)
    topology.UpdateAngles(self)
    topology.UpdateTorsions(self)
    topology.UpdateOutofplanes(self)

  def GetTemperature(self):
    """Calculate instantaneous kinetic temperature [K] of system."""
    energy.GetTemperature(self)

  def GetPressure(self):
    """Calculate instantaneous kinetic pressure [Pa] of system."""
    gradient.GetPressure(self)

  def GetVolume(self):
    """Caclculate approximate volume [A^3] of system."""
    geomcalc.GetVolume(self)

  def PrintData(self):
    """Print energy / geometry / topology data of molecule to screen."""
    fileio.PrintEnergy(self)
    fileio.PrintGeom(self)
    fileio.PrintBonds(self)
    fileio.PrintAngles(self)
    fileio.PrintTorsions(self)
    fileio.PrintOutofplanes(self)

  def PrintEnergy(self):
    """Print energy and component data of molecule to screen."""
    fileio.PrintEnergy(self)

  def PrintGeom(self):
    """Print geometry data of molecule to screen."""
    fileio.PrintGeom(self)

  def PrintBonds(self):
    """Print bond data of molecule to screen."""
    fileio.PrintBonds(self)

  def PrintAngles(self):
    """Print angle data of molecule to screen."""
    fileio.PrintAngles(self)

  def PrintTorsions(self):
    """Print torsion data of molecule to screen."""
    fileio.PrintTorsions(self)

  def PrintOutofplanes(self):
    """Print outofplane data of molecule to screen."""
    fileio.PrintOutofplanes(self)

  def PrintGradient(self):
    """Print gradient data to screen."""
    fileio.PrintGradient(self.g_total, mol.grad_type + ' total gradient')
