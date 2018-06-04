"""Classes and functions for unit testing the mmlib geomcalc module."""

import numpy
import unittest

from mmlib import geomcalc
from mmlib import test

# Input 3d Cartesian points.

ORIGIN = numpy.array([0.0, 0.0, 0.0])

POSITIVE_UNIT_X = numpy.array([1.0, 0.0, 0.0])
POSITIVE_UNIT_Y = numpy.array([0.0, 1.0, 0.0])
POSITIVE_UNIT_Z = numpy.array([0.0, 0.0, 1.0])
NEGATIVE_UNIT_X = numpy.array([-1.0,  0.0,  0.0])
NEGATIVE_UNIT_Y = numpy.array([ 0.0, -1.0,  0.0])
NEGATIVE_UNIT_Z = numpy.array([ 0.0,  0.0, -1.0])

POSITIVE_ARBITRARY_X = numpy.array([ 5.606651, 0.0, 0.0])
NEGATIVE_ARBITRARY_X = numpy.array([-6.703290, 0.0, 0.0])
POSITIVE_ARBITRARY_Y = numpy.array([0.0,  7.711671, 0.0])
NEGATIVE_ARBITRARY_Y = numpy.array([0.0, -4.529674, 0.0])
POSITIVE_ARBITRARY_Z = numpy.array([0.0, 0.0,  2.916089])
NEGATIVE_ARBITRARY_Z = numpy.array([0.0, 0.0, -9.708677])

ARBITRARY_XY1 = numpy.array([ 2.836500,  5.236433,  0.000000])
ARBITRARY_XY2 = numpy.array([ 8.717305, -6.188432,  0.000000])
ARBITRARY_XY3 = numpy.array([-8.733996,  5.445714,  0.000000])
ARBITRARY_XY4 = numpy.array([-5.084499, -1.011108,  0.000000])
ARBITRARY_XZ1 = numpy.array([ 2.975653,  0.000000,  5.139741])
ARBITRARY_XZ2 = numpy.array([ 7.726639,  0.000000, -2.505934])
ARBITRARY_XZ3 = numpy.array([-7.996581,  0.000000,  2.001329])
ARBITRARY_XZ4 = numpy.array([-2.788969,  0.000000, -7.492930])
ARBITRARY_YZ1 = numpy.array([ 0.000000,  2.268236,  8.391077])
ARBITRARY_YZ2 = numpy.array([ 0.000000,  1.520357, -3.379237])
ARBITRARY_YZ3 = numpy.array([ 0.000000, -0.601497,  2.640605])
ARBITRARY_YZ4 = numpy.array([ 0.000000, -5.658143, -0.942184])

ARBITRARY_UNIT_XY1 = numpy.array([ 0.476296,  0.879285,  0.000000])
ARBITRARY_UNIT_XY2 = numpy.array([ 0.815421, -0.578869,  0.000000])
ARBITRARY_UNIT_XY3 = numpy.array([-0.848567,  0.529088,  0.000000])
ARBITRARY_UNIT_XY4 = numpy.array([-0.980795, -0.195042,  0.000000])
ARBITRARY_UNIT_XZ1 = numpy.array([ 0.501038,  0.000000,  0.865425])
ARBITRARY_UNIT_XZ2 = numpy.array([ 0.951223,  0.000000, -0.308504])
ARBITRARY_UNIT_XZ3 = numpy.array([-0.970080,  0.000000,  0.242785])
ARBITRARY_UNIT_XZ4 = numpy.array([-0.348833,  0.000000, -0.937185])
ARBITRARY_UNIT_YZ1 = numpy.array([ 0.000000,  0.260949,  0.965352])
ARBITRARY_UNIT_YZ2 = numpy.array([ 0.000000,  0.410297, -0.911958])
ARBITRARY_UNIT_YZ3 = numpy.array([ 0.000000, -0.222098,  0.975024])
ARBITRARY_UNIT_YZ4 = numpy.array([ 0.000000, -0.986418, -0.164257])

ARBITRARY_XYZ1 = numpy.array([ 3.688426,  7.267425,  5.136008])
ARBITRARY_XYZ2 = numpy.array([ 2.679605,  2.861745, -6.162862])
ARBITRARY_XYZ3 = numpy.array([ 4.999002, -9.340230,  0.395498])
ARBITRARY_XYZ4 = numpy.array([ 9.780405, -9.433397, -0.080764])
ARBITRARY_XYZ5 = numpy.array([-2.967049,  8.213730,  3.731393])
ARBITRARY_XYZ6 = numpy.array([-9.652579,  6.223187, -8.962634])
ARBITRARY_XYZ7 = numpy.array([-3.162726, -1.052111,  6.240231])
ARBITRARY_XYZ8 = numpy.array([-4.986630, -0.831000, -9.742808])

ARBITRARY_UNIT_XYZ1 = numpy.array([ 0.382886,  0.754414,  0.533157])
ARBITRARY_UNIT_XYZ2 = numpy.array([ 0.366860,  0.391797, -0.843747])
ARBITRARY_UNIT_XYZ3 = numpy.array([ 0.471549, -0.881051,  0.037307])
ARBITRARY_UNIT_XYZ4 = numpy.array([ 0.719747, -0.694211, -0.005943])
ARBITRARY_UNIT_XYZ5 = numpy.array([-0.312421,  0.864881,  0.392904])
ARBITRARY_UNIT_XYZ6 = numpy.array([-0.662584,  0.427179, -0.615224])
ARBITRARY_UNIT_XYZ7 = numpy.array([-0.447052, -0.148716,  0.882059])
ARBITRARY_UNIT_XYZ8 = numpy.array([-0.454308, -0.075709, -0.887622])

# Output 3d Cartesian vectors.

UNIT_VECTOR_12 = numpy.array([-0.08289876, -0.36203192, -0.92847223])
UNIT_VECTOR_21 = numpy.array([ 0.08289876,  0.36203192,  0.92847223])

CROSS_UNIT_VECTOR_12 = numpy.array([-0.84550446,  0.51870233, -0.12676282])
CROSS_UNIT_VECTOR_21 = numpy.array([ 0.84550446, -0.51870233,  0.12676282])

CROSS_VECTOR_12 = numpy.array([-59.48608258,  36.49373315, -8.9184937])
CROSS_VECTOR_21 = numpy.array([ 59.48608258, -36.49373315,  8.9184937])

CROSS_VECTOR_XY12 = numpy.array([  0.00000000,  0.00000000, -63.20107094])
CROSS_VECTOR_XZ12 = numpy.array([ -0.00000000, 47.16971329,   0.00000000])
CROSS_VECTOR_YZ12 = numpy.array([-20.42233967,  0.00000000,   0.00000000])

class TestGetR2ij(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetR2ij method."""

  def testSamePoint(self):
    """Asserts zero distance between identical points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ1
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 0.0)

  def testUnitLength(self):
    """Asserts unit value for points one distance unit apart."""
    params = ORIGIN, NEGATIVE_UNIT_X
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 1.0)

  def testOnAxisX(self):
    """Asserts correct value for points separated along x-axis."""
    params = POSITIVE_ARBITRARY_X, NEGATIVE_ARBITRARY_X
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 151.5346474)

  def testOnAxisY(self):
    """Asserts correct value for points separated along y-axis."""
    params = NEGATIVE_ARBITRARY_Y, POSITIVE_ARBITRARY_Y
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 149.8505274)

  def testOnAxisZ(self):
    """Asserts correct value for points separated along z-axis."""
    params = POSITIVE_ARBITRARY_Z, NEGATIVE_ARBITRARY_Z
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 159.3847166)

  def testArbitrary(self):
    """Asserts correct value for arbitrary points in 3d space."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 148.0921993)

  def testReflexive(self):
    """Asserts same value for inverted order of inputs."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ1
    self.assertAlmostEqual(geomcalc.GetR2ij(*params), 148.0921993)


class TestGetRij(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetRij method."""

  def testSamePoint(self):
    """Asserts zero distance between identical points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ1
    self.assertAlmostEqual(geomcalc.GetRij(*params), 0.0)

  def testUnitLength(self):
    """Asserts unit length for points one distance unit apart."""
    params = ORIGIN, NEGATIVE_UNIT_X
    self.assertAlmostEqual(geomcalc.GetRij(*params), 1.0)

  def testOnAxisX(self):
    """Asserts correct value for points separated along x-axis."""
    params = POSITIVE_ARBITRARY_X, NEGATIVE_ARBITRARY_X
    self.assertAlmostEqual(geomcalc.GetRij(*params), 12.3099410)

  def testOnAxisY(self):
    """Asserts correct value for points separated along y-axis."""
    params = NEGATIVE_ARBITRARY_Y, POSITIVE_ARBITRARY_Y
    self.assertAlmostEqual(geomcalc.GetRij(*params), 12.2413450)

  def testOnAxisZ(self):
    """Asserts correct value for points separated along z-axis."""
    params = POSITIVE_ARBITRARY_Z, NEGATIVE_ARBITRARY_Z
    self.assertAlmostEqual(geomcalc.GetRij(*params), 12.6247660)

  def testArbitrary(self):
    """Asserts correct distance between arbitrary points in 3d space."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2
    self.assertAlmostEqual(geomcalc.GetRij(*params), 12.1693138)

  def testReflexive(self):
    """Asserts same value for inverted order of inputs."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ1
    self.assertAlmostEqual(geomcalc.GetRij(*params), 12.1693138)


class TestGetUij(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetUij method."""

  def assertListAlmostEqual(self, test_vector, reference_vector):
    """Aliases vector float near equality function from mmlib.test module."""
    test.assertListAlmostEqual(self, test_vector, reference_vector)

  def testSamePoint(self):
    """Asserts zero vector between identical points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ1
    self.assertListAlmostEqual(geomcalc.GetUij(*params), ORIGIN)

  def testUnitLength(self):
    """Asserts correct direction for points separated by unit distance."""
    params = POSITIVE_UNIT_X, ORIGIN
    self.assertListAlmostEqual(geomcalc.GetUij(*params), NEGATIVE_UNIT_X)

  def testOnAxisX(self):
    """Asserts correct on-axis unit vector between points on the x-axis."""
    params = POSITIVE_ARBITRARY_X, NEGATIVE_ARBITRARY_X
    self.assertListAlmostEqual(geomcalc.GetUij(*params), NEGATIVE_UNIT_X)

  def testOnAxisY(self):
    """Asserts correct on-axis unit vector between points on the y-axis."""
    params = NEGATIVE_ARBITRARY_Y, POSITIVE_ARBITRARY_Y
    self.assertListAlmostEqual(geomcalc.GetUij(*params), POSITIVE_UNIT_Y)

  def testOnAxisZ(self):
    """Asserts correct on-axis unit vector between points on the z-axis."""
    params = POSITIVE_ARBITRARY_Z, NEGATIVE_ARBITRARY_Z
    self.assertListAlmostEqual(geomcalc.GetUij(*params), NEGATIVE_UNIT_Z)

  def testArbitrary(self):
    """Asserts correct unit vector between arbitrary points in 3d space."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2
    self.assertListAlmostEqual(geomcalc.GetUij(*params), UNIT_VECTOR_12)

  def testAntiReflexive(self):
    """Asserts opposite value for inverted order of inputs."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ1
    self.assertListAlmostEqual(geomcalc.GetUij(*params), UNIT_VECTOR_21)


class TestGetUdp(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetUdp method."""

  def testZeroVectorOne(self):
    """Asserts zero output with zero vector first input."""
    params = ORIGIN, ARBITRARY_UNIT_XYZ1
    self.assertAlmostEqual(geomcalc.GetUdp(*params), 0.0)

  def testZeroVectorTwo(self):
    """Asserts zero output with zero vector second input."""
    params = ARBITRARY_UNIT_XYZ2, ORIGIN
    self.assertAlmostEqual(geomcalc.GetUdp(*params), 0.0)

  def testUnitX(self):
    """Asserts correct value with one vector along positive x-axis."""
    params = POSITIVE_UNIT_X, ARBITRARY_UNIT_XYZ1
    self.assertAlmostEqual(geomcalc.GetUdp(*params), 0.382886)

  def testUnitY(self):
    """Asserts correct value with one vector along negative y-axis."""
    params = ARBITRARY_UNIT_XYZ2, NEGATIVE_UNIT_Y
    self.assertAlmostEqual(geomcalc.GetUdp(*params), -0.391797)

  def testUnitZ(self):
    """Asserts correct value with one vector along negative z-axis."""
    params = NEGATIVE_UNIT_Z, ARBITRARY_UNIT_XYZ2
    self.assertAlmostEqual(geomcalc.GetUdp(*params), 0.843747)

  def testSameVector(self):
    """Asserts unit value with same vector in each input."""
    params = ARBITRARY_UNIT_XYZ3, ARBITRARY_UNIT_XYZ3
    self.assertAlmostEqual(geomcalc.GetUdp(*params), 1.0)

  def testOppositeVector(self):
    """Asserts negative unit value with opposite vectors in inputs."""
    params = POSITIVE_UNIT_X, NEGATIVE_UNIT_X
    self.assertAlmostEqual(geomcalc.GetUdp(*params), -1.0)

  def testArbitrary(self):
    """Asserts correct value with two arbitrary unit vectors."""
    params = ARBITRARY_UNIT_XYZ1, ARBITRARY_UNIT_XYZ2
    self.assertAlmostEqual(geomcalc.GetUdp(*params), -0.0138069)

  def testReflexive(self):
    """Asserts same value with inverted order of inputs."""
    params = ARBITRARY_UNIT_XYZ2, ARBITRARY_UNIT_XYZ1
    self.assertAlmostEqual(geomcalc.GetUdp(*params), -0.0138069)


class TestGetUcp(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetUcp method."""

  def assertListAlmostEqual(self, test_vector, reference_vector):
    """Aliases vector float near equality function from mmlib.test module."""
    test.assertListAlmostEqual(self, test_vector, reference_vector)

  def testZeroVector(self):
    """Asserts zero vector result for zero vector input."""
    params = ORIGIN, ARBITRARY_UNIT_XYZ1
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), ORIGIN)

  def testUnitXY(self):
    """Asserts positive z-axis as result of x-axis and y-axis cross."""
    params = POSITIVE_UNIT_X, POSITIVE_UNIT_Y
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), POSITIVE_UNIT_Z)

  def testUnitYX(self):
    """Asserts negative z-axis as result of y-axis and x-axis cross."""
    params = POSITIVE_UNIT_Y, POSITIVE_UNIT_X
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), NEGATIVE_UNIT_Z)

  def testUnitYZ(self):
    """Asserts positive x-axis as result of y-axis and z-axis cross."""
    params = POSITIVE_UNIT_Y, POSITIVE_UNIT_Z
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), POSITIVE_UNIT_X)

  def testUnitZX(self):
    """Asserts positive y-axis as result of z-axis and x-axis cross."""
    params = POSITIVE_UNIT_Z, POSITIVE_UNIT_X
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), POSITIVE_UNIT_Y)

  def testInPlaneXY(self):
    """Asserts z-axis as result of two xy-plane vectors."""
    params = ARBITRARY_UNIT_XY1, ARBITRARY_UNIT_XY2
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), NEGATIVE_UNIT_Z)

  def testInPlaneXZ(self):
    """Asserts y-axis as result of two xz-plane vectors."""
    params = ARBITRARY_UNIT_XZ1, ARBITRARY_UNIT_XZ2
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), POSITIVE_UNIT_Y)

  def testInPlaneYZ(self):
    """Asserts x-axis as result of two xz-plane vectors."""
    params = ARBITRARY_UNIT_YZ3, ARBITRARY_UNIT_YZ4
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), POSITIVE_UNIT_X)

  def testArbitrary(self):
    """Asserts correct value with two arbitrary unit vectors."""
    params = ARBITRARY_UNIT_XYZ1, ARBITRARY_UNIT_XYZ2
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), CROSS_UNIT_VECTOR_12)

  def testAntiReflexive(self):
    """Asserts opposite value with inverted input order."""
    params = ARBITRARY_UNIT_XYZ2, ARBITRARY_UNIT_XYZ1
    self.assertListAlmostEqual(geomcalc.GetUcp(*params), CROSS_UNIT_VECTOR_21)


class TestGetCp(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetCp method."""

  def assertListAlmostEqual(self, test_vector, reference_vector):
    """Aliases vector float near equality function from mmlib.test module."""
    test.assertListAlmostEqual(self, test_vector, reference_vector)

  def testZeroVector(self):
    """Asserts zero vector result for zero vector input."""
    params = ORIGIN, ARBITRARY_XYZ1
    self.assertListAlmostEqual(geomcalc.GetCp(*params), ORIGIN)

  def testUnitXY(self):
    """Asserts positive z-axis as result of x-axis and y-axis cross."""
    params = POSITIVE_UNIT_X, POSITIVE_UNIT_Y
    self.assertListAlmostEqual(geomcalc.GetCp(*params), POSITIVE_UNIT_Z)

  def testUnitYX(self):
    """Asserts negative z-axis as result of y-axis and x-axis cross."""
    params = POSITIVE_UNIT_Y, POSITIVE_UNIT_X
    self.assertListAlmostEqual(geomcalc.GetCp(*params), NEGATIVE_UNIT_Z)

  def testUnitYZ(self):
    """Asserts positive x-axis as result of y-axis and z-axis cross."""
    params = POSITIVE_UNIT_Y, POSITIVE_UNIT_Z
    self.assertListAlmostEqual(geomcalc.GetCp(*params), POSITIVE_UNIT_X)

  def testUnitZX(self):
    """Asserts positive y-axis as result of z-axis and x-axis cross."""
    params = POSITIVE_UNIT_Z, POSITIVE_UNIT_X
    self.assertListAlmostEqual(geomcalc.GetCp(*params), POSITIVE_UNIT_Y)

  def testInPlaneXY(self):
    """Asserts z-axis as result of two xy-plane vectors."""
    params = ARBITRARY_XY1, ARBITRARY_XY2
    self.assertListAlmostEqual(geomcalc.GetCp(*params), CROSS_VECTOR_XY12)

  def testInPlaneXZ(self):
    """Asserts y-axis as result of two xz-plane vectors."""
    params = ARBITRARY_XZ1, ARBITRARY_XZ2
    self.assertListAlmostEqual(geomcalc.GetCp(*params), CROSS_VECTOR_XZ12)

  def testInPlaneYZ(self):
    """Asserts x-axis as result of two xz-plane vectors."""
    params = ARBITRARY_YZ1, ARBITRARY_YZ2
    self.assertListAlmostEqual(geomcalc.GetCp(*params), CROSS_VECTOR_YZ12)

  def testArbitrary(self):
    """Asserts correct value with two arbitrary unit vectors."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2
    self.assertListAlmostEqual(geomcalc.GetCp(*params), CROSS_VECTOR_12)

  def testAntiReflexive(self):
    """Asserts opposite value with inverted input order."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ1
    self.assertListAlmostEqual(geomcalc.GetCp(*params), CROSS_VECTOR_21)


class TestGetAijk(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetAijk method."""

  def testZeroAngle(self):
    """Asserts zero angle between identical vectors."""
    params = ORIGIN, NEGATIVE_UNIT_Z, ORIGIN
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 0.0)

  def testRightAngle(self):
    """Asserts ninety degree angle between orthogonal vectors."""
    params = POSITIVE_UNIT_X, ORIGIN, POSITIVE_UNIT_Y
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 90.0)

  def testStraightAngle(self):
    """Asserts 180 degree angle between oppositie vectors."""
    params = POSITIVE_UNIT_Z, ORIGIN, NEGATIVE_UNIT_Z
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 180.0)

  def testInPlaneXY(self):
    """Asserts correct value within the xy-plane."""
    params = ARBITRARY_XY1, ARBITRARY_XY2, ARBITRARY_XY3
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 29.07348364)

  def testInPlaneXZ(self):
    """Asserts correct value within the xz-plane."""
    params = ARBITRARY_XZ1, ARBITRARY_XZ2, ARBITRARY_XZ3
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 42.14774907)

  def testInPlaneYZ(self):
    """Asserts correct value within the yz-plane."""
    params = ARBITRARY_YZ1, ARBITRARY_YZ2, ARBITRARY_YZ3
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 23.05201986)

  def testArbitraryAcute(self):
    """Asserts correct value between arbitrary acute 3d points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2, ARBITRARY_XYZ3
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 82.37365716)

  def testArbitraryObtuse(self):
    """Asserts correct value between arbitrary obtuse 3d points."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ3, ARBITRARY_XYZ4
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 97.75040394)

  def testReflexive(self):
    """Asserts same value when inverting the order of points."""
    params = ARBITRARY_XYZ4, ARBITRARY_XYZ3, ARBITRARY_XYZ2
    self.assertAlmostEqual(geomcalc.GetAijk(*params), 97.75040394)


class TestGetTijkl(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetTijkl method."""

  def testZeroTorsion(self):
    """Asserts zero torsion within the xy-plane."""
    params = POSITIVE_UNIT_X, ORIGIN, POSITIVE_UNIT_Y, ARBITRARY_UNIT_XY1
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), 0.0)

  def testPositiveRightTorsion(self):
    """Asserts positive ninety degree torsion combination."""
    params = POSITIVE_UNIT_X, ORIGIN, POSITIVE_UNIT_Y, NEGATIVE_UNIT_Z
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), 90.0)

  def testNegativeRightTorsion(self):
    """Asserts negative ninety degree torsion combination."""
    params = POSITIVE_UNIT_X, ORIGIN, POSITIVE_UNIT_Y, POSITIVE_UNIT_Z
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), -90.0)

  def testPositiveRightTorsion(self):
    """Asserts positive ninety degree torsion combination."""
    params = POSITIVE_UNIT_X, ORIGIN, POSITIVE_UNIT_Y, ARBITRARY_UNIT_XY4
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), 180.0)

  def testArbitraryNegative(self):
    """Asserts correct negative value between arbitrary 3d points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2, ARBITRARY_XYZ3, ARBITRARY_XYZ4
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), -92.01024227)

  def testArbitraryPositive(self):
    """Asserts correct negative value between arbitrary 3d points."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ3, ARBITRARY_XYZ4, ARBITRARY_XYZ5
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), 37.31802340)

  def testReflixe(self):
    """Asserts same value when inverting the order of points."""
    params = ARBITRARY_XYZ5, ARBITRARY_XYZ4, ARBITRARY_XYZ3, ARBITRARY_XYZ2
    self.assertAlmostEqual(geomcalc.GetTijkl(*params), 37.31802340)


class TestGetOijkl(unittest.TestCase):
  """Unit tests for mmlib.geomcalc.GetOijkl method."""

  def testInPlane(self):
    """Asserts zero angle for in-plane points."""
    params = ARBITRARY_XY1, ARBITRARY_XY2, ARBITRARY_XY3, ARBITRARY_XY4
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), 0.0)

  def testMaxAngle(self):
    """Asserts ninety degrees as maximum out-of-plane angle."""
    params = ARBITRARY_XY2, ARBITRARY_XY1, ORIGIN, POSITIVE_ARBITRARY_Z
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), 90.0, places=5)

  def testMinAngle(self):
    """Asserts negative ninety degrees as minimum out-of-plane angle."""
    params = ARBITRARY_XY2, ARBITRARY_XY1, ORIGIN, NEGATIVE_ARBITRARY_Z
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), -90.0, places=5)

  def testArbitraryNegative(self):
    """Asserts correct negative value between arbitrary 3d points."""
    params = ARBITRARY_XYZ1, ARBITRARY_XYZ2, ARBITRARY_XYZ3, ARBITRARY_XYZ4
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), -81.99467957)

  def testArbitraryPositive(self):
    """Asserts correct negative value between arbitrary 3d points."""
    params = ARBITRARY_XYZ2, ARBITRARY_XYZ3, ARBITRARY_XYZ4, ARBITRARY_XYZ5
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), 28.81956723)

  def testAntiReflexive(self):
    """Asserts correct negative value between arbitrary 3d points."""
    params = ARBITRARY_XYZ3, ARBITRARY_XYZ2, ARBITRARY_XYZ4, ARBITRARY_XYZ5
    self.assertAlmostEqual(geomcalc.GetOijkl(*params), -28.81956723)


def suite():
  """Builds a test suite of all unit tests in geomcalc_test module."""
  test_classes = (
      TestGetR2ij,
      TestGetRij,
      TestGetUij,
      TestGetUdp,
      TestGetUcp,
      TestGetCp,
      TestGetAijk,
      TestGetTijkl,
      TestGetOijkl)
  
  suite = unittest.TestSuite()
  for test_class in test_classes:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
  return suite

