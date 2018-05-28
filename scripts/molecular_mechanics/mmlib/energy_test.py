"""Classes and functions for unit testing the mmlib energy module."""

import unittest

from mmlib import energy

class TestGetEBond(unittest.TestCase):
  """Unit tests for mmlib.energy.GetEBond method."""

  def testNoBondSpring(self):
    """Asserts no bond energy for zero spring constant."""
    params = r_ij, r_eq, k_b = 2.0, 1.0, 0.0
    self.assertAlmostEqual(energy.GetEBond(*params), 0.0)

  def testEquilibriumBond(self):
    """Asserts no bond energy for equilibrium length bond."""
    params = r_ij, r_eq, k_b = 2.0, 2.0, 1.0
    self.assertAlmostEqual(energy.GetEBond(*params), 0.0)

  def testLongBond(self):
    """Asserts exact value for greater than equilibrium bond length."""
    params = r_ij, r_eq, k_b = 2.5, 2.0, 5.0
    self.assertAlmostEqual(energy.GetEBond(*params), 1.25)

  def testShortBond(self):
    """Asserts exact value for smaller than equilibrium bond length."""
    params = r_ij, r_eq, k_b = 1.5, 2.0, 5.0
    self.assertAlmostEqual(energy.GetEBond(*params), 1.25)

  def testInfiniteBond(self):
    """Asserts infinite energy for infinite bond length."""
    params = r_ij, r_eq, k_b = float('inf'), 2.0, 1.0
    self.assertAlmostEqual(energy.GetEBond(*params), float('inf'))


class TestGetEAngle(unittest.TestCase):
  """Unit tests for mmlib.energy.GetEAngle method."""

  def testNoAngleSpring(self):
    """Asserts no angle energy for zero spring constant."""
    params = a_ijk, a_eq, k_a = 180.0, 90.0, 0.0
    self.assertAlmostEqual(energy.GetEAngle(*params), 0.0)

  def testEquilibriumAngle(self):
    """Asserts no angle energy for equilibrium bond angle."""
    params = a_ijk, a_eq, k_a = 146.6, 146.6, 1.0
    self.assertAlmostEqual(energy.GetEAngle(*params), 0.0)

  def testLargeAngle(self):
    """Asserts correct value for greater than equilibrium bond angle."""
    params = a_ijk, a_eq, k_a = 138.5, 93.5, 4.6
    self.assertAlmostEqual(energy.GetEAngle(*params), 2.8375113)

  def testSmallAngle(self):
    """Asserts correct value for smaller than equilibrium bond angle."""
    params = a_ijk, a_eq, k_a = 48.5, 93.5, 4.6
    self.assertAlmostEqual(energy.GetEAngle(*params), 2.8375113)


class TestGetETorsion(unittest.TestCase):
  """Unit tests for mmlib.energy.GetETorsion method."""

  def testNoTorsionBarrier(self):
    """Asserts no torsion energy for zero rotation barrier."""
    params = t_ijkl, v_n, gamma, nfold, paths = 93.5, 0.0, 4.7, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.0)

  def testMinTorsionBarrier(self):
    """Asserts no torsion energy at maximum angle."""
    params = t_ijkl, v_n, gamma, nfold, paths = 180.0, 1.0, 0.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.0)

  def testMaxTorsionBarrier(self):
    """Asserts maximum torsion energy at minimum angle."""
    params = t_ijkl, v_n, gamma, nfold, paths = 0.0, 1.0, 0.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testPositiveGamma(self):
    """Asserts maximum torsion energy at positive barrier offset angle."""
    params = t_ijkl, v_n, gamma, nfold, paths = 15.0, 1.0, 15.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testNegativeGamma(self):
    """Asserts maximum torsion energy at negative barrier offset angle."""
    params = t_ijkl, v_n, gamma, nfold, paths = -165.0, 1.0, -165.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testThreefold(self):
    """Asserts maximum torsion energy at maximum of threefold barrier angle."""
    params = t_ijkl, v_n, gamma, nfold, paths = 120.0, 1.0, 0.0, 3, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testFivePaths(self):
    """Asserts one fifth of regular torsion energy with five paths."""
    params = t_ijkl, v_n, gamma, nfold, paths = 0.0, 1.0, 0.0, 1, 5
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.4)

  def testAllTorsionParams(self):
    """Asserts correct value for an arbitrary combination of parameters."""
    params = t_ijkl, v_n, gamma, nfold, paths = 136.7, 1.8, -14.4, 3, 9
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.2861022)


class TestGetEOutofplane(unittest.TestCase):
  """Unit tests for mmlib.energy.GetEOutofplane method."""

  def testNoOutofplaneBarrier(self):
    """Asserts no outofplane energy for zero barrier."""
    params = o_ijkl, v_n = 45.0, 0.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.0)

  def testInPlaneAngle(self):
    """Asserts no outofplane energy for in-plane angle."""
    params = o_ijkl, v_n = 0.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.0)

  def testMaxOutofplaneAngle(self):
    """Asserts maximum outofplane energy for ninety degree angle."""
    params = o_ijkl, v_n = 90.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 2.0)

  def testMinOutofplaneAngle(self):
    """Asserts maximum outofplane energy for negative ninety degree angle."""
    params = o_ijkl, v_n = -90.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 2.0)

  def testAllOutofplaneParams(self):
    """Asserts correct value for an arbitrary combination of parameters."""
    params = o_ijkl, v_n = 16.7, 1.3
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.2146978)


class TestGetEVdwIJ(unittest.TestCase):
  """Unit tests for mmlib.energy.GetEVdwIJ method."""

  def testNoEpsilon(self):
    """Asserts no vdw energy for zero epsilon."""
    params = r_ij, eps_ij, ro_ij = 2.5, 0.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), 0.0)

  def testEquilibriumDistance(self):
    """Asserts negative epsilon energy for equilibrium separation."""
    params = r_ij, eps_ij, ro_ij = 2.5, 2.0, 2.5
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), -2.0)

  def testLargeDistance(self):
    """Asserts correct value for greater than equilibrium separation."""
    params = r_ij, eps_ij, ro_ij = 4.0, 3.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), -0.0930176)

  def testSmallDistance(self):
    """Asserts correct value for smaller than equilibrium separation."""
    params = r_ij, eps_ij, ro_ij = 1.6, 3.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), 20.7675621)

  def testInfiniteDistance(self):
    """Asserts zero energy for infinite separation."""
    params = r_ij, eps_ij, ro_ij = float('inf'), 1.0, 1.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), 0.0)


class TestGetEElstIJ(unittest.TestCase):
  """Unit tests for mmlib.energy.GetEElstIJ method."""

  def testNoCharge1(self):
    """Asserts no elst energy when first charge is zero."""
    params = r_ij, q_i, q_j, epsilon = 1.0, 0.0, 1.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 0.0)

  def testNoCharge2(self):
    """Asserts no elst energy when second charge is zero."""
    params = r_ij, q_i, q_j, epsilon = 1.0, 1.0, 0.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 0.0)

  def testUnitValues(self):
    """Asserts correct value for unit charges and separation."""
    params = r_ij, q_i, q_j, epsilon = 1.0, 1.0, 1.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 332.0637500)

  def testLargeDielectric(self):
    """Asserts correct value for non-unit dielectric constant."""
    params = r_ij, q_i, q_j, epsilon = 1.0, 1.0, 1.0, 10.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 33.2063750)

  def testInfiniteDistance(self):
    """Asserts zero elst energy when separation is infinite."""
    params = r_ij, q_i, q_j, epsilon = float('inf'), 1.0, 1.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 0.0)

  def testAllElstParams(self):
    """Asserts correct value for arbitrary combination of parameters."""
    params = r_ij, q_i, q_j, epsilon = 1.5, 0.2, -0.4, 2.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), -8.8550333)


def suite():
  """Builds a test suite of all unit tests in energy_test module."""
  test_classes = (
      TestGetEBond,
      TestGetEAngle,
      TestGetETorsion,
      TestGetEOutofplane,
      TestGetEVdwIJ,
      TestGetEElstIJ)
  
  suite = unittest.TestSuite()
  for test_class in test_classes:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
  return suite

