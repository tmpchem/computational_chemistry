"""Classes and functions for unit testing the mmlib energy module."""

import unittest

from mmlib import energy

class TestGetEBond(unittest.TestCase):

  def testNoBondSpring(self):
    params = r_ij, r_eq, k_b = 2.0, 1.0, 0.0
    self.assertAlmostEqual(energy.GetEBond(*params), 0.0)

  def testEquilibriumBond(self):
    params = r_ij, r_eq, k_b = 2.0, 2.0, 1.0
    self.assertAlmostEqual(energy.GetEBond(*params), 0.0)

  def testLongBond(self):
    params = r_ij, r_eq, k_b = 2.5, 2.0, 5.0
    self.assertAlmostEqual(energy.GetEBond(*params), 1.25)

  def testShortBond(self):
    params = r_ij, r_eq, k_b = 1.5, 2.0, 5.0
    self.assertAlmostEqual(energy.GetEBond(*params), 1.25)

class TestGetEAngle(unittest.TestCase):

  def testNoAngleSpring(self):
    params = a_ijk, a_eq, k_a = 180.0, 90.0, 0.0
    self.assertAlmostEqual(energy.GetEAngle(*params), 0.0)

  def testEquilibriumAngle(self):
    params = a_ijk, a_eq, k_a = 146.6, 146.6, 1.0
    self.assertAlmostEqual(energy.GetEAngle(*params), 0.0)

  def testLargeAngle(self):
    params = a_ijk, a_eq, k_a = 138.5, 93.5, 4.6
    self.assertAlmostEqual(energy.GetEAngle(*params), 2.8375113)

  def testSmallAngle(self):
    params = a_ijk, a_eq, k_a = 48.5, 93.5, 4.6
    self.assertAlmostEqual(energy.GetEAngle(*params), 2.8375113)

class TestGetETorsion(unittest.TestCase):

  def testNoTorsionBarrier(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 93.5, 0.0, 4.7, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.0)

  def testMinTorsionBarrier(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 180.0, 1.0, 0.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.0)

  def testMaxTorsionBarrier(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 0.0, 1.0, 0.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testPositiveGamma(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 15.0, 1.0, 15.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testNegativeGamma(self):
    params = t_ijkl, v_n, gamma, nfold, paths = -165.0, 1.0, -165.0, 1, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testThreefold(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 120.0, 1.0, 0.0, 3, 1
    self.assertAlmostEqual(energy.GetETorsion(*params), 2.0)

  def testFivePaths(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 0.0, 1.0, 0.0, 1, 5
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.4)

  def testAllTorsionParams(self):
    params = t_ijkl, v_n, gamma, nfold, paths = 136.7, 1.8, -14.4, 3, 9
    self.assertAlmostEqual(energy.GetETorsion(*params), 0.2861022)

class TestGetEOutofplane(unittest.TestCase):

  def testNoOutofplaneBarrier(self):
    params = o_ijkl, v_n = 45.0, 0.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.0)

  def testInPlaneAngle(self):
    params = o_ijkl, v_n = 0.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.0)

  def testMaxOutofplaneAngle(self):
    params = o_ijkl, v_n = 90.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 2.0)

  def testMinOutofplaneAngle(self):
    params = o_ijkl, v_n = -90.0, 1.0
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 2.0)

  def testAllOutofplaneParams(self):
    params = o_ijkl, v_n = 16.7, 1.3
    self.assertAlmostEqual(energy.GetEOutofplane(*params), 0.2146978)

class TestGetEVdwIJ(unittest.TestCase):

  def testNoEpsilon(self):
    params = r_ij, eps_ij, ro_ij = 2.5, 0.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), 0.0)

  def testLargeDistance(self):
    params = r_ij, eps_ij, ro_ij = 4.0, 3.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), -0.0930176)

  def testSmallDistance(self):
    params = r_ij, eps_ij, ro_ij = 1.6, 3.0, 2.0
    self.assertAlmostEqual(energy.GetEVdwIJ(*params), 20.7675621)

class TestGetEElstIJ(unittest.TestCase):

  def testNoCharge1(self):
    params = r_ij, q_i, q_j, epsilon = 1.0, 0.0, 1.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 0.0)

  def testNoCharge2(self):
    params = r_ij, q_i, q_j, epsilon = 1.0, 1.0, 0.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 0.0)

  def testUnitValues(self):
    params = r_ij, q_i, q_j, epsilon = 1.0, 1.0, 1.0, 1.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), 332.0637500)

  def testAllElstParams(self):
    params = r_ij, q_i, q_j, epsilon = 1.5, 0.2, -0.4, 2.0
    self.assertAlmostEqual(energy.GetEElstIJ(*params), -8.8550333)

def suite():
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

