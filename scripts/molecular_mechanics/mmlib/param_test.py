"""Classes and functions for unit testing the mmlib geomcalc module."""

import unittest

from mmlib import param

# Long parameter array for outofplane tests.
_CT_CT_N_C_TORSION_PARAMS = [
    (0.00,   0.0,  2, 6),
    (0.50, 180.0, -4, 1),
    (0.15, 180.0, -3, 1),
    (0.00, 180.0, -2, 1),
    (0.53,   0.0,  1, 1)]

class TestGetMass(unittest.TestCase):
  """Unit tests for mmlib.param.GetMass method."""

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetMass('')

    self.assertEqual(str(e.exception), 'No atomic mass found for element: ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing string input."""
    with self.assertRaises(ValueError) as e:
      param.GetMass('AA')

    self.assertEqual(str(e.exception), 'No atomic mass found for element: AA')

  def testDummyInput(self):
    """Asserts zero mass for dummy atom input."""
    self.assertEqual(param.GetMass('X'), 0.0)

  def testHydrogen(self):
    """Asserts correct mass for hydrogen atom."""
    self.assertEqual(param.GetMass('H'), 1.00794)

  def testHeaviestAtom(self):
    """Asserts correct mass for heaviest present atom."""
    self.assertEqual(param.GetMass('Kr'), 83.7980)


class TestGetCovRad(unittest.TestCase):
  """Unit tests for mmlib.param.GetCovRad method."""

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetCovRad('')
    
    self.assertEqual(str(e.exception),
                     'No covalent radius found for element: ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing element input."""
    with self.assertRaises(ValueError) as e:
      param.GetCovRad('AA')

    self.assertEqual(str(e.exception),
                         'No covalent radius found for element: AA')

  def testDummyInput(self):
    """Asserts zero radius for dummy atom input."""
    self.assertEqual(param.GetCovRad('X'), 0.0)

  def testHydrogen(self):
    """Asserts correct radius for hydrogen atom."""
    self.assertEqual(param.GetCovRad('H'), 0.37)

  def testHeaviestAtom(self):
    """Asserts correct radius for heaviest present atom."""
    self.assertEqual(param.GetCovRad('Kr'), 1.03)


class TestGetVdwParam(unittest.TestCase):
  """Unit tests for mmlib.param.GetVdwParam method."""

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetVdwParam('')
    
    self.assertEqual(str(e.exception), 'No vdw param found for atom type: ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing atom type input."""
    with self.assertRaises(ValueError) as e:
      param.GetVdwParam('AA')

    self.assertEqual(str(e.exception), 'No vdw param found for atom type: AA')

  def testDummyInput(self):
    """Asserts zero parameters for dummy atom type."""
    self.assertEqual(param.GetVdwParam('X'), (0.0000, 0.00000))

  def testAromaticCarbon(self):
    """Asserts correct parameters for carbon atom type."""
    self.assertEqual(param.GetVdwParam('CA'), (1.9080, 0.0860))

  def testWaterOxygen(self):
    """Asserts correct parameters for water oxygen atom type."""
    self.assertEqual(param.GetVdwParam('OW'), (1.7683, 0.1520))


class TestGetBondParam(unittest.TestCase):
  """Unit tests for mmlib.param.GetBondParam method."""

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetBondParam('', '')
    
    self.assertEqual(str(e.exception),
                     'No bond parameters found for atom type pair: , ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing atom types input."""
    with self.assertRaises(ValueError) as e:
      param.GetBondParam('AA', 'BB')

    self.assertEqual(str(e.exception),
                     'No bond parameters found for atom type pair: AA, BB')

  def testDummyInput(self):
    """Asserts zero parameters for dummy atom types input."""
    self.assertEqual(param.GetBondParam('X', 'X'), (0.0, 0.0))

  def testOxygenHydrogenBond(self):
    """Asserts correct parameters for water hydrogen-oxygen bond."""
    self.assertEqual(param.GetBondParam('OW', 'HW'), (553.0, 0.9572))

  def testCarbonHydrogenReversedBond(self):
    """Asserts correct parameters for aromatic reversed carbon-hydrogen bond."""
    self.assertEqual(param.GetBondParam('HA', 'CA'), (367.0, 1.080))


class TestGetAngleParam(unittest.TestCase):
  """Unit tests for mmlib.param.GetAngleParam method."""
  pass

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetAngleParam('', '', '')
    
    self.assertEqual(str(e.exception),
                     'No angle parameters found for atom type triplet: , , ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing string input."""
    with self.assertRaises(ValueError) as e:
      param.GetAngleParam('AA', 'BB', 'CC')

    self.assertEqual(
        str(e.exception),
        'No angle parameters found for atom type triplet: AA, BB, CC')

  def testDummyInput(self):
    """Asserts zero mass for dummy atom input."""
    self.assertEqual(param.GetAngleParam('X', 'X', 'X'), (0.0, 0.0))

  def testHOHWater(self):
    """Asserts correct mass for hydrogen atom."""
    self.assertEqual(param.GetAngleParam('HW', 'OW', 'HW'), (100.0, 104.52))

  def testCCHBenzene(self):
    """Asserts correct mass for heaviest present atom."""
    self.assertEqual(param.GetAngleParam('HA', 'CA', 'CA'), (35.0, 120.0))


class TestGetTorsionParam(unittest.TestCase):
  """Unit tests for mmlib.param.GetTorsionParam method."""

  def testEmptyInput(self):
    """Asserts raise of ValueError for empty string input."""
    with self.assertRaises(ValueError) as e:
      param.GetTorsionParam('', '', '', '')

    self.assertEqual(str(e.exception),
                     'No torsion parameters found for central atom pair: , ')

  def testBadInput(self):
    """Asserts raise of ValueError for missing string input."""
    with self.assertRaises(ValueError) as e:
      param.GetTorsionParam('AA', 'BB', 'CC', 'DD')

    self.assertEqual(
        str(e.exception),
        'No torsion parameters found for central atom pair: BB, CC')

  def testCentral(self):
    """Asserts correct parameters for given central atom pair."""
    self.assertEqual(param.GetTorsionParam('X', 'C', 'CA', 'X'),
                     [(14.50, 180.0, 2, 4)])

  def testCentralReversed(self):
    """Asserts correct parameters for given reversed central atom pair."""
    self.assertEqual(param.GetTorsionParam('X', 'CA', 'C', 'X'),
                     [(14.50, 180.0, 2, 4)])

  def testFourAtom(self):
    """Asserts correct parameters for given atom quartet."""
    self.assertEqual(param.GetTorsionParam('CT', 'CT', 'N', 'C'),
                     _CT_CT_N_C_TORSION_PARAMS)

  def testFourAtomReversed(self):
    """Asserts correct parameters for given reversed atom quartet."""
    self.assertEqual(param.GetTorsionParam('C', 'N', 'CT', 'CT'),
                     _CT_CT_N_C_TORSION_PARAMS)


class TestGetOutofplaneParam(unittest.TestCase):
  """Unit tests for mmlib.param.GetOutofplaneParam method."""

  def testEmptyInput(self):
    """Asserts zero value for empty string input."""
    self.assertEqual(param.GetOutofplaneParam('', '', '', ''), 0.0)

  def test34(self):
    """Asserts correct parameter for ending atom pair."""
    self.assertEqual(param.GetOutofplaneParam('', '', 'C', 'O'), 10.5)

  def test234(self):
    """Asserts correct parameter for ending atom triplet."""
    self.assertEqual(param.GetOutofplaneParam('', 'CT', 'N', 'CT'), 1.0)

  def test1234(self):
    """Asserts correct paramter for complete atom quartet."""
    self.assertEqual(param.GetOutofplaneParam('CA', 'CA', 'C', 'OH'), 1.1)


def suite():
  """Builds a test suite of all unit tests in param_test module."""
  test_classes = (
      TestGetVdwParam,
      TestGetMass,
      TestGetCovRad,
      TestGetBondParam,
      TestGetAngleParam,
      TestGetTorsionParam,
      TestGetOutofplaneParam)
  
  suite = unittest.TestSuite()
  for test_class in test_classes:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
  return suite

