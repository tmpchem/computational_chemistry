"""Classes and functions for unit testing the mmlib fileio module."""

import unittest

from mmlib import fileio
from mmlib import molecule
from mmlib import test

class TestGetAtomFromXyzq(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAtomsFromXyzq method."""

  def testEmptyArray(self):
    """Asserts error raise for empty string input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAtomFromXyzq([])
    self.assertEqual(
        str(e.exception),
        'Insufficient columns to parse xyzq file row into Atom object: ')

  def testShortArray(self):
    """Asserts error raise for insufficient string size."""
    param = ['H', '0.0', '0.0', '0.0']
    with self.assertRaises(IndexError) as e:
      fileio.GetAtomFromXyzq(param)
    self.assertEqual(
        str(e.exception),
        ('Insufficient columns to parse xyzq file row into Atom object: '
         'H 0.0 0.0 0.0'))

  def testBadAtomType(self):
    """Asserts error raise for non-alphanumeric atom type input."""
    param = ['_', '0.0', '0.0', '0.0', '0.0']
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromXyzq(param)
    self.assertEqual(str(e.exception), 'Atom type must begin with letter: _')

  def testBadCoordinate(self):
    """Asserts error raise for non-numeric coordinate input."""
    param = ['C', '0.0', '0.0', 'A', '0.0']
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromXyzq(param)
    self.assertEqual(str(e.exception),
                     'Atomic coordinates must be numeric values: 0.0 0.0 A')

  def testBadCharge(self):
    """Asserts error raise for non-numeric charge input."""
    param = ['Ne', '-1.0', '0.0', '+1.0', '?']
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromXyzq(param)
    self.assertEqual(str(e.exception), 'Atomic charge must be numeric value: ?')

  def testZeroValues(self):
    """Asserts correct output for zero value numeric inputs."""
    param = ['OW', '0.0', '0.0', '0.0', '0.0']
    output = molecule.Atom('OW', [0.0, 0.0, 0.0], 0.0)
    test.assertObjectEqual(self, fileio.GetAtomFromXyzq(param), output)

  def testArbitrary(self):
    """Asserts correct output for arbitrary Atom inputs."""
    param = ['N*', '999.9', '-0.0001', '2.0', '-0.52']
    output = molecule.Atom('N*', [999.9, -0.0001, 2.0], -0.52)
    test.assertObjectEqual(self, fileio.GetAtomFromXyzq(param), output)


class TestGetAtomsFromXyzq(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAtomsFromXyzq method."""
  
  def testEmptyArray(self):
    """Asserts error raise for empty input."""
    with self.assertRaises(EOFError) as e:
      fileio.GetAtomsFromXyzq([])
    self.assertEqual(str(e.exception), 'XYZQ file is empty.')

  def testEmptyFirstRow(self):
    """Asserts error raise for empty first row of input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAtomsFromXyzq([[], ['X', '0.0', '0.0', '0.0', '0.0']])
    self.assertEqual(str(e.exception),
                     'First line of XYZQ file must be number of atoms.')

  def testZeroAtoms(self):
    """Asserts error raise for non-positive number of atoms."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomsFromXyzq([['0']])
    self.assertEqual(str(e.exception), 'Number of atoms must be positive: 0')\

  def testTooFewRows(self):
    """Asserts error raise for too few rows for atom number."""
    with self.assertRaises(EOFError) as e:
      fileio.GetAtomsFromXyzq([['1'], ['comment', 'row']])
    self.assertEqual(
        str(e.exception), 
        'XYZQ file does not contain enough lines for stated number of atoms: 1')

  def testOneAtom(self):
    """Asserts correct atom returned from one atom input."""
    param = [['1'], ['ignored', 'comment'], ['H', '0.0', '0.0', '0.0', '0.0']]
    test_output = fileio.GetAtomsFromXyzq(param)
    ref_output = [molecule.Atom('H', [0.0, 0.0, 0.0], 0.0)]
    test.assertObjectArrayEqual(self, test_output, ref_output)

  def testTwoAtoms(self):
    """Asserts correct atoms returned from two atom input."""
    param = [['2'], [], ['N*', '1.0', '0.1', '-.2', '+0.553'],
             ['H4', '0.523', '-100.12', '5.5', '-0.27']]
    test_output = fileio.GetAtomsFromXyzq(param)
    ref_output = [molecule.Atom('N*', [1.0, 0.1, -0.2], 0.553),
                  molecule.Atom('H4', [0.523, -100.12, 5.5], -0.27)]
    test.assertObjectArrayEqual(self, test_output, ref_output)


class TestGetAtomFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAtomFromPrm method."""

  def testEmptyArray(self):
    """Asserts error raise for empty input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAtomFromPrm([])
    self.assertEqual(
        str(e.exception),
        'Insufficient columns to parse prm file row into Atom object: ')

  def testShortArray(self):
    """Asserts error raise for not enough fields in row."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAtomFromPrm(['ATOM', '1', 'CT', '0', '0', '0', '0', '1'])
    self.assertEqual(
        str(e.exception),
        ('Insufficient columns to parse prm file row into Atom object: '
         'ATOM 1 CT 0 0 0 0 1'))

  def testBadAtomicType(self):
    """Asserts error raise for improper atom type input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromPrm(
          ['ATOM', '1', '/N', '0.0', '0.0', '0.0', '0.0', '1.0', '0.0'])
    self.assertEqual(str(e.exception), 'Atom type must begin with letter: /N')

  def testBadCoordinates(self):
    """Asserts error raise for improper atomic coordinate input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromPrm(
          ['ATOM', '2', 'O', '0.0', '&', '0.0', '0.0', '0.1', '0.0'])
    self.assertEqual(str(e.exception),
                     'Atomic coordinates must be numeric values: 0.0 & 0.0')

  def testBadCharge(self):
    """Asserts error raise for improper atomic charge input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromPrm(
          ['ATOM', '3', 'CT', '0.1', '0.1', '0.1', '0.0.0', '10', '10'])
    self.assertEqual(str(e.exception),
                     'Atomic charge must be numeric value: 0.0.0')

  def testBadRadius(self):
    """Asserts error raise for improper atomic vdw radius input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromPrm(
          ['ATOM', '4', 'HC', '-0.1', '-0.0', '2.0', '.1', '(1.0)', '0.0'])
    self.assertEqual(str(e.exception),
                     'Atomic vdw radius must be positive numeric value: (1.0)')

  def testBadEpsilon(self):
    """Asserts error raise for improper atomic epsilon input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAtomFromPrm(
          ['ATOM', '5', 'CC', '1.', '2.', '3.', '4.', '5.', '.'])
    self.assertEqual(str(e.exception),
                     'Atomic epsilon must be non-negative numeric value: .')

  def testZeroValues(self):
    """Asserts correct Atom object for zero value numeric inputs."""
    param = ['ATOM', '6', 'X', '0.', '0.', '0.', '0.', '0.1', '0.']
    test.assertObjectEqual(self, fileio.GetAtomFromPrm(param),
                           molecule.Atom('X', [0., 0., 0.], 0., 0.1, 0.))

  def testArbitrary(self):
    """Asserts correct Atom object for arbitrary input parameters."""
    param = ['ATOM', '7', 'NA', '1.1', '2.2', '3.3', '4.4', '5.5', '6.6']
    test.assertObjectEqual(self, fileio.GetAtomFromPrm(param),
                           molecule.Atom('NA', [1.1, 2.2, 3.3], 4.4, 5.5, 6.6))


class TestGetBondFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetBondFromPrm method."""

  def testEmptyArray(self):
    """Asserts error raise for empty array is input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetBondFromPrm([])
    self.assertEqual(
        str(e.exception),
        'Insufficient columns to parse prm file row into Bond object: ')

  def testShortArray(self):
    """Asserts error raise for not enough fields in input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetBondFromPrm(['BOND', '1', '2', '1.0'])
    self.assertEqual(
        str(e.exception),
        ('Insufficient columns to parse prm file row into Bond object: '
         'BOND 1 2 1.0'))

  def testBadIndex1(self):
    """Asserts error raise for non-positive-integer atomic index 1 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetBondFromPrm(['BOND', '0', '1', '1.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Bond atom1 index must be positive integer: 0')

  def testBadIndex2(self):
    """Asserts error raise for non-positive-integer atomic index 2 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetBondFromPrm(['BOND', '2', 'A', '1.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Bond atom2 index must be positive integer: A')

  def testBadBondSpringConstant(self):
    """Asserts error raise for non-positive numeric bond spring constant."""
    with self.assertRaises(ValueError) as e:
      fileio.GetBondFromPrm(['BOND', '100', '101', '0.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Bond spring constant must be positive numeric value: 0.0')

  def testBadEquilibriumBondLength(self):
    """Asserts error raise for negative or non-numeric equilibrium length."""
    with self.assertRaises(ValueError) as e:
      fileio.GetBondFromPrm(['BOND', '10000', '3', '0.1', '#'])
    self.assertEqual(
        str(e.exception),
        'Equilibrium bond length must be positive numeric value: #')

  def testSmallestValues(self):
    """Asserts correct Bond object for smallest allowed values."""
    param = ['BOND', '1', '2', '0.000001', '0.000001']
    test.assertObjectEqual(self, fileio.GetBondFromPrm(param),
                           molecule.Bond(0, 1, 0.000001, 0.000001))

  def testArbitrary(self):
    """Asserts correct Bond object for arbitrary input parameters."""
    param = ['BOND', '455', '27', '8.484', '6.63']
    test.assertObjectEqual(self, fileio.GetBondFromPrm(param),
                           molecule.Bond(454, 26, 8.484, 6.63))


class TestGetAngleFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAngleFromPrm method."""

  def testEmptyArray(self):
    """Asserts error raise for empty array is input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAngleFromPrm([])
    self.assertEqual(
        str(e.exception),
        'Insufficient columns to parse prm file row into Angle object: ')

  def testShortArray(self):
    """Asserts error raise for not enough fields in input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '1', '2', '3', '1.0'])
    self.assertEqual(
        str(e.exception),
        ('Insufficient columns to parse prm file row into Angle object: '
         'ANGLE 1 2 3 1.0'))

  def testBadIndex1(self):
    """Asserts error raise for non-positive-integer atomic index 1 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '0', '1', '2', '1.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Angle atom1 index must be positive integer: 0')

  def testBadIndex2(self):
    """Asserts error raise for non-positive-integer atomic index 2 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '2', 'A', '1', '1.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Angle atom2 index must be positive integer: A')

  def testBadIndex3(self):
    """Asserts error raise for non-positive-integer atomic index 3 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '2', '1', '1^', '1.0', '1.0'])
    self.assertEqual(str(e.exception),
                     'Angle atom3 index must be positive integer: 1^')

  def testBadAngleSpringConstant(self):
    """Asserts error raise for non-positive numeric angle spring constant."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '100', '101', '102', '-1.0', '1.0'])
    self.assertEqual(
        str(e.exception),
        'Angle spring constant must be positive numeric value: -1.0')

  def testBadEquilibriumBondAngle(self):
    """Asserts error raise for negative or non-numeric equilibrium angle."""
    with self.assertRaises(ValueError) as e:
      fileio.GetAngleFromPrm(['ANGLE', '9999', '2', '34', '0.1', '180.1'])
    self.assertEqual(
        str(e.exception),
        'Equilibrium bond angle must be between 0.0 and 180.0: 180.1')

  def testSmallestValues(self):
    """Asserts correct Angle object for smallest allowed values."""
    param = ['ANGLE', '1', '2', '3', '0.000001', '0.000001']
    test.assertObjectEqual(self, fileio.GetAngleFromPrm(param),
                           molecule.Angle(0, 1, 2, 0.000001, 0.000001))

  def testArbitrary(self):
    """Asserts correct Angle object for arbitrary input parameters."""
    param = ['ANGLE', '62', '64', '67', '179.838', '17.76']
    test.assertObjectEqual(self, fileio.GetAngleFromPrm(param),
                           molecule.Angle(61, 63, 66, 179.838, 17.76))


class TestGetTorsionFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetTorsionFromPrm method."""

  def testEmptyArray(self):
    """Asserts error raise for empty array is input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetTorsionFromPrm([])
    self.assertEqual(
        str(e.exception),
        'Insufficient columns to parse prm file row into Torsion object: ')

  def testShortArray(self):
    """Asserts error raise for not enough fields in input."""
    with self.assertRaises(IndexError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '1', '2', '3', '4', '1.0', '0.0', '1'])
    self.assertEqual(
        str(e.exception),
        ('Insufficient columns to parse prm file row into Torsion object: '
         'TORSION 1 2 3 4 1.0 0.0 1'))

  def testBadIndex1(self):
    """Asserts error raise for non-positive-integer atomic index 1 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '0', '1', '2', '3', '1.0', '0.0', '1', '1'])
    self.assertEqual(str(e.exception),
                     'Torsion atom1 index must be positive integer: 0')

  def testBadIndex2(self):
    """Asserts error raise for non-positive-integer atomic index 2 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '2', 'A', '1', '4', '1.0', '0.0', '1', '1'])
    self.assertEqual(str(e.exception),
                     'Torsion atom2 index must be positive integer: A')

  def testBadIndex3(self):
    """Asserts error raise for non-positive-integer atomic index 3 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '2', '1', '1^', '3', '1.0', '0.0', '1', '1'])
    self.assertEqual(str(e.exception),
                     'Torsion atom3 index must be positive integer: 1^')

  def testBadIndex4(self):
    """Asserts error raise for non-positive-integer atomic index 4 input."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '10', '11', '12', 'i', '1.0', '0.0', '1', '1'])
    self.assertEqual(str(e.exception),
                     'Torsion atom4 index must be positive integer: i')

  def testBadTorsionBarrier(self):
    """Asserts error raise for non-positive numeric torsion half-barrier."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '100', '101', '102', '103', '-1.0', '1.0', '2', '6'])
    self.assertEqual(
        str(e.exception),
        'Torsion half-barrier must be positive numeric value: -1.0')

  def testBadTorsionOffsetAngle(self):
    """Asserts error raise for non-numeric or out-of-bounds offset angle."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '9999', '2', '34', '6', '0.1', '180.1', '1', '3'])
    self.assertEqual(
        str(e.exception),
        'Torsion offset angle must be between -180.0 and 180.0: 180.1')

  def testBadTorsionFrequency(self):
    """Asserts error raise for non-postive-integer torsion frequency."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '1', '2', '3', '4', '1.0', '0.0', '0', '3'])
    self.assertEqual(str(e.exception),
                     'Torsion barrier frequency must be positive integer: 0')

  def testBadTorsionPaths(self):
    """Asserts error raise for non-positive-integer torsion path number."""
    with self.assertRaises(ValueError) as e:
      fileio.GetTorsionFromPrm(
          ['TORSION', '5', '4', '3', '1', '2.5', '30.0', '3', '-2'])
    self.assertEqual(str(e.exception),
                     'Torsion path number must be positive integer: -2')

  def testSmallestValues(self):
    """Asserts correct Torsion object for smallest allowed values."""
    param = ['TORSION', '1', '2', '3', '4', '0.000001', '-180.0', '1', '1']
    test.assertObjectEqual(self, fileio.GetTorsionFromPrm(param),
                           molecule.Torsion(0, 1, 2, 3, 0.000001, -180.0, 1, 1))

  def testArbitrary(self):
    """Asserts correct Torsion object for arbitrary input parameters."""
    param = ['ANGLE', '62', '64', '67', '51', '10.0', '90.0', '3', '9']
    test.assertObjectEqual(self, fileio.GetTorsionFromPrm(param),
                           molecule.Torsion(61, 63, 66, 50, 10.0, 90.0, 3, 9))


class TestGetOutofplaneFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetOutofplaneFromPrm method."""
  pass


class TestGetAtomsFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAtomsFromPrm method."""
  pass


class TestGetBondsFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetBondsFromPrm method."""
  pass


class TestGetAnglesFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetAnglesFromPrm method."""
  pass


class TestGetTorsionsFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetTorsionsFromPrm method."""
  pass


class TestGetOutofplanesFromPrm(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetOutofplanesFromPrm method."""
  pass


class TestGetPrintCoordsXyzString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintCoordsXyzString method."""
  pass


class TestGetPrintGradientString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintGradientString method."""
  pass


class TestGetPrintGeomString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintGeomString method."""
  pass


class TestGetPrintBondsString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintBondsString method."""
  pass


class TestGetPrintAnglesString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintAnglesString method."""
  pass


class TestGetPrintTorsionsString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintTorsionsString method."""
  pass


class TestGetPrintOutofplanesString(unittest.TestCase):
  """Unit tests for mmlib.fileio.GetPrintOutofplanesString method."""
  pass


class TestValidateInput(unittest.TestCase):
  """Unit tests for mmlib.fileio.TestValidateInput method."""
  pass


def suite():
  """Builds a test suite of all unit tests in fileio_test module."""
  test_classes = (
      TestGetAtomFromXyzq,
      TestGetAtomsFromXyzq,
      TestGetAtomFromPrm,
      TestGetBondFromPrm,
      TestGetAngleFromPrm,
      TestGetTorsionFromPrm,
      TestGetOutofplaneFromPrm,
      TestGetBondsFromPrm,
      TestGetAnglesFromPrm,
      TestGetTorsionsFromPrm,
      TestGetOutofplanesFromPrm,
      TestGetPrintCoordsXyzString,
      TestGetPrintGradientString,
      TestGetPrintGeomString,
      TestGetPrintBondsString,
      TestGetPrintAnglesString,
      TestGetPrintTorsionsString,
      TestGetPrintOutofplanesString,
      TestValidateInput)
  
  suite = unittest.TestSuite()
  for test_class in test_classes:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
  return suite
