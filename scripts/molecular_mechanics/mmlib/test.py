"""Classes for testing proper function of all modules within mmlib.

WARNING: Incomplete. Work in progress. Subject to change at any time.
"""

import numpy
import unittest

from mmlib import energy_test
from mmlib import fileio_test
from mmlib import geomcalc_test
from mmlib import param_test

# Message to print at conclusion of test suite.
_SUCCESS_MESSAGE = '\nAll tests succeeded. mmlib is ready for use.'
_FAILURE_MESSAGE = '\nSome tests failed. Results may not be reliable.'

def RunTests():
  """Executes test suites for all mmlib modules."""
  test_suites = [
      fileio_test.suite(),
      geomcalc_test.suite(),
      energy_test.suite(),
      param_test.suite()]

  combo_suite = unittest.TestSuite(test_suites)
  result = unittest.TextTestRunner().run(combo_suite)

  if (result.wasSuccessful()):
    print(_SUCCESS_MESSAGE)
  else:
    print(_FAILURE_MESSAGE)


def assertListAlmostEqual(test_case, test_vector, reference_vector):
  """Supplemental function for near equality comparison of float vectors.

  Args:
    test_case (unittest.TestCase): Unit test class instance.
    test_vector (type*): Array of numeric floating-point values to test.
    reference_vector (type*): Reference to assert equality of test_vector.
  """
  test_case.assertEqual(len(test_vector), len(reference_vector))

  for test_value, reference_value in zip(test_vector, reference_vector):
    if hasattr(test_value, '__len__'):
      assertListAlmostEqual(test_case, test_value, reference_value)
    else:
      test_case.assertAlmostEqual(test_value, reference_value, places=6)

def assertObjectEqual(test_case, test_object, reference_object):
  """Supplemental function for equality of all object attributes.

  Args:
    test_case (unittest.TestCase): Unit test class instance.
    test_object (type): Object of type with attributes to test.
    reference_object (type): Reference to assert equality of test_object.
  """
  test_case.assertEqual(test_object.__dict__.keys(),
                        reference_object.__dict__.keys())

  for attribute in test_object.__dict__:
    test_value = getattr(test_object, attribute)
    reference_value = getattr(reference_object, attribute)
    
    if isinstance(test_value, numpy.ndarray):
      assertListAlmostEqual(test_case, test_value, reference_value)
    else:
      test_case.assertAlmostEqual(test_value, reference_value)

def assertObjectArrayEqual(test_case, test_array, reference_array):
  """Supplemental function for equality of all objects in an array.

  Args:
    test_case (unittest.TestCase): Unit test class instance.
    test_array (type*): Array of objects of type with attributes to test.
    reference_array (type*): Reference to assert equality of test_array.
  """
  test_case.assertEqual(len(test_array), len(reference_array))

  for test_object, reference_object in zip(test_array, reference_array):
    assertObjectEqual(test_case, test_object, reference_object)
