"""Classes for testing proper function of all modules within mmlib.

WARNING: Incomplete. Work in progress. Subject to change at any time.
"""

import unittest

from mmlib import energy_test
from mmlib import geomcalc_test
from mmlib import param_test

# Message to print at conclusion of test suite.
_SUCCESS_MESSAGE = '\nAll tests succeeded. mmlib is ready for use.'
_FAILURE_MESSAGE = '\nSome tests failed. Results may not be reliable.'

def RunTests():
  """Executes test suites for all mmlib modules."""
  test_suites = [
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
  """Supplemental function for near equality comparison of float vectors."""
  test_case.assertEqual(len(test_vector), len(reference_vector))

  for test_value, reference_value in zip(test_vector, reference_vector):
    test_case.assertAlmostEqual(test_value, reference_value, places=6)
