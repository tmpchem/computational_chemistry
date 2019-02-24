"""Classes and functions for unit tetsing the mmlib molecule module."""

import unittest

class Test(unittest.TestCase):
  """Unit tests for mmlib. method."""


def suite():
  """Builds a test suite of all unit test in _test module."""
  test_classes = (
    Test)

  suite = unittest.TestSuite()
  for test_class in test_classes:
    tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
    suite.addTests(tests)
  return suite
