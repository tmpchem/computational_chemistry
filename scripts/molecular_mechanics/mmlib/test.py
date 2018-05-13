"""Classes for testing proper function of all modules within mmlib.

WARNING: Incomplete. Work in progress. Subject to change at any time.
"""

import unittest

from mmlib import energy_test

def RunTests():
  test_suites = [
      energy_test.suite()]

  combo_suite = unittest.TestSuite(test_suites)
  unittest.TextTestRunner().run(combo_suite)

