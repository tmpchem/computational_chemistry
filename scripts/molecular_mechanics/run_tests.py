"""Main file for executing the test suite of mmlib directory modules.

This program imports all test classes / functions from each unit test module
in the mmlib directory and executes them as a test suite to verify proper
functionality and the absence of regression erros in the working copy of this
directory.

WARNING: The test suite is incomplete and subject to change at any time as test
coverage is updated.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import test

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Development'
__date__ = '2018-05-20'

if __name__ == '__main__':
  # Run all unit tests.
  test.RunTests()
