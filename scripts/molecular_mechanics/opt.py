"""Main file for executing molecular mechanics energy minimization.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
components and gradients, and performs an energy minimization algorithm to
achieve a minimum energy coordinate configuration.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Initial version as of 2017-02-16'

import mmlib

if __name__ == '__main__':
  # check input syntax
  infile_name = mmlib.fileio.GetInput()

  # read in molecular and optimization data
  opt = mmlib.optimize.Optimization(infile_name)

  # run energy minimization
  opt.Optimize()
