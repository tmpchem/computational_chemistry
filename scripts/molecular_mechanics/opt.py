"""Main file for executing molecular mechanics energy minimization.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
components and gradients, and performs an energy minimization algorithm to
achieve a minimum energy coordinate configuration.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import fileio
from mmlib import optimize

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2017-02-16'

if __name__ == '__main__':
  # Check input syntax
  input_file_name = fileio.ValidateInput(__file__)

  # Read in molecular and optimization data
  optimization = optimize.Optimization(input_file_name)

  # Run energy minimization.
  optimization.Optimize()
