"""Main file for executing Monte Carlo molecular simulations.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
comopnents, and performs Metropolis Monte Carlo configuration propogation for a
specified duration.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import fileio
from mmlib import simulate

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2017-02-08'

if __name__ == '__main__':
  # Check input syntax.
  input_file_name = fileio.ValidateInput(__file__)

  # Read in molecular and simulation data.
  simulation = simulate.MonteCarlo(input_file_name)

  # Run molecular dynamics simulation.
  simulation.Run()
