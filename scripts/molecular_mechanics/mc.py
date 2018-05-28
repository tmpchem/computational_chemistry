"""Main file for executing Monte Carlo molecular simulations.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
comopnents, and performs Metropolis Monte Carlo configuration propogation for a
specified duration.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

import mmlib

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2017-02-08'

if __name__ == '__main__':
  # check input syntax
  infile_name = mmlib.fileio.ValidateInput(__file__)

  # read in molecular and simulation data
  sim = mmlib.simulate.MonteCarlo(infile_name)

  # run molecular dynamics
  sim.Run()
