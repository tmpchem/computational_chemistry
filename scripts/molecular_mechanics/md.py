"""Main file for executing molecular dynamics simulations.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
components and gradients, and performs molecular dynamics configuration
propogation for a specified duration.

No guarantees are made the the results of this program are correct and the
author assumes no liability for their reliability.
"""

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Initial version as of 2016-05-31'

import mmlib

if __name__ == '__main__':
  # check input syntax
  infile_name = mmlib.fileio.GetInput()

  # read in molecular and simulation data
  sim = mmlib.simulate.MolecularDynamics(infile_name)

  # run molecular dynamics
  sim.Run()
