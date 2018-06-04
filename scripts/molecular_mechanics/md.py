"""Main file for executing molecular dynamics simulations.

This program reads in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), calculates AMBER94 molecular mechanics energy
components and gradients, and performs molecular dynamics configuration
propogation for a specified duration.

No guarantees are made the the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import fileio
from mmlib import simulate

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2016-05-31'

if __name__ == '__main__':
  # Check input syntax
  infile_name = fileio.ValidateInput(__file__)

  # Read in molecular and simulation data
  simulation = simulate.MolecularDynamics(infile_name)

  # Run molecular dynamics simulation.
  simulation.Run()
