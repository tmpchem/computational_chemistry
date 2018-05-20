"""Main file for computing AMBER94 molecular mechanics energies.

This program takes in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), computes AMBER94 molecular mechanics energy
components, and outputs the resulting values to screen.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Initial version as of 2016-02-15'

import mmlib

if __name__ == '__main__':
  # check input syntax
  infile_name = mmlib.fileio.GetInput()

  # read in molecular geometry and topology
  mol = mmlib.molecule.Molecule(infile_name)

  # calculate potential energy
  mol.GetEnergy('nokinetic')

  # print results to screen
  mol.PrintData()
