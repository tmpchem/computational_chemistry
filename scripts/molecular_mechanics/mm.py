"""Main file for computing AMBER94 molecular mechanics energies.

This program takes in a set of molecular coordinates and parameters, infers
bonded topology (if unspecified), computes AMBER94 molecular mechanics energy
components, and outputs the resulting values to screen.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import fileio
from mmlib import molecule

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2016-02-15'

if __name__ == '__main__':
  # Check input syntax.
  input_file_name = fileio.ValidateInput(__file__)

  # Read in molecular geometry and topology.
  molecule = molecule.Molecule(input_file_name)

  # Calculate potential energy.
  molecule.GetEnergy('nokinetic')

  # Print results to screen.
  molecule.PrintData()
