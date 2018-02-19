import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mmlib

#############################################################################
#                              Welcome to mm.py                             #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and charges,     #
# determines bonded topology and parameters, calculates amber94 molecular   #
# mechanics energy components, and prints resulting values.                 #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 02/15/2016                                #
#############################################################################

# check input syntax
infile_name = mmlib.fileio.GetInput()

# read in molecular geometry and topology
mol = mmlib.molecule.Molecule(infile_name)

# calculate energy
mol.GetEnergy('nokinetic')

# print results to screen
mol.PrintData()
