import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mmlib

#############################################################################
#                             Welcome to opt.py                             #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and charges,     #
# determines bonded topology and parameters, calculates amber94 molecular   #
# mechanics energy components, gradients, and displaces the atomic          #
# coordinates to a potential energy local minimum.                          #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 02/16/2017                                #
#############################################################################

# check input syntax
infile_name = mmlib.fileio.GetInput()

# read in molecular and optimization data
opt = mmlib.optimize.Optimization(infile_name)

# run energy minimization
opt.Optimize()
