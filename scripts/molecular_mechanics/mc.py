import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mmlib

#############################################################################
#                             Welcome to mc.py                              #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and charges,     #
# determines bonded topology and parameters, calculates amber94 molecular   #
# mechanics energy components and Metropolis Monte Carlo configurations.    #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 02/08/2017                                #
#############################################################################

# check input syntax
infile_name = mmlib.fileio.GetInput()

# read in molecular and simulation data
sim = mmlib.simulate.MonteCarlo(infile_name)

# run molecular dynamics
sim.Run()
