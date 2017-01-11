import os, sys
sys.path.append('/'.join(os.path.realpath(__file__).split('/')[:-1]))
from mmlib import *

#############################################################################
#                              Welcome to md.py                             #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and charges,     #
# determines bonded topology and parameters, calculates amber94 molecular   #
# mechanics energy components, gradients, and molecular dynamics            #
# trajectories                                                              #
#                                                                           #
# AMBER FF94 from Cornell et. al, J. Am. Chem. Soc. 1995, 117, 5179-5197.   #
# Atom types in Table 1. Bonded and Non-bonded parameters in Table 14.      #
# Energy function in Equation 1.                                            #
#                                                                           #
# Download AmberTools15 from "http://ambermd.org/AmberTools15-get.html".    #
# After unzipping, parameters located in "amber14/dat/leap/parm/parm94.dat" #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 05/31/2016                                #
#############################################################################

## MAIN BLOCK ##

# check input syntax
infile_name = fileio.get_input()

# read in molecular and simulation data
sim = dynamics.simulation(infile_name)

# run molecular dynamics
sim.run_dynamics('leapfrog')

# end of program

