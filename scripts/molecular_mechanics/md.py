import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import mmlib

#############################################################################
#                              Welcome to md.py                             #
#                                                                           #
# This program takes in a set of molecular xyz coordinates and charges,     #
# determines bonded topology and parameters, calculates amber94 molecular   #
# mechanics energy components, gradients, and molecular dynamics            #
# trajectories                                                              #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 05/31/2016                                #
#############################################################################

# check input syntax
infile_name = mmlib.fileio.get_input()

# read in molecular and simulation data
sim = mmlib.simulate.MolecularDynamics(infile_name)

# run molecular dynamics
sim.run()

# end of program

