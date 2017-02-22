import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import mmlib

#############################################################################
#                           Welcome to analysis.py                          #
#                                                                           #
# This program read in a set of molecular geometry and energy data          #
# previously output from a molecular mechanics simulation and computes      #
# and plots resulting dynamic and ensemble data.                            #
#                                                                           #
# No guarantees are made that the results of this program are correct       #
# and the author assumes no liability for their reliability.                #
#                                                                           #
#                              Trent M. Parker                              #
#                                 02/22/2017                                #
#############################################################################

# check input syntax
infile_name = mmlib.fileio.get_input()

# read in ensemble geometry and energy data
ana = mmlib.analyze.Analysis(infile_name)

# compute and plot ensemble properties
ana.run_analysis()

# end of program

