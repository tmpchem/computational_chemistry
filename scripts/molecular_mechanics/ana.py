"""Main file for analyzing and plotting molecular dynamics simulation data.

This program reads in a set of molecular geometry and energy data output by a
molecular simulation and computes and plots resulting ensemble data.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

import mmlib

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2017-02-22'

if __name__ == '__main__':
  # check input syntax
  infile_name = mmlib.fileio.GetInput()

  # read in ensemble geometry and energy data
  ana = mmlib.analyze.Analysis(infile_name)

  # compute and plot ensemble properties
  ana.RunAnalysis()
