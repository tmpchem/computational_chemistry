"""Main file for analyzing and plotting molecular dynamics simulation data.

This program reads in a set of molecular geometry and energy data output by a
molecular simulation and computes and plots resulting ensemble data.

No guarantees are made that the results of this program are correct and the
author assumes no liability for their reliability.
"""

from mmlib import analyze
from mmlib import fileio

__author__ = 'Trent M. Parker'
__email__ = 'tmpchemistry@gmail.com'
__status__ = 'Prototype'
__date__ = '2017-02-22'

if __name__ == '__main__':
  # Check input syntax.
  input_file_name = fileio.ValidateInput(__file__)

  # Read in ensemble geometry and energy analysis data.
  analysis = analyze.Analysis(input_file_name)

  # Compute and plot ensemble properties.
  analysis.RunAnalysis()
