import os
import sys

sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mmlib

# amount of printing
printval = 3

# run all tests
mmlib.test.RunTests(printval)
