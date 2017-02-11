import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import mmlib

# amount of printing
printval = 3

# run all tests
mmlib.test.run_tests(printval)

# end of program

