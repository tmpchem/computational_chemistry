"""Functions for testing proper performance of mmlib modules.

Warning: work in progress. Module is incomplete, and is subject to substantial
change without notice until it achieves full testing coverage and a stable form.
"""

import ast
import math
import numpy
import os
import sys

from mmlib import energy
from mmlib import fileio
from mmlib import geomcalc
from mmlib import topology

# absolute and relative thresholds for float equality comparisons
comp_coef = 1.0
comp_exp = -6
abstol = comp_coef * 10**(comp_exp) 
reltol = comp_coef * 10**(comp_exp) 

# boolean equality strings
bool_eq_char = {True: '==', False: '!='}

# boolean pass strings
bool_pass = {True: 'PASS', False: 'FAIL'}

# directory with test files
mm_dir = os.path.dirname(os.path.realpath(__file__))
test_dir = os.path.join(mm_dir, 'tests')

# names of functions to be unit tested
test_funcs = os.listdir(test_dir)
test_funcs = ''.join(test_funcs).replace('-', '.').split('.dat')[:-1]

# file with test case reference values
test_files = []
for i in range(len(test_funcs)):
  test_files.append(test_dir + os.sep + test_funcs[i].replace('.', '-')
      + '.dat')

# read in test case inputs and references from file
def read_in_tests(file_name):
  infile = open(file_name, "r")
  infile_lines = infile.readlines()
  infile.close()
  tests = []
  for line in infile_lines:
    if (not line[0] == '#' and len(line.split(';')) > 1):
      inputs = line.split(';')[0]
      ref = ast.literal_eval(line.split(';')[1].replace(' ',''))
      tests.append([inputs, ref])
  return tests

# equality comparison of generic data types
def equality_comparison(val, ref):
  same = True
  if   ('float' in str(type(val))):
    same *= math.isclose(ref, val, rel_tol=reltol, abs_tol=abstol)
  elif (type(val) == int or type(val) == str):
    same *= (val == ref)
  else:
    if (not len(val) == len(ref)):
      return False
    for i in range(len(val)):
      same *= equality_comparison(val[i], ref[i])
  return same

# print result of an individual unit test
def print_success_test(index, n_tests, success, val, ref, printval):
  if (printval >= 2):
    n_dig = int(math.ceil(math.log10(n_tests)))
    print('Test %0*i -> ' % (n_dig, index), end='')
    print('%s(%i)' % (bool_pass[success], success), end='')
  if (printval >= 3):
    eq = bool_eq_char[success]
    if (type(val) in [float, numpy.float64]):
      ref_dig = int(math.floor(math.log10(max(ref, abstol))))
      print_dig = max(0, -comp_exp - ref_dig*(ref_dig>0))
      print(', %10.*f %s %10.*f' % (print_dig, val, eq, print_dig, ref), end='')
    elif (type(val) in [int, str]):
      print(', %10s %s %10s' % (val, eq, ref), end='')
  if (printval >= 2):
    print('')

# print result of an individual function test
def print_success_function(index, n_pass, n_fail, name, printval):
  n_tests = n_pass + n_fail
  success = (n_pass == n_tests)
  if (printval >= 1):
    print('--- Function Test %i -> ' % (index), end='')
    print('%s(%i) ---' % (bool_pass[success], success))
    print('(%s), %i/%i Subtests Passed\n' % (name, n_pass, n_tests))

# run unit tests for a given function
def run_tests(printval):
  print('')
  nf_tests = len(test_files)
  nf_pass, nf_fail = 0, 0
  for j in range(len(test_files)):
    tests = read_in_tests(test_files[j])
    name = test_funcs[j]
    n_tests = len(tests)
    n_pass, n_fail = 0, 0
    for i in range(n_tests):
      inputs = tests[i][0]
      ref = tests[i][1]
      evalstr = '%s(%s)' % (test_funcs[j], inputs)
      val = eval(evalstr)
      success = equality_comparison(val, ref)
      print_success_test(i+1, n_tests, success, val, ref, printval)
      n_pass += int(success)
      n_fail += int(not success)
    print_success_function(j+1, n_pass, n_fail, name, printval)
    nf_pass += int(n_pass == n_tests)
    nf_fail += int(not n_pass == n_tests)
  all_success = (nf_pass == nf_tests)
  if (printval >= 0):
    print('------------------------')
    print('| All Tests -> ', end='')
    print('%s(%i) |' % (bool_pass[all_success], all_success))
    print('------------------------')
    print('%i/%i Function Tests Passed\n' % (nf_pass, nf_tests))
  return all_success
