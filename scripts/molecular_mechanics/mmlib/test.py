from mmlib import geomcalc, fileio, topology, energy
import os, sys, ast, math

# absolute and relative thresholds for float equality comparisons
comp_coef = 1.0
comp_exp = -7
abstol = comp_coef * 10**(comp_exp) 
reltol = comp_coef * 10**(comp_exp) 

# boolean equality strings
bool_eq_char = {True: '==', False: '!='}

# boolean pass strings
bool_pass = {True: 'PASS', False: 'FAIL'}

# directory with test files
mm_dir = '/'.join(os.path.abspath(energy.__file__).split('/')[:-1])
test_dir = mm_dir + '/test_files/'

# names of functions to be unit tested
test_funcs = ['geomcalc.get_r2_ij']

# file with test case reference values
test_files = []
for i in range(len(test_funcs)):
    test_files.append(test_dir + test_funcs[i].replace('.', '-') + '.dat')

# read in test case inputs and references from file
def read_in_tests(file_name):
    infile = open(file_name, "r")
    infile_lines = infile.readlines()
    infile.close()
    tests = []
    for line in infile_lines:
        if (not line[0] == '#'):
            inputs = line.split(';')[0]
            ref = ast.literal_eval(line.split(';')[1].replace(' ',''))
            tests.append([inputs, ref])
    return tests

# print result of an individual unit test
def print_success_test(index, n_tests, success, val, ref, printval):
    if (printval >= 2):
        n_dig = int(math.ceil(math.log10(n_tests)))
        print('Test %0*i -> ' % (n_dig, index), end='')
        print('%s(%i)' % (bool_pass[success], success), end='')
    if (printval >= 3):
        ref_dig = int(math.floor(math.log10(max(ref, abstol))))
        print_dig = max(0, -comp_exp - ref_dig*(ref_dig>0))
        eq = bool_eq_char[success]
        print(' %10.*f %s %-10.*f' % (print_dig, val, eq, print_dig, ref), end='')
    if (printval >= 2):
        print('')

# print result of an individual function test
def print_success_function(index, n_pass, n_fail, name, printval):
    n_tests = n_pass + n_fail
    success = (n_pass == n_tests)
    if (printval >= 1):
        print('  %i/%i Subtests Passed,' % (n_pass, n_tests), end='')
        print(' Func Test %i\n  (%s) ' % (index, name), end='')
        print('-> %s(%i)' % (bool_pass[success], int(success)))

# run unit tests for a given function
def run_tests(printval):
    for j in range(len(test_files)):
        tests = read_in_tests(test_files[j])
        name = test_funcs[j]
        n_tests = len(tests)
        n_pass = 0
        n_fail = 0
        for i in range(n_tests):
            inputs = tests[i][0]
            ref = tests[i][1]
            evalstr = '%s(%s)' % (test_funcs[j], inputs)
            val = eval(evalstr)
            success = math.isclose(ref, val, rel_tol=reltol, abs_tol=abstol)
            print_success_test(i+1, n_tests, success, val, ref, printval)
            n_pass += int(success)
            n_fail += int(not success)
        print_success_function(j+1, n_pass, n_fail, name, printval)

