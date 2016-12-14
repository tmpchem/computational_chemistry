from mmlib import fileio, topology, energy

# threshold energy exponent for passing tests
test_thresh = -8

# directory structure for test files
head_dir = 'C:\Users\Trent\Documents\TMPChem\Coding'
test_dir = head_dir + '\tests'

# compare if two values are equal within a threshold
def are_equal(a, b, cutoff):
    are_equal = (abs(a - b) < cutoff)
    return are_equal

# run tests to check for proper program function
def run_tests():
    tests = ['Methane Dimer', 'Water Dimer', 'Benzene dimer']
    N_tests = len(tests)

    test_energy_file = test_dir + '/energy_tests.dat'
    test_energies = fileio.get_file_string_array(test_energy_file)
    for i in range(N_tests):
        test_energies[i] = float(test_energies[i][0])
  
    ALL_SUCCESS = 1

    for i in range(N_tests):
        test_prefix = test_dir + '/test%i' % (i + 1)
        dimer_file = test_prefix + '_dimer.xyzq'
        monoA_file = test_prefix + '_monoA.xyzq'
        monoB_file = test_prefix + '_monoB.xyzq'
    
        dimer_geom = fileio.get_geom(dimer_file)
        monoA_geom = fileio.get_geom(monoA_file)
        monoB_geom = fileio.get_geom(monoB_file)
    
        dimer_topology = topology.get_topology(dimer_geom)
        monoA_topology = topology.get_topology(monoA_geom)
        monoB_topology = topology.get_topology(monoB_geom)
    
        dimer_energy = energy.get_e_nonbond(dimer_geom, dimer_topology)
        monoA_energy = energy.get_e_nonbond(monoA_geom, monoA_topology)
        monoB_energy = energy.get_e_nonbond(monoB_geom, monoB_topology)
    
        int_energy = dimer_energy - (monoA_energy + monoB_energy)
        reference_energy = test_energies[i]
        TEST_PASS = are_equal(int_energy, reference_energy, 10**test_thresh)

        if (TEST_PASS):
            print('Test %i PASS %15.10f (%s)' % (i + 1, int_energy, tests[i]))
        else:
            print('Test %i FAIL %15.10f (%s)' % (i + 1, int_energy, tests[i]))

    ALL_SUCCESS *= TEST_PASS

    if (ALL_SUCCESS):
        print('All Tests PASS')
    else:
        print('Some Tests FAIL')

    return ALL_SUCCESS
