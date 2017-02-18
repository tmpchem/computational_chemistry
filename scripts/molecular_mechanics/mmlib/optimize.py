
"""Classes and functions for optimizing molecular coordinates."""

import os, sys, math, numpy
from mmlib import energy, gradient, fileio

class Trajectory:
    """Trajectory class for optimization step history.
    
    Keeps an array each for the history of the total energy, total
    energy gradient, and molecular coordinates at each step of an
    optimization procedure.
    
    Args:
        mol (mmlib.molecule.Molecule): Molecule object with molecular
            geometry, energy, and gradient data.
    
    Attributes:
        n_atoms (int): number of atoms in molecule.
        n_steps (int): number of steps in optimization history.
        energy (float*): array of total energy history.
        coords (float***): array of molecular coordinate history.
        grad (float***): array of total energy gradient history.
    """
    def __init__(self, mol):
        self.n_atoms = mol.n_atoms
        self.n_steps = 0
        self.energy = []
        self.coords = []
        self.grad = []
        self.append_step(mol)
    """Append current molecule data to Trajectory object."""
    def append_step(self, mol):
        self.n_steps += 1
        self.coords.append(numpy.zeros((self.n_atoms, 3)))
        self.grad.append(numpy.zeros((self.n_atoms, 3)))
        self.energy.append(mol.e_total)
        for i in range(self.n_atoms):
            for j in range(3):
                self.coords[-1][i][j] = mol.atoms[i].coords[j]
                self.grad[-1][i][j] = mol.g_total[i][j]

class Optimization:
    """Optimization class for finding extremes of molecular energy.
    
    Propogates molecular geometry towards total energy local minimum
    based on total energy gradient. Convergence is reached when a set
    number of optimization criteria are reached. Histories of energies,
    gradients, and coordinates are kept for printing and adjustments.
    
    Args:
        infile_name (str): Path of optimization input file. May be
            relative or absolute, though absolute is safer.
    
    Attributes:
        infile(str): Input file (see Args).
        indir (str): Directory of input file.
        mol (mmlib.molecule.Molecule): Molecule object with molecular
            geometry, energy, and gradient data.
        traj (mmlib.optimize.Trajectory): Trajectory object with energy,
            coordinate, and gradient history.
        opt_type (str): Optimization algorithm:
            `sd`: Steepest descent.
            `cg`: Conjugate gradient.
            `nr`: Newton-Raphson.
        opt_str (float* or str): Criteria for optimization
            convergence. May either be explicitly enumerated array (float):
                [delta_e, grad_max, grad_rms, disp_max, disp_rms]
            or standard reference values from dictionary (str):
                (see ref_opt_criteria method).
        geomout (str): Geometry printing output file path.
        energyout (str): Energy printing output file path.

        delta_e (float): Change in total energy [kcal/mol] from previous
            to current geometry configuration.
        grad_rms (float): Root-mean-square of total energy gradient
            matrix [kcal/(mol*A)].
        grad_max (float): Maximum total energy gradient matrix element
            [kcal/(mol*A)].
        disp_rms (float): Root-mean-square of coordinate displacement
            matrix from previous to current configuration [Angstrom].
        disp_max (float): Maximum total coordinate displacement matrix
            element [Angstrom].
        conv_delta_e (float): `delta_e` needed for convergence.
        conv_grad_rms (float): `grad_rms` needed for convergence.
        conv_grad_max (float): `grad_max` needed for convergence.
        conv_disp_rms (float): `disp_rms` needed for convergence.
        conv_disp_max (float): `disp_max` needed for convergence.

        must_converge (bool*): Array denoting which of the 5 convergence
            criteria must be met for success.
        are_converged (bool*): Array denoting which of the 5 convergence
            criteria have been met.
        is_converged (bool): Whether sufficient convergence criteria
            have been met.
        n_maxiter (int): Largest possible number of optimization iterations.
        n_iter (int): Current number of optimization iterations.
        n_subiter (int): Number of subiterations within current iteration.
        disp_mag (float): Scalar multiple [Angstrom] of gradient vector
            for displacement step.
        disp_deriv (float): Numerical energy derivative in direction of
            specified displacement vector [kcal/(mol*A)].
        ccoords (float**): Array with copy of molecular coordinates
            to revert from displacement [Angstrom].
    """
    def __init__(self, infile_name):
        self.infile = os.path.realpath(infile_name)
        self.indir = os.path.dirname(self.infile)
        self.name = '.'.join(self.infile.split('/')[-1].split('.')[:-1])
        self.opt_type = 'sd' # steepest descent by default
        self.opt_str = ''
        self.mol = []
        self.geomout = 'geom.xyz'
        self.energyout = 'energy.dat'

        self.delta_e = float('inf')
        self.grad_rms = float('inf')
        self.grad_max = float('inf')
        self.disp_rms = float('inf')
        self.disp_max = float('inf')

        self.conv_delta_e  = 1.0 * 10**-6
        self.conv_grad_rms = 1.0 * 10**-4
        self.conv_grad_max = 2.0 * 10**-4
        self.conv_disp_rms = 1.0 * 10**-3
        self.conv_disp_max = 2.0 * 10**-3

        self.must_converge = [True for i in range(5)]
        self.are_converged = [False for i in range(5)]
        self.is_converged = False
        self.n_maxiter = 200
        self.n_iter = 0
        self.n_subiter = 0

        self.disp_mag = 1.0 * 10**-4
        self.disp_deriv = 0.0

        self.read_in_data()
        self.get_opt_criteria()
        self.update_energy()
        self.update_gradient()
        self.traj = Trajectory(self.mol)

    def read_in_data(self):
        """Read in optimization data from input file."""
        fileio.get_opt_data(self)

    def optimize(self):
        """Displace molecule to minimum energy molecular coordinates."""
        self.open_output_files()
        while (self.n_iter < self.n_maxiter and not self.is_converged):
            self.n_iter += 1
            self.choose_step_direction(self.opt_type)
            self.line_search(-1.0 * self.step_dir)
            self.update_energy()
            self.update_gradient()
            self.traj.append_step(self.mol)
            self.update_criteria()
            self.check_convergence()
            self.print_status()
        self.close_output_files()

    def choose_step_direction(self, opt_type):
        """Choose step direction based on energy minimization algorithm.

        Args:
            opt_type (str): Specific optimization algorithm.     
                `sd` (steepest descent): Travel along gradient.
                `cg` (conjugate gradient): Improve gradient using gradient
                    history.
        """
        if   (opt_type == 'sd'):
            self.get_sd_step_dir()
        elif (opt_type == 'cg'):
            self.get_cg_step_dir()
        else:
            print('Error: optimization type (%s) not recognized!' % (opt_type))
            sys.exit()

    def get_sd_step_dir(self):
        """Steepest descent optimization step direction vector."""
        self.step_dir = self.mol.g_total

    def get_cg_step_dir(self):
        """Conjugate gradient optimization step direction vector."""
        if (self.n_iter <= 1):
            self.hvec = self.mol.g_total
            gamma = 0.0
        else:
            v1 = self.traj.grad[-1] - self.traj.grad[-2]
            v1 = v1.reshape((1, 3*self.mol.n_atoms))
            v2 = self.traj.grad[-1].reshape((3*self.mol.n_atoms, 1))
            gamma  = numpy.linalg.norm(numpy.dot(v1, v2))
            gamma *= 1.0 / numpy.linalg.norm(self.traj.grad[-1])**2
            self.hvec = self.mol.g_total + gamma * self.hvec
        self.step_dir = self.hvec

    def update_criteria(self):
        """Update values of the 5 optimization convergence criteria."""
        grad = self.traj.grad[-1]
        disp = self.traj.coords[-1] - self.traj.coords[-2]
        self.delta_e = self.traj.energy[-1] - self.traj.energy[-2]
        self.grad_max = numpy.amax(grad)
        self.disp_max = numpy.amax(disp)
        self.grad_rms = math.sqrt(numpy.mean(grad**2))
        self.disp_rms = math.sqrt(numpy.mean(disp**2))

    def check_convergence(self):
        """Determine if all required convergence criteria are met.
        
        Each of 5 values (energy change, rms and max gradient, rms and
        max displacement) can set to False if they must converge for
        success.
        """
        self.is_converged = True
        self.are_converged[0] = (abs(self.delta_e) < self.conv_delta_e)
        self.are_converged[1] = (self.grad_rms < self.conv_grad_rms)
        self.are_converged[2] = (self.grad_max < self.conv_grad_max)
        self.are_converged[3] = (self.disp_rms < self.conv_disp_rms)
        self.are_converged[4] = (self.disp_max < self.conv_disp_max)
        for i in range(5):
            if (self.must_converge[i] and not self.are_converged[i]):
                self.is_converged = False

    def line_search(self, disp_vector):
        """Perform line search along vector to find numerical minimum.
        
        Displacement vector provides 1d potential energy surface
        manifold. This function displaces the molecular coordinates
        to the nearest local minimum on that surface through binary
        searches.
        
        Args:
            disp_vector (float**): Displacement vector for all 3N atomic
                coordinates.
        """
        self.get_disp_deriv(self.disp_mag, disp_vector)
        disp_mag = self.disp_mag
        disp_sign = 1.0 if (self.disp_deriv <= 0.0) else -1.0
        disp_mag *= disp_sign
        disp_sign_same = True
        self.n_subiter = 0
        # binary search tree to find upper bound on displacement magnitude
        while (disp_sign_same):
            self.n_subiter += 1
            self.displace_coords(+1.0 * disp_mag, disp_vector)
            self.get_disp_deriv(disp_mag, disp_vector)
            self.displace_coords(-1.0 * disp_mag, disp_vector)
            old_disp_sign = disp_sign
            disp_sign = 1.0 if (self.disp_deriv <= 0.0) else -1.0
            disp_sign_same = bool(disp_sign == old_disp_sign)
            disp_mag *= 2.0
        self.get_disp_deriv(disp_mag, disp_vector)
        self.adjust_disp_mag(self.n_subiter)
        numer = 1.0
        denom = 2.0
        # binary search tree to find value of displacement within bounds
        for i in range(7):
            self.n_subiter += 1
            test_disp = disp_mag * numer / denom
            self.displace_coords(+1.0 * test_disp, disp_vector)
            self.get_disp_deriv(disp_mag / (2**(-i)), disp_vector)
            self.displace_coords(-1.0 * test_disp, disp_vector)
            direc = 1.0 if (self.disp_deriv < 0.0) else -1.0
            numer = 2*numer + direc
            denom = 2*denom
        disp_mag *= numer / denom
        # final line search energy minimized molecular coordinates
        self.displace_coords(+1.0 * disp_mag, disp_vector)

    def get_opt_criteria(self):
        """Dictionary of reference values of convergence criteria.

        The opt_str member is set to null string by default, and doesn't
        execute the conditional. If overridden from input file to value
        in dictionary, reset all 5 convergence criteria to specified
        value."""
        print(self.opt_str)
        opt_criteria_refs = {
            'loose'     : [1.0 * 10**(-4), 1.0 * 10**(-3), 2.0 * 10**(-3),
                           1.0 * 10**(-2), 2.0 * 10**(-2)],
            'default'   : [1.0 * 10**(-6), 1.0 * 10**(-4), 2.0 * 10**(-4),
                           1.0 * 10**(-3), 2.0 * 10**(-3)],
            'tight'     : [1.0 * 10**(-8), 1.0 * 10**(-5), 2.0 * 10**(-5),
                           1.0 * 10**(-4), 2.0 * 10**(-4)],
            'verytight' : [1.0 * 10**(-9), 1.0 * 10**(-6), 2.0 * 10**(-6),
                           1.0 * 10**(-5), 2.0 * 10**(-5)]}
        if (self.opt_str in opt_criteria_refs):
            opt_vals = opt_criteria_refs[self.opt_str]
            self.conv_delta_e  = opt_vals[0]
            self.conv_grad_rms = opt_vals[1] 
            self.conv_grad_max = opt_vals[2] 
            self.conv_disp_rms = opt_vals[3] 
            self.conv_disp_max = opt_vals[4] 

    def get_test_energy(self, test_coords):
        """Determine energy of molecule at set of test coordinates.
        
        Args:
            test_coords (float**): Array of atomic coordinates [Angstrom]
                at which to determine molecular energy.
        """
        self.copy_coords()
        self.update_coords(test_coords)
        self.update_energy()
        self.test_energy = self.mol.e_total
        self.update_coords(self.ccoords)

    def displace_coords(self, disp_mag, disp_vector):
        """Dispace coordinates by disp_mag in disp_vector direction.
        
        Args:
            disp_mag (float): Magnitude of displacement [Angstrom].
            disp_vector (float**): Displacement direction vector.
        """
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].coords[j] += disp_mag * disp_vector[i][j]
    def copy_coords(self):
        """Create a copy of current molecular coordinates."""
        self.ccoords = numpy.zeros((self.mol.n_atoms, 3))
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.ccoords[i][j] = self.mol.atoms[i].coords[j]
    
    def get_disp_deriv(self, disp_mag, disp_vector):
        """Numerical energy derivative in direction of displacement."""
        num_disp = 0.01 * disp_mag
        self.update_energy()
        e_neg_total = self.mol.e_total
        self.displace_coords(+1.0 * num_disp, disp_vector)
        self.update_energy()
        e_pos_total = self.mol.e_total
        self.disp_deriv = (e_pos_total - e_neg_total) / num_disp
        self.displace_coords(-1.0 * num_disp, disp_vector)

    def adjust_disp_mag(self, n_subiter):
        """Adjust starting guess for line search displacement.
        
        Args:
            n_subiter (int): Previous number of initial of doublings
                needed to switch the sign of the line derivative.
        """
        if (n_subiter == 1):
            self.disp_mag *= 0.7
        else:
            self.disp_mag *= 1.4

    def update_energy(self):
        """Update energy at current molecular coordinates."""
        self.mol.get_energy('nokinetic')

    def update_gradient(self):
        """Update energy gradient at current molecular coordinates."""
        self.mol.get_gradient('analytic')

    def update_coords(self, new_coords):
        """Update atomic coordinates to values in a given vector."""
        for i in range(self.mol.n_atoms):
            for j in range(3):
                self.mol.atoms[i].coords[j] = new_coords[i][j]

    def print_energy_header(self):
        """Print header of convergence output columns to file."""
        self.efile.write('# iter          energy    delta_e    grad_max')
        self.efile.write('    grad_rms    disp_max    disp_rms\n')
        self.efile.write('#                      ')
        self.efile.write('%10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n' % (
            self.conv_delta_e, self.conv_grad_max, self.conv_grad_rms,
            self.conv_disp_max, self.conv_disp_rms))

    def print_energy(self, n_iter):
        """Print convergence information to file."""
        e = self.efile
        t = self.traj
        grad = t.grad[n_iter]
        disp = t.coords[n_iter] - t.coords[max(0, n_iter-1)]
        delta_e = t.energy[n_iter] - t.energy[max(0, n_iter-1)]
        gmax = numpy.amax(grad)
        dmax = numpy.amax(disp)
        grms = math.sqrt(numpy.mean(grad**2))
        drms = math.sqrt(numpy.mean(disp**2))
        conv_str = [' ' for i in range(5)]
        if (n_iter > 0):
            conv_str[0] = '*' if (abs(delta_e) < self.conv_delta_e) else ' '
            conv_str[1] = '*' if (gmax < self.conv_grad_max) else ' '
            conv_str[2] = '*' if (grms < self.conv_grad_rms) else ' '
            conv_str[3] = '*' if (dmax < self.conv_disp_max) else ' '
            conv_str[4] = '*' if (drms < self.conv_disp_rms) else ' '
        pstr  = '%3i' % (n_iter)
        pstr += ' %18.12f' % (t.energy[n_iter])
        pstr += ' %10.3e%s' % (delta_e, conv_str[0])
        pstr += ' %10.3e%s' % (gmax, conv_str[1])
        pstr += ' %10.3e%s' % (grms, conv_str[2])
        pstr += ' %10.3e%s' % (dmax, conv_str[3])
        pstr += ' %10.3e%s' % (drms, conv_str[4])
        e.write('%s\n' % (pstr))

    def print_status(self):
        """Print xyz-format geometry of system to trajectory file."""
        comment = '%s, iteration %i' % (self.name, self.n_iter)
        fileio.print_coords_file(self.gfile, self.mol, comment)
        self.print_energy(self.n_iter)
        self.gfile.flush()
        self.efile.flush()

    def open_output_files(self):
        """Open output files for convergence and geometry data printing."""
        self.gfile = open(self.geomout, "w")
        self.efile = open(self.energyout, "w")
        self.print_energy_header()

    def close_output_files(self):
        """Close output files for convergence and geometry data printing."""
        self.gfile.close()
        self.efile.close()

# end of module

