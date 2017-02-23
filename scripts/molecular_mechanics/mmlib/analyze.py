
"""Classes and functions for analyzing and plotting simulation data.

Contains classes to plot dynamic variation of energy and other property
data over the course of a simulation. Also contains classes to compute
expectation values and various other ensemble properties.
"""

import os, sys, math, numpy
import matplotlib.pyplot as plt
from mmlib import fileio, simulate

def point_per_inch():
    """Unit conversion between points and inches."""
    return 72

def percent_image_plot():
    """Fraction of image width which is covered by the plot field."""
    return 0.75

def property_dict():
    """Obtain legend labels, line colors, and priority for properties."""
    pdict = {
        'e_total'   : ['Total',      12, '#000000',  1, ],
        'e_kin'     : ['Kinetic',    11, '#007D34',  2, ],
        'e_pot'     : ['Potential',   1, '#C10020',  3, ],
        'e_nonbond' : ['Non-bonded',  7, '#0000FF',  4, ],
        'e_bonded'  : ['Bonded',      2, '#FF6800',  5, ],
        'e_boundary': ['Boundary',   10, '#551A8B',  6, ],
        'e_vdw'     : ['Vdw',         8, '#00BFFF',  7, ],
        'e_elst'    : ['Elst',        9, '#EEC900',  8, ],
        'e_bond'    : ['Bonds',       3, '#F08080',  9, ],
        'e_angle'   : ['Angles',      4, '#90EE90', 10, ],
        'e_tors'    : ['Torsions',    5, '#FF83FA', 11, ],
        'e_oop'     : ['Outofplanes', 6, '#A9A9A9', 12, ]}
    return pdict

class Plot:
    """Plot class for plotting simulation data.
    
    Base class for all specific plot types used in this module. Derived
    classes include TrajectoryPlot (dynamic energy data variation).
    
    Args:
        ana (mmlib.analysis.Analysis): Analysis object with dynamic
            energy and geometry data.
 
    Attributes:
        simtype (str): Type of molecular simulation:
            `md`: Molecular dynamics.
            `mc`: Metropolis Monte-Carlo.
        data (float**): Dictionary of energy component arrays [kcal/mol].
        plotout (str): File path to output file for plotting image.
        pstart (float): Percent of simulation completed at starting point
            of data analysis (0-100, default=0).
        pstop (float): Analogous value for end point of analysis.
        pdict (list*): Dictionary of energy component plotting parameters.

        figwidth (float): Width of output figure [inches].
        figheight (float): Height of output figure [inches].`
        figsize (float*): Tuple of figure dimensions [inches].
        img_format (str): Format of plotting image file ('pdf' default).
        n_maxpoints (int): Maximum number of data points per component
            in final plotted image.
        n_terms (int): Number of energy component terms to plot.

        line_colors (str*): Array of rgb hex color codes [#rrggbb].
        line_priors (int*): Array of energy component plotting priorities.
        leg_labels (str*): Array of energy component legend labels.
        leg_priors (int*): Array of energy component legend priorities.
        ekeys (str*): Dictionary of energy component string labels.
        xticchars (int*): Dictionary of order-of-magnitude xtic labels.
    """
    def __init__(self, ana):
        self.simtype = ana.simtype
        self.data = ana.prop
        self.plotout = ana.plotout
        self.pstart = ana.percent_start
        self.pstop = ana.percent_stop
        self.pdict = ana.pdict
    
        self.figwidth = 7.5
        self.figheight = 4.5
        self.figsize = (self.figwidth, self.figheight)
        self.img_format = 'pdf'

        ppi = point_per_inch()
        p_im = percent_image_plot()
        self.n_maxpoints = int(math.floor(2.0*ppi*self.figwidth*p_im))
        self.n_terms = len(self.data.keys()) - 1

        self.line_colors = {}
        self.line_priors = {}
        self.leg_labels = {}
        self.leg_priors = {}

        self.ekeys = []
        for key in self.data:
            if ('e_' in key):
                self.line_colors[key] = self.pdict[key][2]
                self.line_priors[key] = self.pdict[key][1]
                self.leg_labels[key] = self.pdict[key][0]
                self.leg_priors[key] = self.pdict[key][3]
                self.ekeys.append(key)
        self.ekeys = sorted(list(self.ekeys), key=lambda x: self.pdict[x][1])
        
        self.ticchars = {0: '', 1: 'k', 2: 'M', 3: 'B', 4: 'T', 5: 'P'}

class TrajectoryPlot(Plot):
    """Plot class for plotting simulation trajectory data.
    
    Molecular dynamics, Monte Carlo, and optimizations produce energy
    and property data at each configuration which is best presented in
    a visual format. This class produces a line graph of that information
    over the course of the trajectory.
    
    To decrease file sizes, once a threshold number of data points is
    reached, the data array is downsampled while trying to maintain an
    equivalent image to plotting all data points.

    Args:
        ana (mmlib.analysis.Analysis): Analysis object with dynamic
            energy and geometry data.

    Attributes:
        name (str): "Name" of simulation from simulation input file name.
        title (str): Title string for plotting image.
        ylabel (str): Y-axis label string for plotting image.
        xlabel (str): X-axis label string for plotting image.
        xvar (str): Data file label for x-axis plotting variable.

        line_width (float): Plotting line width [pts].
        line_alpha (float): Opaqueness of line (0-1).
        line_label (str): Legend label for line (null, handled elsewhere).

        leg_linewidth (float): Linewidth in legend [pts].
        leg_col (int): Number of columns in legend.
        leg_fontsize (float): Font size of legend text [pts].
        leg_pos (float*): Tuple of legend position.
        leg_alpha (float): Opaqueness of legend (0-1).
        leg_framecolor (rgbhex): Color of legend outline (white).

        grid_color (rgbhex): Color of grid lines (black).
        grid_alpha (float): Opaqueness of grid lines (0-1).
        grid_linestyle (str): Style of grid lines (continuous).

        title_fontsize (float): Font size of title text [pts].
        yaxis_fontsize (float): Font size of yaxis label text [pts].
        xaxis_fontsize (float): Font size of xaxis label text [pts].

        xminpad (float): Pad on -x axis bound [relative].
        xmaxpad (float): Pad on +x axis bound [relative].
        yminpad (float): Pad on -y axis bound [relative].
        ymaxpad (float): Pad on +y axis bound [relative].
    """
    def __init__(self, ana):
        Plot.__init__(self, ana)
        self.name = ana.name
        self.ylabel = 'Energy Terms (kcal/mol)'
        if (self.simtype == 'md'):
            self.title = 'Molecular Dynamics Simualtion of '
            self.xlabel = 'Time (ps)'
            self.xvar = 'time'
        elif (self.simtype == 'mc'):
            self.title = 'Metropolis Monte-Carlo Simulation of '
            self.xlabel = 'Configuration Number'
            self.xvar = 'conf'
        self.title += self.name

        self.line_width = 1.0
        self.line_alpha = 1.0
        self.line_label = ''

        self.leg_linewidth = 4.0
        self.leg_col = 4
        self.leg_fontsize = 10.5
        self.leg_pos = (0.985, 1.000)
        self.leg_alpha = 1.0
        self.leg_framecolor = '#FFFFFF'

        self.grid_color = '#000000'
        self.grid_alpha = 0.10
        self.grid_linestyle = '-'        

        self.title_fontsize = 12
        self.yaxis_fontsize = 12
        self.xaxis_fontsize = 12

        self.xminpad = 0.001
        self.xmaxpad = 0.001
        self.yminpad = 0.040
        self.ymaxpad = 0.320

    def make_plot(self):
        """Execute functions needed to plot with appropriate parameters."""
        plt.figure(figsize=self.figsize)
        self.get_lines()
        self.get_axes()
        self.get_xticks()
        self.get_yticks()
        self.get_labels()
        self.get_grid()
        self.get_legend()
        self.output_plot()

    def output_plot(self):
        """Finalize image and save to specified output file."""
        plt.savefig(self.plotout, format=self.img_format)
        plt.close()

    def get_axes(self):
        """Set the boundaries of the axes based on data values."""
        self.get_axis_bounds()
        plt.axis([self.xlow, self.xhigh, self.ylow, self.yhigh])

    def get_labels(self):
        """Set labels to x-axis, y-axis, and title."""
        plt.title(self.title)
        plt.ylabel(self.ylabel)
        plt.xlabel(self.xlabel)

    def get_grid(self):
        """Set grid parameters (transparent gray lines by default)."""
        plt.grid(
            color = self.grid_color,
            alpha = self.grid_alpha,
            linestyle = self.grid_linestyle)

    def get_legend(self):
        """Set top-center legend for each energy component."""
        self.ekeys = sorted(list(self.ekeys), key=lambda x: self.pdict[x][3])
        for key in self.ekeys:
            plt.plot(0, 0,
                linewidth = self.leg_linewidth,
                color = self.line_colors[key],
                alpha = self.leg_alpha,
                label = self.leg_labels[key])
        legend = plt.legend(
            ncol = self.leg_col,
            fontsize = self.leg_fontsize,
            bbox_to_anchor = self.leg_pos,
            framealpha = self.leg_alpha)
        frame = legend.get_frame()
        frame.set_edgecolor(self.leg_framecolor)

    def get_axis_bounds(self):
        """Determine the extrema of x and y axes from data values."""
        self.xmindat = self.data[self.xvar][0]
        self.xmaxdat = self.data[self.xvar][len(self.data[self.xvar])-1]
        self.xrandat = self.xmaxdat - self.xmindat
        self.xmin = self.xmindat + (self.pstart*self.xrandat)/100.0
        self.xmax = self.xmindat + (self.pstop*self.xrandat)/100.0
        self.xran = self.xmax - self.xmin
        self.xlow  = self.xmin - self.xminpad * self.xran
        self.xhigh = self.xmax + self.xmaxpad * self.xran
        self.xplotran = self.xhigh - self.xlow
        self.ymin =  float('inf')
        self.ymax = -float('inf')
        for key in self.ekeys:
            self.ymin = min(self.ymin, numpy.amin(self.data[key]))
            self.ymax = max(self.ymax, numpy.amax(self.data[key]))
        self.yran = self.ymax - self.ymin
        self.ylow  = self.ymin - self.yminpad * self.yran
        self.yhigh = self.ymax + self.ymaxpad * self.yran
        self.yplotran = self.yhigh - self.ylow

    def get_lines(self):
        """Plot lines for each energy component during simulation."""
        self.get_point_indices()
        self.get_xvals()
        self.get_yvals()
        for key in self.ekeys:
            plt.plot(self.xvals, self.yvals[key],
                linewidth=self.line_width,
                color=self.line_colors[key],
                alpha=self.line_alpha,
                label=self.line_label)

    def get_point_indices(self):
        """Compute values needed for data down-sampling algorithm."""
        self.n_confs = len(self.data[self.xvar])
        self.n_start = int(math.floor(self.pstart * self.n_confs)/100.0)
        self.n_stop = int(math.ceil(self.pstop * self.n_confs)/100.0)
        self.n_points = self.n_stop - self.n_start
        self.n_points = min(self.n_points, self.n_maxpoints)
        self.pranges = numpy.linspace(self.n_start, self.n_stop,
            self.n_points + 1)
        for i in range(len(self.pranges)):
            self.pranges[i] = round(self.pranges[i])
        self.pranges = self.pranges.astype(int)

    def get_xvals(self):
        """Compute down-sampled x-axis data points.
        
        Point index function gives a range for each plotted point.
        X-axis values are chosen to be the median within the range.
        """
        self.xvals = numpy.zeros(self.n_points)
        for i in range(self.n_points):
            n_start = self.pranges[max(i-1,0)]
            n_stop = self.pranges[min(i+1, self.n_points)]
            val_array = self.data[self.xvar][n_start:n_stop]
            val = numpy.median(val_array)
            self.xvals[i] = val
    
    def get_yvals(self):
        """Compute down-sampled y-axis data points.
        
        Point index function gives a arnage for each plotted point.
        Y-axis values are chosen to sequentially oscillate between
        minima and maxima of the values within the range.

        This produces a nearly indistinguishable plot from one with all
        data points, due to sampling 2 points for each unit of linewidth.
        """
        self.yvals = {}
        for key in self.ekeys:
            self.yvals[key] = numpy.zeros(self.n_points)
            for i in range(self.n_points):
                start_vals = self.pranges[max(0, i-1):i+1]
                stop_vals  = self.pranges[i:min(i+2, self.n_points+2)]
                n_start = int(math.ceil(numpy.average(start_vals)))
                n_stop = int(math.ceil(numpy.average(stop_vals)))
                val_array = self.data[key][n_start:n_stop]
                if (i%2 == 0):
                    val = numpy.amax(val_array)
                else:
                    val = numpy.amin(val_array)
                self.yvals[key][i] = val

    def get_xticks(self):
        """Determine array of X-axis tick mark values."""
        self.xdig = int(math.floor(math.log10(self.xplotran)))
        self.xlead = self.xran * 10**(-self.xdig)
        self.xbase = 10.0
        if (self.xlead <= 2.0):
            self.xbase = 2.0
        elif (self.xlead <= 5.0):
            self.xbase = 5.0
        self.xres = self.xbase * 10**(self.xdig - 1)
        self.xdelta = 0.1 * self.xres
        xminarr = self.xres * math.floor(self.xmin / self.xres)
        xmaxarr = self.xres * math.ceil(self.xmax / self.xres)
        xmaxarr += self.xmaxpad * (xmaxarr - xminarr)
        self.xtics = numpy.arange(xminarr, xmaxarr, self.xres)
        self.xticlabels = ['' for i in range(len(self.xtics))]
        for i in range(len(self.xtics)):
            exp = int(math.floor(math.log10(max(1, abs(self.xtics[i]))))/3)
            val = self.xtics[i] / 10**(3*exp)
            if (int(val) == val):
                val = int(val)
            self.xticlabels[i] = '%s%s' % (val, self.ticchars[exp])
        plt.xticks(self.xtics, self.xticlabels)

    def get_yticks(self):
        """Determine array of Y-axis tick mark values."""
        self.ydig = int(math.floor(math.log10(self.yplotran)))
        self.ylead = (self.yhigh - self.ylow) * 10**(-self.ydig)
        self.ybase = 10.0
        if (self.ylead <= 2.0):
            self.ybase = 2.0
        elif (self.ylead <= 5.0):
            self.ybase = 5.0
        self.yres = self.ybase * 10**(self.ydig - 1)
        self.ydelta = 0.1 * self.yres
        yminarr = self.yres * math.floor(self.ymin / self.yres)
        ymaxarr = self.yres * math.ceil(self.ymax / self.yres)
        ymaxarr += self.ymaxpad * (ymaxarr - yminarr)
        self.ytics = numpy.arange(yminarr, ymaxarr, self.yres)
        self.yticlabels = ['' for i in range(len(self.ytics))]
        for i in range(len(self.ytics)):
            exp = int(math.floor(math.log10(max(1, abs(self.ytics[i]))))/3)
            val = self.ytics[i] / 10**(3*exp)
            if (int(val) == val):
                val = int(val)
            self.yticlabels[i] = '%s%s' % (val, self.ticchars[exp])
        plt.yticks(self.ytics, self.yticlabels)
    
class Analysis:
    """Analysis class for computing ensemble simulation data.
    
    Molecular dynamics, Monte Carlo, and optimizations produce dynamic
    data at a large number of configurations. This module contains many
    functions for processing data and computing ensemble properties.
    
    Args:
        infile_name (str): Path (relative or absolute) to input file
            from completed md or mc simulation.
    
    Attributes:
        infile (str): Absolute path to input file with input parameters.
        indir (str): Absolute path to directory of input file.
        pdict (list*): Dictionary with energy component data.
        simtype (str): Simulation type, options:
            `md`: Molecular dynamics.
            `mc`: Metropolis Monte Carlo.
        simdir (str): File path to simulation input file.
        plotout (str): File path to trajectory plot output file.
        energyin (str): File path to input file with energy trajectory data.
        geomin (str): File path to input file with geometry trajectory data.
        percent_start (float): Starting point of data analysis relative to
            number of simulations (0-100, default = 0).
        percent_stop (float): Analogous stopping point (default = 100).
        name (str): "name" of simulation, from simulation input file name.
        prop (float**): Dictionary of energy component data arrays.
        traj (float***): Array of geometry arrays / coordinate trajectory.
        tplt (mmlib.analyze.TrajectoryPlot): TrajectoryPlot object with
            dynamic trajectory energy data for plotting to file.    
        eavg (float*): Dictionary of energy component averages [kcal/mol].
        estd (float*): Ibid for standard deviations.
        emin (float*): Ibid for minima.
        emax (float*): Ibid for maxima.
    """
    def __init__(self, infile_name):
        self.infile = os.path.realpath(infile_name)
        self.indir = os.path.dirname(self.infile)
        self.pdict = property_dict()

        self.simtype = ''
        self.simfile = ''
        self.simdir = ''
        self.plotout = 'plot.pdf'
        self.energyin = ''
        self.geomin = ''
        self.percent_start = 0.0
        self.percent_stop = 100.0

        self.read_in_data()

    def read_in_data(self):
        """Read in input parameters from input file."""
        fileio.get_analysis_data(self)
        self.read_in_files()
        self.read_in_prop()

    def read_in_files(self):
        """Read in property and trajectory file names from input file."""
        cwd = os.getcwd()
        os.chdir(self.simdir)
        self.sim = simulate.Simulation(self.simfile, self.simtype)
        os.chdir(cwd)
        self.energyin = self.sim.energyout
        self.geomin = self.sim.geomout
        self.name = self.sim.mol.name
        delattr(self, 'sim')

    def read_in_prop(self):
        """Read in molecular properties from simulation energy file."""
        self.prop = fileio.get_properties(self.energyin)
    
    def read_in_geom(self):
        """Read in molecular trajectory from simulation geometry file."""
        self.traj = fileio.get_trajectory(self.geomin)

    def run_analysis(self):
        """Plot dynamic energy data and print expectation values."""
        self.tplt = TrajectoryPlot(self)
        self.tplt.make_plot()
        self.get_energy_stats()
        fileio.print_averages(self)

    def get_energy_stats(self):
        """Compute average, stdev, min, and max of each energy term."""
        self.eavg = {}
        self.estd = {}
        self.emin = {}
        self.emax = {}
        for key in self.pdict:
            if (key in self.prop):
                keydata = self.prop[key]
                self.eavg[key] = numpy.average(keydata)
                self.estd[key] = numpy.std(keydata)
                self.emin[key] = numpy.amin(keydata)
                self.emax[key] = numpy.amax(keydata)

# end of module

