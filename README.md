TMP Chem Computational Chemistry Source Code
============================================

This repository contains scripts, programs, and data files used with the
"Computational Chemistry" playlist on the [TMP Chem YouTube channel][yt]. 
Primary language is Python3, with toy programs to demonstrate modeling and 
analysis methods. Repository is always subject to change, and no guarantees 
are made of the correctness of output.

[yt]: https://www.youtube.com/tmpchem

Table of contents
-----------------
* [Getting started](#getting-started)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Running the tests](#running-tests)
* [Running the scripts](#running-scripts)
	* [Geometry analysis](#geometry-analysis)
	* [Molecular mechanics](#molecular-mechanics)
* [Author](#author)
* [Acknowledgments](#acknowledgments)

<h2 id="getting-started">Getting started</h2>

Most recent version of project is located in 
[TMP Chem GitHub account][github]. Download by following instructions for 
Git clone from GitHub. Requires a terminal, IDE, etc. to execute Python 
scripts and ability to write to file system within project directories.

[github]: https://www.github.com/tmpchem/computational_chemistry

<h2 id="prerequisites">Prerequisites</h2>

Requires Python 3.5 or greater for script execution. Requires access
to numpy and matplotlib modules. All prerequisites can be met by
downloading and using Python from 
[most recent Anaconda distribution][anaconda].

[anaconda]: https://www.anaconda.com/download/

<h2 id="installation">Installation</h2>

No additional installation necessary after cloning the repository.

<h2 id="running-tests">Running the tests</h2>

WARNING: Test suite is incomplete and subject to change without notice.

To run tests for a project, go to the script directory for that project,
`[top_level_path]/scripts/[project]`. If present, execute `run_tests.py` 
for the project.

    python run_tests.py

This command executes a test suite of unit tests from each test module present
in the mmlib directory. Each test module contains a set of unit tests of methods
within the module, confirming proper behavior and protecting against regression
errors and system misconfigurations.

Once executed, standard output will indicate success with an 'OK' message and
the number of executed unit tests as well as total run time. Any other message
indicates failure and will include the nature of the failing tests.

Open an issue to give feedback or report bugs.

<h2 id="running-scripts">Running the scripts</h2>

As of 16 Feb 2017, repository contains two projects: `geometry_analysis`,
and `molecular_mechanics`.

<h4 id="geometry-analysis">Geometry analysis</h4>

The `geometry_analysis` project contains scripts which take an xyz-format
molecular geometry file as input, and output to screen associated
geometry data, including bond lengths, bond angles, torsion angles,
out-of-plane angles, center of mass, and/or moment of inertia, etc. Sample
xyz files are located in `[top_level_path]/geom/xyz` directory.

<h4 id="molecular-mechanics">Molecular mechanics</h4>

The `molecular_mechanics` project contains scripts to compute molecular
mechanics energy of a system (`mm.py`), molecular dynamics trajectories
(`md.py`), Metropolis Monte Carlo ensembles (`mc.py`), and optimize
molecular coordinates to potential energy minima (`opt.py`).

The energy function and parameters in all cases is based on AMBER FF94.

    Cornell et. al, 
    J. Am. Chem. Soc. 1995, 117, 5179-5197.
    doi.org/10.1021/ja00124a002

Energy function in Equation 1. Atom types in Table 1. Parameter values 
in Table 14. Download AmberTools15 from [here][amber]. After unzipping, 
parameters located in `amber14/dat/leap/parm/parm94.dat`.

[amber]: http://ambermd.org/AmberTools15-get.html

Sample input files for mm are located in
`[top_level_path]/geom/[file_type]` directories, where `file_type]` is
`xyzq` or `prm`. Sample input files for md and mc are located in
`[top_level_path]/geom/sim` directory. Output files for each are
demonstrated in samples, and may be written to any accessible file
name.

<h2 id="author">Author</h2>

The sole author of this package is Trent M. Parker.

- Email: `tmpchemistry@gmail.com`
- LinkedIn: https://www.linkedin.com/in/tmpchem
- Youtube: https://www.youtube.com/tmpchem

<h2 id="acknowledgments">Acknowledgments</h2>

The author wishes to thank the following individuals:
**Dr. Michael S. Marshall**, for encouraging him to learn the Python
language;
**Dr. Lori S. Burns**, for encouraging him to conform to style guidelines
for readable Python code;
**Dr. Michael A. Lewis**, for providing the inspiration for the author to
initiate studies in this field; and
**Dr. C. David Sherrill**, for providing an enormous number of opportunities
to continue to learn and grow as a scientist and a person.
