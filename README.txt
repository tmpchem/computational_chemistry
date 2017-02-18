# Computational Chemistry - by TMP Chem (2014-2017)

This repository contains scripts, programs, and data files used with the
"Computational Chemistry" playlist on the TMP Chem YouTube channel
(youtube.com/tmpchem). Primary language is Python3, with toy programs to
demonstrate modeling and analysis methods. Repository is always subject
to change, and no guarantees are made of the correctness of output.

## Getting Started ##

Most recent version of project is located in TMP Chem GitHub account
(github.com/tmpchem/computational_chemistry). Download by following
instructions for Git clone from GitHub. Requires a terminal, IDE, etc.
to execute Python scripts and ability to write to file system within
project directories.

## Prerequisites ##

Requires Python 3.5 or greater for script execution. Requires access
to numpy and matplotlib modules. All prerequisites can be met by
downloading and using Python from most recent Anaconda package
(www.continuum.io) for appropriate system.

## Installing ##

No additional installation necessary after downloading.

## Running the tests ##

To run tests for a project, go to the script directory for that project
([top_level_path]/scripts/[project]). If present, run `run_tests.py` for
the project.

python run_tests.py

Depending on print level setting in run_tests, each subtest may print
success or failure message (with or without values and reference). If
all tests pass for a function, function receives a test pass. If all
functions pass, the overall unit test receives a pass, and the scripts
are ready to execute. If failed, search recursively for failure source.
Tests are a work in progress, and may not be present or complete.

Bugs or other feedback may be sent via email to `tmpchemistry@gmail.com`.

## Running the scripts ##

As of 16 Feb 2017, repository contains two projects: geometry_analysis,
and molecular_mechanics.

- Geometry Analysis -

The geometry_analysis project contains scripts which take an xyz-format
molecular geometry file as input, and output to screen associated
geometry data, including bond lengths, bond angles, torsion angles,
outofplane angles, center of mass, and/or moment of inertia, etc. Sample
xyz files are located in `[top_level_path]/geom/xyz` directory.

- Molecular Mechanics -

The molecular_mechanics project contain scripts to compute molecular
mechanics energy of a system (mm.py), molecular dynamics trajectories
(md.py), Metropolis Monte Carlo ensembles (mmc.py), and optimize
molecular coordinates to potential energy minima (opt.py).

The energy function and parameters in all cases is based on AMBER FF94
(Cornell et. al, J. Am. Chem. Soc. 1995, 117, 5179-5197.
doi.org/10.1021/ja00124a0002). Energy function in Equation 1. Atom types
in Table 1. Parameter values in Table 14. Download AmberTools15 from
"http://ambermd.org/AmberTools15-get.html". After unzipping, parameters
located in "amber14/dat/leap/parm/parm94.dat".

Sample input files for mm are located in
`[top_level_path]/geom/[file_type]` directories, where [file_type] =
xyzq or prm. Sample input files for md and mmc are located in
`[top_level_path]/geom/sim` directory. Output files for each are
demonstrated in samples, and may be written to any accessible file
name.

## Author ##

The sole author of this package is Trent M. Parker
(tmpchemistry@gmail.com, linkedin.com/in/tmpchem, youtube.com/tmpchem).

## Acknowledgments ##

The author wishes to thank the following individuals:

Dr. Michael S. Marshall - for encouraging him to learn the Python
language.

Dr. Lori S. Burns - for encouraging him to conform to style guidelines
for readable Python code.

Dr. Michael A. Lewis - for providing the inspiration for the author to
initiate studies in this field.

Dr. C. David Sherrill - for providing an enormous number of opportunities
to continue to learn and grow as a scientist and a person.

# end #

