[![Build Status](https://travis-ci.org/kamran-haider/SSTMap.svg?branch=master)](https://travis-ci.org/kamran-haider/SSTMap)
[![Anaconda-Server Badge](https://anaconda.org/solvationtools/sstmap/badges/installer/conda.svg)](https://conda.anaconda.org/solvationtools)
[![Anaconda-Server Badge](https://anaconda.org/solvationtools/sstmap/badges/downloads.svg)](https://anaconda.org/solvationtools/sstmap)
[![DOI for Citing SSTMap](https://img.shields.io/badge/DOI-10.1021%2Fj.jctc.2017.11.021-blue.svg)](http://doi.org/10.1021/acs.jctc.7b00592)

Welcome to SSTMap Release version 1.1.4
Thank you for using SSTMap, we've gone through a great deal of work trying to ensure compatibility with various operating systems and MD packages though we have not covered all the possible combinations.  If you run into a bug, please report it on the issues and we will work to quickly resolve it.

Thank you,

Kamran and the Development Team.


SSTMap
======

SSTMap is a tool to study structural and thermodynamic properties of water molecules on solute surfaces. It combines grid inhomogeneous solvation theory (IST) with measures of water structure to produce mapping of solvation structure and thermodynamics in protein binding sites and on the surfaces of other molecules of interest, such as small-molecules or host-guest systems. It provides both site-based and grid-based calculations in one package, with support for multiple MD packages and can be integrated into Python’s scientific computing environment for advanced applications. The cross-platform support is enabled by trajectory and topology parsers of <a href="http://mdtraj.org">MDTraj</a> and <a href="http://parmed.github.io/ParmEd/html/index.html">ParmEd</a>.

Installation
------------
### Conda Installation
The anaconda python distribution comes with `conda` package manager, which can then be used to install `SSTMap` with the following commands.

```
conda config --add channels omnia
conda config --add channels solvationtools
conda install sstmap
```

### GitHub Source Code
You can install the latest development version from the GitHub repository by executing

```
pip install git+git://github.com/kurtzmanlab/sstmap.git@1.1.4#egg=sstmap
```

You can also download the release package manually from GitHub, unzip it, navigate into the directory, and execute the command:

```bash
https://github.com/KurtzmanLab/SSTMap/archive/1.1.4.tar.gz or https://github.com/KurtzmanLab/SSTMap/archive/1.1.4.zip
tar -xvf 1.1.4.tar.gz or unzip 1.1.4.zip
cd SSTMap-1.1.4
python setup.py install
```
Or you can clone the GitHub repository, navigate into the directory, and execute the command: 

```bash
git clone git@github.com:KurtzmanLab/SSTMap.git
cd SSTMap
```
**For the release version:**
```
git checkout tags/1.1.4
python setup.py install
```
**For the developmental version:**

```
python setup.py install
```
When building from the source code or using the release package, make sure that you manually install the dependencies: `mdtraj` and `parmed`. You can do this by:
```
conda config --add channels omnia
conda install mdtraj parmed
``` 
Or using pip:
```
pip install ParmEd 
pip install mdtraj 
```  
Usage
-----

SSTMap provides command-line tools for hydration site analysis (HSA) and Grid Inhomogeneous Solvation Theory (GIST), `run_hsa` and `run_gist`, respectively. The functionality of SSTMap is also available through its Python API, available as the `sstmap` module. For more details, please visit [sstmap.org](sstmap.org).

An example of running a quick HSA and GIST calculation on a test system available in [sstmap_test_suite](https://github.com/KurtzmanLab/sstmap_test_suite). You can download the full test suite from [here](https://www.dropbox.com/sh/hrijgk8n5z12bgi/AABSigcBf9PN_7-Z26VCCPePa?dl=0) (since Github repository doesn't contain trajectory files). For a given platform, `cd` to its sub-directory and run the commands as shown below. 
```bash
cd sstmap_test_suite/platforms/amber
$ run_hsa -i testcase.prmtop -t md100ps.nc -l ligand.pdb -s 0 -f 100 -o testcase
$ run_gist -i testcase.prmtop -t md100ps.nc -l ligand.pdb -g 20 20 20 -s 0 -f 100 -o testcase
```
For examples using MD simulations generated from other packages, such as [Amber](http://ambermd.org/), [Charmm](https://www.charmm.org), [Gromacs](http://www.gromacs.org/), [NAMD](http://www.ks.uiuc.edu/Research/namd/), [OpenMM](http://openmm.org/) and [Desmond](https://www.deshawresearch.com/resources_desmond.html), please follow [this tutorial](http://sstmap.org/2017/05/03/simple-examples/) on [sstmap.org](sstmap.org). SSTMap can also be used as a Python module:

```python
import sstmap as sm
# Example 1: Run a grid-based calculation
# with all quantities.
gist = sm.GridWaterAnalysis(
        "casp3.prmtop", 
        "casp3.netcdf",
        ligand_file="casp3_ligand.pdb", 
        grid_dimensions=[40, 40, 40], 
        prefix="casp3")
gist.calculate_grid_quantities()
# Example 2: Run gist with only energy calculations.
gist.calculate_grid_quantities(energy=True)
# Example 3: Initialize gist with a grid center position.
gist = sm.GridWaterAnalysis(
        "casp3.prmtop", 
        "casp3.netcdf",
        grid_center=[35.33, 52.23, 54.96], 
        grid_dimensions=[40, 40, 40], 
        prefix="casp3")
```

Principal Developer(s)
---------------------
* Kamran Haider <kamranhaider.mb@gmail.com>

Co-Developers
-------------
* Steven Ramsey <vpsramsey@gmail.com>
* Anthony Cruz Balberdy <anthonycruzpr@gmail.com>
* Tobias Wulsdorf <tobias.wulsdorf@pharmazie.uni-marburg.de>

Principal Investigators
---------------------
* Tom Kurtzman
* Michael Gilson

References
----------
Please cite the following when you use this software.  

[1] Haider K, Cruz A, Ramsey S, Gilson M. and Kurtzman T. Solvation Structure and Thermodynamic Mapping (SSTMap): An Open-Source, Flexible Package for the Analysis of Water in Molecular Dynamics Trajectories. J. Chem. Theory Comput. (10.1021/acs.jctc.7b00592) 2017.
[2] Crystal N. Nguyen, Michael K. Gilson, Tom Young. Structure and Thermodynamics of Molecular Hydration via Grid Inhomogeneous Solvation Theory. eprint arXiv:1108.4876, (2011).

[3] Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson. Grid inhomogeneous solvation theory: hydration structure and thermodynamics of the miniature receptor cucurbit[7]uril. J. Chem. Phys. 137, 044101 (2012)

[4] Haider K, Wickstrom L, Ramsey S, Gilson MK and Kurtzman T. Enthalpic Breakdown of Water Structure on Protein Active Site Surfaces. J Phys Chem B. 120:8743-8756, (2016). http://dx.doi.org/10.1021/acs.jpcb.6b01094.

[5] Themis Lazaridis. Inhomogeneous Fluid Approach to Solvation Thermodynamics. 1. Theory. The Journal of Physical Chemistry B 102 (18), 3531-3541, (1998). DOI: 10.1021/jp9723574


License
-------

`SSTMap` is free software and is licensed under the MIT license.


Acknowledgements
--------
This project is funded through the NIH Grant: R01-GM100946 

