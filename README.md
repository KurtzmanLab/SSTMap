[![Build Status](https://travis-ci.org/choderalab/pymbar.png)](https://travis-ci.org/choderalab/pymbar)
[![Anaconda Cloud Downloads Badge](https://anaconda.org/omnia/pymbar/badges/downloads.svg)](https://anaconda.org/omnia/pymbar)
[![Anaconda Cloud Badge](https://anaconda.org/omnia/pymbar/badges/installer/conda.svg)](https://anaconda.org/omnia/pymbar)
[![PyPI Version](https://badge.fury.io/py/pymbar.png)](https://pypi.python.org/pypi/pymbar)

SSTMap
======

SSTMap is a tool to study structural and thermodynamic properties of water molecules on solute surfaces. It combines grid inhomogeneous solvation theory (IST) with measures of water structure to produce mapping of solvation structure and thermodynamics in protein binding sites and on the surfaces of other molecules of interest, such as small-molecules or host-guest systems. The spatial decomposition of IST integrals allow determination of contribution of speicific regions towards solvation enthalpy and entropy. Alongside thermodynamic information, SSTMap calculates structural properties of water molecules that can aid in developing an understanding of water behvior on the surface and evaluating its displacement.


Installation
------------
__The software is in final stages of release, once released, the following installation commands would work__
The easiest way to install the `SSTMap` release is via [conda](http://conda.pydata.org):
```bash
conda install -c kurtzman_lab sstmap
```
The development version of SSTMap can be download directly from the Github repository and installed as follows:
```bash
git clone git@github.com:KurtzmanLab/SSTMap.git
cd SSTMap
python setup.py install
```
Usage
-----

SSTMap usage involves either unning command-line tools or importing the `sstmap` module. Below are some examples. Stay tuned for tutorials and examples.
```bash
run_hsa -i topology.top -t trajectory.traj -l casp3_ligand.pdb -f 10000 -s 0 -o casp3
run_gist -i topology.top -t trajectory.traj -l casp3_ligand.pdb -d 40 40 40 -f  10000 -s 0 -o casp3
# Replace topology.top and trajectory.traj with files corresponding 
# to the MD package used for simulation, e.g.,
# casp3.prmtop & casp3.netcdf or casp3.gro & casp3.xtc.
```

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
----------------
* Steven Ramsey <vpsramsery@gmail.com>
* Anthony Cruz Balberdi <anthonycruzpr@gmail.com>

References
----------

[1] Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson, J. Chem. Phys. 137, 044101 (2012)

[2] Haider K, Wickstrom L, Ramsey S, Gilson MK and Kurtzman T. Enthalpic Breakdown of Water Structure on Protein Active Site Surfaces. J Phys Chem B. 120:8743-8756, 2016. http://dx.doi.org/10.1021/acs.jpcb.6b01094.

License
-------

`SSTMap` is free software and is licensed under the LGPLv2 license.


Thanks
------
Thanks to Prof. Tom Kurtzman and Prof. Michael Gilson for their guidance in the development of this package.

Notes
-----
