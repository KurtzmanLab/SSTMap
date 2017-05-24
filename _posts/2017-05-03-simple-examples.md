---
layout: post
author: kamran-haider
title: Running SSTMap Calculations
date: 2017-05-03
published: true
---

`SSTMap` provides two main approaches for mapping water structure and thermodynamics on to solute surfaces, such as protein binding sites, the hydration site analysis (HSA) and Grid Inhomogeneous Solvation Theory (GIST). The theory behind these approaches are described in several publications (see References). Here we provide selected examples and discuss some practical aspects of running HSA and GIST calculations, through `run_hsa` and `run_gist`, respectively, which are the main command-line tools in `SSTMap`. For a detail list of command-line arguments to these programs and the MD trajectory requirements, see the bottom of this page. This tutorial focuses the commands for running calculations. For a detailed description of the outputs genrated by these programs, see this [post](http://sstmap.org/2017/05/09/undestanding-output/). 
<!--more-->
### Amber
```bash
$ run_hsa -i testcase.prmtop -t md100ps.nc -l ligand.pdb -s 0 -f 100 -o testcase

$ run_gist -i testcase.prmtop -t md100ps.nc -l ligand.pdb -g 20 20 20 -s 0 -f 100 -o testcase
```
Since amber prmtop file contains both non-bonded parameters and topology information for the system, therefore, you can leave out the `-p` flag, which is used to specify additional parameter files.
### Charmm/NAMD
The `-p` flag is mandatory for (Charmm, NAMD, Gromacs). For these package, in order to obtain non-bonded parameters for energy calculation, additional files are required by the underlying `Parmed` parsers. For example, the commands for a Charmm simulation would be as follows: 
```
run_hsa -i testcase.psf -t md100ps.dcd -l ligand.pdb -p toppar/ -s 0 -f 100 -o testcase
$ run_gist -i testcase.psf -t md100ps.dcd -l ligand.pdb -g 10 10 10 -p toppar/ -s 0 -f 100 -o testcase
```

When not sure of what additional files to provide, just run the `run_hsa` or `run_gist` on your system, without the `-p` flag and the help message will let you know what files to provide, e.g.,

```
$ run_hsa -i testcase.psf -t md100ps.dcd -l ligand.pdb -s 0 -f 100 -o testcase
    SSTMap requires toppar as a supporting file/data for psf parameter format. Please provide it as an argument to supporting_file argument or if you are running run_gist or run_hsa, provide it as an argument to -p flag.
    More specifically, Please provide a folder named toppar that contains charmm parameter/topology files.
```   
NAMD uses the same forcefileds as Charmm. Hence the same commands, as shown for as Charmm, are valid for NAMD.

### Gromacs
```
$ run_hsa -i testcase.gro -t md100ps.xtc -l ligand.pdb -p params.top -s 0 -f 100 -o testcase
$ run_gist -i testcase.gro -t md100ps.xtc -l ligand.pdb -p params.top -g 20 20 20 -s 0 -f 100 -o testcase
```
### OpenMM
OpenMM is an MD toolkit that allows running simulations using Amber99SB forcefield or custom forcefields (written as xml file formats as described [here](http://docs.openmm.org/7.1.0/userguide/application.html#creating-force-fields)). Indeed, systems can also be build from coordinates and topology files for Charmm and Gromacs. Depending on which of the three forcefields was used to build the system for the OpenMM simulation, the correspnding `run_hsa` and `run_gist` commands for the packages will be valid.

### Desmond/Others

The desmond `cms` file format, which contains forcefield and topology information for simulations run in Desmond, is not supported by `ParmEd`, the underlying parameter and topology parser program. However, with three additional steps, it is possible to use Desmond simulations for SSTMap calculations.

* <strong>Step 1</strong>: Convert `cms` file into `pdb` file. This can be done in Schrodinger's Maestro or VMD.

* <strong>Step 2</strong>: Convert dtr file into a netcdf file. Although MFTraj allows reading dtr file format but the cross-platform trajectory loaders of mdtraj (e.g., `mdtraj.load()` or `mdtraj.load_frame()`) do not accept dtr file format. Only a specific loader `mdtraj.load_dtr()` is available for this purpose. In SSTMap code, the trajectory reading is only done through MDTraj's cross-platform loaders. As a result, to process Desmond simulation through SSTMap it is necessary to convert them into netcdf format (although alternative format, dcd, xtc, H5 etc could also be used but netcdf is recommended).Th dtr to netcdf conversion can be be done in VMD using a Tcl script or using a Python script that uses MDTraj's dtr loader. An example Tcl script for this purpose `dtr_to_netcdf.tcl` is shown below:

```tcl
mol new input.cms type mae first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 1traj/clickme.dtr type dtr first 0 last -1 step 10 filebonds 1 autobonds 1 waitfor all
animate write output.nc beg 1 end -1
```

This cane be run on the command-line as:

```bash
$ vmd -dispdev text -eoftext < dtr_to_netcdf.tcl
```

Similarly, an example Python script for tihs purposes, `dtr_to_netcdf.py` is given below:

```python
traj = md.load_dtr("traj/clickme.dtr", top=input.pdb)
traj.save_netcdf("converted.nc")
```
* <strong>Step 3</strong>: Finally, extract non-bonded parameters for each atom in the system. This can be done by running an auxiliary Python script provided in SSTMap repository under [scripts](https://github.com/KurtzmanLab/SSTMap/tree/master/sstmap/scripts) folder: `desmond_extract_nbparams.py`. This script uses Schrodinger Suite's Python API to parse cms file and extract non-bonded parameters. In ordder to run this script, make sure $SCHRODINGER enviornment variable is set on your terminal (which is set by default when a Desmond or Schrodinger installation is done.). For example:
```bash
$SCHRODINGER/run desmond_extract_nbparams.py input.cms
```
will generate a text file called `inputcms_nb_parms.txt`. This text file is a $$N \times 3$$ matrix of non-bonded prameters, where N is the number of atoms and three columns correspond to sigma, epsilon and charge parameters for each atom in the system.

Once the convetred pdb file, netcdf file and parameters text file is available, run_hsa ad run_gist can be run as:  
```
$ run_hsa -i testcase.pdb -t md100ps.nc -l ligand.pdb -p params.txt -s 0 -f 100 -o testcase

run_gist -i testcase.pdb -t md100ps.nc -l ligand.pdb -p params.txt -g 20 20 20 -s 0 -f 100 -o testcase
```
This approach is applicable to any MD package that's not supported. SSTMap calculations are feasibale as long as you can convert your topology and trajectory files to suitable format and extract non-bonded paramters for every atom of your system into a text file with the same format as given above. 

## MD Trajectory Requirements for `SSTMap`
SSTMap calculations require an explicit solvent molecular dynamics trajectory with restrained solute. While `SSTMap` is  agnostic of the water model used, however, it has been tested only for `TIP3P`, `TIP4P`, `TIP4P-Ewald`, `TIP5P` and `OPC` models. This list will be updated as we test more water models. Also note that currently only simulation run in orthorhombic periodic boundary conditions are supported. We intend to provide support for non-orthorhombic periodic boundary conditions in future. The simulations can be generated in one of the following packages: [Amber](http://ambermd.org/), [Charmm](https://www.charmm.org), [Gromacs](http://www.gromacs.org/), [NAMD](http://www.ks.uiuc.edu/Research/namd/), [OpenMM](http://openmm.org/) and [Desmond](https://www.deshawresearch.com/resources_desmond.html). This list is based on simulation packages that are supported by [MDTraj](https://mdtraj.org) and [ParmEd](http://parmed.github.io/ParmEd/html/index.html), both of which are dependecnies of SSTMap. It's however, possible to expand the applicability of SSTMap calculations even beyond these packages, as long as the trajectory and toplogy formats can be converted to those currently supported (See the Desmond example for this apporach). 
 
## Command-line arguments

A list of command-line arguments for `run_hsa` and `run_gist` (with brief explanations) can be obtained at the terminal by running these commands without any argument or with `-h` flag. 
```
$ run_hsa -h
usage: run_hsa [-h] -i INPUT_TOP -t INPUT_TRAJ -l LIGAND [-f NUM_FRAMES]
               [-c CLUSTERS] [-p PARAM_FILE] [-s START_FRAME]
               [-d BULK_DENSITY] [-b CALC_HBONDS] [-o OUTPUT_PREFIX]

Run SSTMap site-based calculations (Hydration Site Analysis) through command-
line.

required arguments:
  -i INPUT_TOP, --input_top INPUT_TOP
                        Input toplogy File.
  -t INPUT_TRAJ, --input_traj INPUT_TRAJ
                        Input trajectory file.
  -l LIGAND, --ligand LIGAND
                        Input ligand PDB file.
  -f NUM_FRAMES, --num_frames NUM_FRAMES
                        Total number of frames to process.

optional arguments:
  -h, --help            show this help message and exit
  -c CLUSTERS, --clusters CLUSTERS
                        PDB file containing cluster centers.
  -p PARAM_FILE, --param_file PARAM_FILE
                        Additional parameter files, specific for MD package
  -s START_FRAME, --start_frame START_FRAME
                        Starting frame.
  -d BULK_DENSITY, --bulk_density BULK_DENSITY
                        Bulk density of the water model.
  -b CALC_HBONDS, --calc_hbonds CALC_HBONDS
                        True or False for whether to calculate h-bonds during
                        calculations.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for all the results files.
```

```
$ run_gist -h
usage: run_gist [-h] -i INPUT_TOP -t INPUT_TRAJ -l LIGAND -g GRID_DIM GRID_DIM
                GRID_DIM [-f NUM_FRAMES] [-p PARAM_FILE] [-s START_FRAME]
                [-d BULK_DENSITY] [-b CALC_HBONDS] [-o OUTPUT_PREFIX]

Run SSTMap grid-based (GIST) calculations through command-line.

required arguments:
  -i INPUT_TOP, --input_top INPUT_TOP
                        Input toplogy File.
  -t INPUT_TRAJ, --input_traj INPUT_TRAJ
                        Input trajectory file.
  -l LIGAND, --ligand LIGAND
                        Input ligand PDB file.
  -g GRID_DIM GRID_DIM GRID_DIM, --grid_dim GRID_DIM GRID_DIM GRID_DIM
                        grid dimensions e.g., 10 10 10
  -f NUM_FRAMES, --num_frames NUM_FRAMES
                        Total number of frames to process.

optional arguments:
  -h, --help            show this help message and exit
  -p PARAM_FILE, --param_file PARAM_FILE
                        Additional parameter files, specific for MD package
  -s START_FRAME, --start_frame START_FRAME
                        Starting frame.
  -d BULK_DENSITY, --bulk_density BULK_DENSITY
                        Bulk density of the water model.
  -b CALC_HBONDS, --calc_hbonds CALC_HBONDS
                        True or False for whether to calculate h-bonds during
                        calculations.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for all the results files.```

The arguments supplied for `-i`, `-t` and `-p` flags vary depending on the MD packages used for simulation. For demonstrative purposes, we use input topology and trajectories from a repository of test cases, which is available on [Github](https://github.com/KurtzmanLab/sstmap_test_suite). You can download the full test suite from [here](https://www.dropbox.com/sh/hrijgk8n5z12bgi/AABSigcBf9PN_7-Z26VCCPePa?dl=0) (since Github repository doesn't contain trajectory files). For a given platform, `cd` to its sub-directory and run the commands as shown below.
<!--more-->
