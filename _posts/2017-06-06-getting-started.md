---
layout: post
author: kamran-haider
title: Installing SSTMap
date: 2017-06-06
published: true
---
The main requirement for installation of `SSTMap` is the anaconda python distribution. If you are not familiar with Anaconda, please have a look at [continuum.io/why-anaconda](https://www.continuum.io/why-anaconda). SSTMap is implemented in python 2.7, so if you do not have the anaconda installation, downlaod the anaconda for python 2.7, as shown below. Altenratively, if you have a an anaconda installation for python 3.6, you can create a python 2.7 enviornment to run sstmap.

<!--more-->
For instructions on how to download and install Anaconda, please go to [Anaconda Download page](https://www.continuum.io/downloads). 
The installation is very easy and can be done with the following two steps. 

For Linux-64
```bash
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
bash Anaconda2-4.3.1-Linux-x86_64.sh
```
For MacOSX
```bash
curl -O https://repo.continuum.io/archive/Anaconda2-4.3.1-MacOSX-x86_64.sh
bash Anaconda2-4.3.1-MacOSX-x86_64.sh
```
If you already have an anaconda installation for python 3.6, you can install and run SSTMap in a seprate python 2.7 envioronment. The enviornment is created as follows:
```
conda create -n py27 python=2.7
source activate py27
```

### Conda Installation
The anaconda python distribution comes with `conda` package manager, which can then be used to install `SSTMap` with the following commands.

```
conda install sstmap -c omnia -c solvationtools
```

### GitHub Source Code
You can install the latest development version from the GitHub repository by executing

```
pip install git+git://github.com/kurtzmanlab/sstmap.git#egg=sstmap
```

You can also download the package manually from GitHub, unzip it, navigate into the package, and execute the command:

```
python setup.py install
```

When building from the source code, make sure that you manually install dependencies. You can do this by:
```
conda config --add channels omnia
conda install mdtraj parmed
``` 
