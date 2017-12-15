---
layout: post
author: kamran-haider
title: Installing SSTMap
date: 2017-06-06
published: true
---
The main requirement for installation of `SSTMap` is the anaconda python distribution. If you are not familiar with Anaconda, please have a look at [continuum.io/why-anaconda](https://www.continuum.io/why-anaconda). SSTMap is implemented in python 2.7, so if you do not have the anaconda installation, downlaod the anaconda for python 2.7, as shown below. Alternatively, if you have a an anaconda installation for python 3.6, you can create a python 2.7 enviornment to run sstmap.

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
pip install git+git://github.com/kurtzmanlab/sstmap.git@v1.0#egg=sstmap
```

You can also download the release package manually from GitHub, unzip it, navigate into the directory, and execute the command:

```bash
https://github.com/KurtzmanLab/SSTMap/archive/v1.0.tar.gz or https://github.com/KurtzmanLab/SSTMap/archive/v1.0.zip
tar -xvf v1.0.tar.gz or unzip v1.0.zip
cd SSTMap-1.0
python setup.py install
```
Or you can clone the GitHub repository, navigate into the directory, and execute the command: 

```bash
git clone git@github.com:KurtzmanLab/SSTMap.git
cd SSTMap
```
**For the release version:**
```
git checkout tags/v1.0
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
