---
layout: post
author: Eric Chen & Kamran-Haider
title: Installing SSTMap
date: 2019-12-06
published: true
---

### SSTMap Installation

##Step 1: Install Anaconda

The main requirement for installation of `SSTMap` is the anaconda python distribution.

<!--more-->


For instructions on how to download and install Anaconda, please go to [Anaconda Download page](https://www.anaconda.com/distribution/).
 

Alternatively, the Anaconda installation can be also done with the following two commands. 

For Linux-64
```bash
wget https://repo.continuum.io/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```
##Step 2: Create a virtual enviroment for SSTMap

<!--more-->

The virtual enviornment can be created as follows:

```
conda create -n sstmap_py36 python=3.6
source activate sstmap_py36
```


##Step 3: Get the gitHub Source Code by git


You can get the source code by git clone ( provided git has been installed on your machine):

```bash 
git clone https://github.com/KurtzmanLab/SSTMap.git

```

or you can also directly download the zip file by click the download botton on [the webpage](https://github.com/KurtzmanLab/SSTMap), then unzip it by:

```bash
unzip SSTMap-master.zip

```
##Step 4: Install the package required by SSTMap 

```bash
conda install -c omnia mdtraj
```
##Step 5: Change to the SSTMap directory where the setup.py locates, and run the installation command:

```bash
python setup.py install 
```
##Step 6: You have installed SSTMap successfully if you can run the following in your terminal (within the sstmap_py36 virtual environment)

```bash
run_hsa
```

If you have question regarding the installation of SSTMap, feel free to contact us by simpleliquid@gmail.com. If you want to run GIST, we recommand you to install CPPTRAJ(https://github.com/Amber-MD/cpptraj), where GIST tool is well maintained.
