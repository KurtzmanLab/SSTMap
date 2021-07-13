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
conda install -c omnia mdtraj numpy=1.17
```
##Step 5: We have tested gcc/g++=7,if you use Ubuntu20, you can use this [link]( https://linuxconfig.org/how-to-switch-between-multiple-gcc-and-g-compiler-versions-on-ubuntu-20-04-lts-focal-fossa) to swich to gcc/g++=7 compilerbefore go to next step. 

##Step 6: Change to the SSTMap directory where the setup.py locates, and run the installation command:

```bash
python setup.py install 
```
##Step 7: Now you can go to the [Ussage page](https://github.com/KurtzmanLab/SSTMap#usage) to run a quick test of your instalation. Of note, you have to run the commands within the sstmap_py36 virtual environment


If you have question regarding the installation of SSTMap, feel free to contact us by simpleliquid@gmail.com. If you want to run GIST, we recommand you to install CPPTRAJ(https://github.com/Amber-MD/cpptraj), where GIST tool is well maintained.
