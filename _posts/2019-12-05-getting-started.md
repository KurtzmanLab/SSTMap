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
##Step 4: Install the dependencies in case you do not have 

```bash
conda install -c omnia mdtraj
```
In Centos:
```bash
 sudo yum install gsl-devel
```
In Ubuntu:
```bash 
sudo apt-get install libgsl-dev  
ln path/to/libgsl.so path/to/libgsl.so.0
```

##Step 5: Install the compiler in case you do not have

In Centos:

```bash
sudo yum install g++
```

In Ubuntu:
```bash
sudo apt-get install g++
```
##Step 6: Install the SSTMap by wheel

if you get the souce code by git clone:
```bash
pip install ./SSTMap/devtools/wheel/sstmap-1.1.4-cp36-cp36m-linux_x86_64.whl
```

if you get the source code by clicking the download button:
```bash
pip install ./SSTmap-master/devtools/wheel/sstmap-1.1.4-cp36-cp36m-linux_x86_64.whl
```

##Step 7: You have installed SSTMap successfully if you can run these two commands in your terminal (within the sstmap_py37 virtual environment)

```bash
run_gist
```

or
```bash
run_hsa

```

If you have question regarding the installation of SSTMap, please contact us eric.clyang521@gmail.com, we will be glad to help!
