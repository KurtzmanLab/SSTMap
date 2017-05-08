---
layout: post
title: Getting Started with SSTMap
---
The main requirement for installation of `SSTMap` is the anaconda python distribution. If you are not familiar with Anaconda, please have a look at [continuum.io/why-anaconda](https://www.continuum.io/why-anaconda).

<!--more-->
For instructions on how to download and install Anaconda, please go to [Anaconda Download page](https://www.continuum.io/downloads). 
The installation is very easy and can be done with the following two steps (applicable to Linux 64-bit systems). 
```bash
wget https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
bash Anaconda2-4.3.1-Linux-x86_64.sh
```
For MacOSX
```bash
curl -O https://repo.continuum.io/archive/Anaconda2-4.3.1-MacOSX-x86_64.sh
bash Anaconda2-4.3.1-MacOSX-x86_64.sh
```

The anaconda python distribution comes with `conda` package manager, which can then be used to install `SSTMap` with the following commands.

```
conda config --add channels omnia
conda config --add channels solvationtools
conda install sstmap
```
