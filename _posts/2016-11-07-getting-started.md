---
layout: post
title: Getting Started with SSTMap
---

## Installing `SSTMap`

The main requirement for installation of `SSTMap` is the anaconda python distribution. If you are not familiar with Anaconda, please have a look at: https://www.continuum.io/why-anaconda.

For instructions on how to download and install Anaconda, please go to: https://www.continuum.io/downloads. 
The installation is very easy and can be done with the following two steps (applicable to Linux 64-bit systems). 

```
wget http://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh
bash Anaconda2-4.0.0-Linux-x86_64.sh
```

The anaconda python distribution comes with `conda` package manager, which can then be used to install `SSTMap` and its dependencies. Once anaconda is installed on your system, run the following command to install SSTMap.

```
conda install sstmap -c kurtzman_lab
```