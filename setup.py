from setuptools import setup, Extension, find_packages
import numpy
import subprocess
import os
# define the extension module
extensions = []
extensions.append(Extension('_sstmap_ext', 
					sources=['sstmap/_sstmap_ext.c'], 
					include_dirs=[numpy.get_include()]))

extensions.append(Extension('_sstmap_bruteclust', 
                              sources=['scripts/make_clust_brute.cpp'], 
                              language='c++'))


setup(name='sstmap',
      version='0.1',
      description='Library for analysis of water molecules in MD trajectories',
      url='https://github.com/KurtzmanLab/SSTMap',
      author='Kamran Haider',
      author_email='kamranhaider.mb@gmail.com',
      license='None',
      ext_modules=extensions,
      zip_safe=False,
      packages=find_packages(),
      entry_points={
      'console_scripts':
      ['run_hsa = sstmap.scripts.run_hsa:entry_point',
      'run_gist = sstmap.scripts.run_gist:entry_point',
      'desmond_extract_nbparams = sstmap.scripts.desmond_extract_nbparams:entry_point',
      'dtr_to_netcdf = sstmap.scripts.dtr_to_netcdf:entry_point',]},
      )

#g++ -o /usr/local/bin/kdhsa102 $RECIPE_DIR/../../scripts/kdhsa102.cpp $RECIPE_DIR/../../scripts/kdhsa102_main.cpp
#g++ -o /usr/local/bin/6dimprobable $RECIPE_DIR/../../scripts/6dimprobable.cpp $RECIPE_DIR/../../scripts/6dim_main.cpp
#g++ -o bruteclust make_clust_brute.cpp
c_prog = os.path.abspath("scripts/make_clust_brute.cpp")
try:
      subprocess.check_call("g++  -o /usr/local/bin/bruteclust " + c_prog)
except Exception as e:
      print(e)
print("g++  -o /usr/local/bin/bruteclust " + c_prog)