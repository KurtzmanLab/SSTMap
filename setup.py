from setuptools import setup, Extension, find_packages
import numpy
# define the extension module
ext_module = Extension('_sstmap_ext', 
					sources=['sstmap/_sstmap_ext.c'], 
					include_dirs=[numpy.get_include()])

setup(name='sstmap',
      version='0.1',
      description='Library for analysis of water molecules in MD trajectories',
      url='https://github.com/KurtzmanLab/SSTMap',
      author='Kamran Haider',
      author_email='kamranhaider.mb@gmail.com',
      license='None',
      ext_modules=[ext_module],
      zip_safe=False,
      packages=find_packages(),
      )
entry_points={
      'console_scripts':
      ['run_hsa = sstmap.scripts.run_hsa:entry_point',
      'run_gist = sstmap.scripts.run_gist:entry_point',
      'desmond_extract_nbparams = sstmap.scripts.desmond_extract_nbparams:entry_point',
      'dtr_to_netcdf = sstmap.scripts.dtr_to_netcdf:entry_point',]},

