from setuptools import setup, Extension, find_packages
import numpy

__version__ = "1.1.4"

# define the extension module
extensions = []
extensions.append(Extension('_sstmap_ext',
                            sources=['sstmap/_sstmap_ext.c'],
                            include_dirs=[numpy.get_include()],
                            extra_link_args=['-lgsl','-lgslcblas']))

extensions.append(Extension('_sstmap_entropy',
                            sources=['sstmap/_sstmap_entropy.cpp', 'sstmap/kdhsa102.cpp'],
                            language="c++"))

extensions.append(Extension('_sstmap_probableconfig',
                            sources=['sstmap/_sstmap_probable.cpp', 'sstmap/probable.cpp'],
                            language="c++"))

setup(name='sstmap',
      author='Kamran Haider',
      author_email='kamranhaider.mb@gmail.com',
      description='SSTMap: A computational tool for studying structure and thermodynamics of water molecules on solute surfaces',
      version=__version__,
      license='MIT',
      url='https://github.com/KurtzmanLab/SSTMap',
      platforms=['Linux', 'Mac OS X',],
      packages=find_packages(),
      ext_modules=extensions,
      zip_safe=False,
      entry_points={
          'console_scripts':
              ['run_hsa = sstmap.scripts.run_hsa:entry_point',
               'run_gist = sstmap.scripts.run_gist:entry_point']}, )
