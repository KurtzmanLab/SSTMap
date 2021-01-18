from setuptools import setup, Extension, find_packages
import sys

if sys.version_info.major == 3 and sys.version_info.minor >= 8:
    sys.stderr.write("\nRuntimeError: Python 3.7.x or lower is required.\n"
                     "You are using Python {}.{}.{}."
                     "Please Install using Python version lower than Python 3.8.\n\n".format(sys.version_info.major,
                                                                                             sys.version_info.minor,
                                                                                             sys.version_info.micro))
    sys.exit(1)

try:
    import numpy
except ImportError as e:
    sys.stderr.write("\n{}: Numpy needed for installation. "
                     "Please install numpy <= 1.17.5.\n\n".format(type(e).__name__))
    sys.exit(1)

if '1.17' not in numpy.__version__ and int(numpy.__version__.split(".")[2]) > 17:
    sys.stderr.write("\nRuntimeError: Numpy <= 1.17.5 needed for installation."
                     "Please downgrade your numpy and install numpy <= 1.17.5.\n\n")
    sys.exit(1)


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
      python_requires='<3.8',
      install_requires=['parmed==3.2.0','matplotlib==2.2.3','mdtraj','numpy<=1.17.5'],
      setup_requires=['numpy<=1.17.5'],
      packages=find_packages(),
      ext_modules=extensions,
      zip_safe=False,
      entry_points={
          'console_scripts':
              ['run_hsa = sstmap.scripts.run_hsa:entry_point',
               'run_gist = sstmap.scripts.run_gist:entry_point']}, )
