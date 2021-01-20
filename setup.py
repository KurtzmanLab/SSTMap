from setuptools import setup, Extension, find_packages
from os import environ, path
from subprocess import check_output
import sys
import shlex

# Check if we are using Python <= 3.7.x
if sys.version_info.major == 3 and sys.version_info.minor >= 8:
    sys.stderr.write("\nRuntimeError: Python 3.7.x or lower is required for SSTMap installation.\n"
                     "You are using Python {}.{}.{}.\n"
                     "Please Install a version of Python lower than Python 3.8.\n\n".format(sys.version_info.major,
                                                                                            sys.version_info.minor,
                                                                                            sys.version_info.micro))
    sys.exit(1)

# Check if Numpy is installed and if is Numpy <= 1.17.5
try:
    import numpy
except ImportError:
    sys.stderr.write("\nModuleNotFoundError: Numpy needed for installation. \n"
                     "Please install Numpy version <= 1.17.5.\n\n")
    sys.exit(1)

if '1.17' not in numpy.__version__ and int(numpy.__version__.split(".")[1]) > 17:
    sys.stderr.write("\nRuntimeError: Numpy version <= 1.17.5 needed for installation.\n"
                     "You are using Numpy version {}.\n"
                     "Please downgrade your numpy and install Numpy version <= 1.17.5.\n\n".format(numpy.__version__))
    sys.exit(1)

# Check C and C+ compiler version < 8.0
if environ.get('CC') is None:
    cc_ver_command_line = "gcc -dumpversion"
else:
    cc_ver_command_line = "{} -dumpversion".format(environ.get('CC'))

cc_ver_command = shlex.split(cc_ver_command_line)
cc_ver = None
try:
    if sys.version_info.major == 3:
        cc_ver = int(check_output(cc_ver_command, encoding='utf-8').strip()[0])
    else:
        cc_ver = int(check_output(cc_ver_command).strip()[0])
except:
    sys.stderr.write("\nCompilerNotFoundError: C and C++ compilers needed for installation.\n"
                     "Please install GNU C (gcc) and GNU C++ (g++) compiler version < 8.0.\n\n")

    sys.exit(1)

if cc_ver > 7:
    sys.stderr.write("\nRuntimeError: GCC and G++ compiler version < 8.0 needed for installation.\n"
                     "Please install GCC and G++  compile version < 8.0.\n"
                     "Set CC environment variable to the gcc executable (version < 8.0) and\n"
                     "set CXX environment variable to the g++ executable (version < 8.0).\n\n"
                     "Ex. If you install gcc and g++ version 7 then: \n"
                     "    export CC=/usr/bin/gcc-7\n"
                     "    export CXX=/usr/bin/g++-7\n\n")

    sys.exit(1)

# Check that the GNU Scientific Library developmental (gsl-dev) is installed
if not path.exists('/usr/include/gsl'):
    sys.stderr.write("\nRuntimeError: GNU Scientific Library (gls) needed for installation.\n"
                     "Please install the GNU Scientific Library development (gsl-dev)\n"
                     "using your Linux distribution package manager or install"
                     " from source (https://www.gnu.org/software/gsl/).\n\n")

    sys.exit(1)

__version__ = "1.1.4"

# define the extension module
extensions = []
extensions.append(Extension('_sstmap_ext',
                            sources=['sstmap/_sstmap_ext.c'],
                            include_dirs=[numpy.get_include()],
                            extra_link_args=['-lgsl', '-lgslcblas']))
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
      platforms=['Linux', 'Mac OS X', ],
      python_requires='<=3.7',
      install_requires=['parmed==3.2.0', 'matplotlib==2.2.3', 'mdtraj', 'numpy<=1.17.5'],
      setup_requires=['numpy<=1.17.5'],
      packages=find_packages(),
      ext_modules=extensions,
      zip_safe=False,
      entry_points={
          'console_scripts':
              ['run_hsa = sstmap.scripts.run_hsa:entry_point',
               'run_gist = sstmap.scripts.run_gist:entry_point']}, )
