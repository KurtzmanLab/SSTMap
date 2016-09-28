from setuptools import setup, Extension, find_packages
import numpy
# define the extension module
ext_module = Extension('_sstmap_ext', 
					sources=['sstmap/core/_sstmap_ext.c'], 
					include_dirs=[numpy.get_include()])
ext_module_old = Extension('_sstmap_ext_old', 
                              sources=['sstmap/core/_sstmap_ext_old.c'], 
                              include_dirs=[numpy.get_include()])
# run the setup
setup(name='sstmap',
      version='0.1',
      description='Library for analysis of water molecules in MD trajectories',
      url='https://github.com/KurtzmanLab/WatAnalTool',
      author='Kamran Haider',
      author_email='kamranhaider.mb@gmail.com',
      license='None',
      ext_modules=[ext_module, ext_module_old],
      zip_safe=False,
      packages=find_packages(),
      entry_points={'console_scripts':
          ['run_hsa = wateranalysistools.scripts.run_hsa:entry_point',
            'run_gist = wateranalysistools.scripts.run_phba:entry_point',]},
      )