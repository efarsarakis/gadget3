from distutils.core import setup, Extension
import glob
import os
import shutil

incl_dirs = ['/usr/local/lib/python2.6/site-packages/numpy/core/include/numpy',
             '/usr/lib/python2.6/dist-packages/numpy/core/include/numpy',
             '/usr/local/include']
libs_dirs = ['/usr/local/lib']
libs = ['m']

calcGrid = Extension(   'calcGrid',
		    include_dirs = incl_dirs,
		    libraries    = libs,
		    library_dirs = libs_dirs,
                    sources = ['calcGrid.c'])

setup ( name = 'Python GADGET lib',
        version = '1.0',
        description = 'Scripts to work with GADGET input & output',
        author = 'Ruediger',
        author_email = 'rpakmor@mpa-garching.mpg.de',
        ext_modules = [calcGrid],
        script_args = ['build_ext', '--inplace'] )

for lib in glob.glob( "*.so" ):
    if os.path.exists( "../" + lib ):
        os.remove( "../" + lib )
    shutil.move( lib, "../" )


