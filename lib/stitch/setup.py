from setuptools import Extension, setup
from os import path
import numpy

libstitch_ext=Extension(
                        name="stitch.libstitch", 
                        include_dirs=[numpy.get_include()],
                        sources=["src/stitch/libstitch/sqlite3.c",
                                 "src/stitch/libstitch/stitch.c",
                                 "src/stitch/libstitch/stitchmodule.c"]
                        )

setup(ext_modules=[libstitch_ext])
