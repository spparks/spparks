#!/usr/bin/env python3
# encoding: utf-8
# Copyright 2019 National Technology & Engineering Solutions of
# Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525
# with NTESS, the U.S. Government retains certain rights in this software.
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# 
# For more information, contact Jay Lofstead (gflofst@sandeia.gov) or
# John Mitchell (jamitch@sandia.gov) for more information.

from os.path import join

def get_libraries_list(lapack_info):
    try:
        libraries=list(lapack_info['libraries'])
    except KeyError:
        libraries=list()
    return libraries

def get_library_dirs_list(lapack_info):
    try:
        library_dirs=list(lapack_info['library_dirs'])
    except KeyError:
        library_dirs=list()
    return library_dirs

def get_include_dirs_list(lapack_info):
    try:
        include_dirs=list(lapack_info['include_dirs'])
    except KeyError:
        include_dirs=list()
    return include_dirs

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    lapack = get_info('lapack_opt')
    blas = get_info('blas_opt')
    config=Configuration('libstitch',parent_package,top_path)

    # Install stitch.h relative to 'install-headers' defined
    #    either in 'setup.cfg' or in command line;
    # python setup.py install --install-headers='path to base'
    config.add_headers(['include/',('stitch.h')])

    # MAC and Linux have different 'lapack' dictionaries that
    #   can cause the build to fail if not properly handled;
    extra_compile_args=[]
    extra_link_args=[]
    import platform
    if 'Darwin'==platform.system():
        try:
            extra_compile_args=lapack['extra_compile_args']
        except KeyError:
            # not compile args
            pass
        try:
            extra_link_args=lapack['extra_link_args']
        except KeyError:
            # not link args
            pass

    # package
    libraries=get_libraries_list(lapack);
    libraries.append('stitch')
    library_dirs=get_library_dirs_list(lapack);
    sources=['stitchmodule.c','sqlite3.c','stitch.c']
    include_dirs=['.']
    extra_compile_args.append('-std=gnu99')
    config.add_extension('libstitch',
                         sources=[join('',x) for x in sources],
                         #define_macros=[("STITCH_N_TO_1",'')],
                         include_dirs=include_dirs,
                         libraries=libraries,
                         library_dirs=library_dirs,
                         extra_compile_args=extra_compile_args,
                         extra_link_args=extra_link_args)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
