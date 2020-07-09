#!/usr/bin/env python3
# encoding: utf-8
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
            extra_compile_args=None
        try:
            extra_link_args=lapack['extra_link_args']
        except KeyError:
            extra_link_args=None

    # package
    libraries=get_libraries_list(lapack);
    libraries.append('stitch')
    library_dirs=get_library_dirs_list(lapack);
    sources=['stitchmodule.c','sqlite3.c','stitch.c']
    include_dirs=['.']
    extra_compile_args.append('-std=gnu99')
    config.add_extension('libstitch',
                         sources=[join('',x) for x in sources],
                         include_dirs=include_dirs,
                         libraries=libraries,
                         library_dirs=library_dirs,
                         extra_compile_args=extra_compile_args,
                         extra_link_args=extra_link_args)

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
