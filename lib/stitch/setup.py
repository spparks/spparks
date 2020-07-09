import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    from distutils.sysconfig import get_python_inc
    lapack = get_info('lapack_opt')
    config=Configuration('stitch',parent_package,top_path)
    config.add_include_dirs(get_python_inc())
    config.add_subpackage('libstitch')

    # Installed libstitch
    sources = ['sqlite3.c','stitch.c']
    config.add_installed_library('stitch',
                                 sources=[join('libstitch',x) for x in sources],
                                 install_dir='stitch/lib',
                                 build_info=dict(extra_compiler_args=['-std=gnu99']))
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
