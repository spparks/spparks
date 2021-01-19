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
    # TODO 
    # Figure out how to configure one of these to bettwer document installation of libstitch
    #config.add_npy_pkg_config('libstitch/libstitch.ini.in','stitch/libstitch',subst_dict={'foo':'stitch'})
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

      
