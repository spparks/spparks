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
# John Mitchell (jamitch@sandia.gov) for more information.

import unittest
import numpy
from math import fabs,pi,cos,sin,log
from numpy import array, ndarray, dot, transpose

import cppyy
from cppyy.gbl.std import vector
import bor_variants as bor
import cpp_quaternion as cppq
quaternion=cppyy.gbl.SPPARKS_NS.quaternion

def copy_vector_to_numpy_array(v):
    n=v.size()
    q=numpy.ndarray(shape=(n,),dtype=numpy.float64)
    for i in range(n):
        q[i]=v[i]
    return q

class unit_hcp_bor_variants(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        # Compute randomly oriented beta prior
        one=1
        q0=quaternion.generate_random_unit_quaternions(one)
        # Compute variants associated with beta prior
        # q0=vector['double']((1,0,0,0))
        cpp_variants=bor.cpp_get_variants(q0=q0)
        python_variants=bor.get_variants(q0=q0)
        # Convert variant cpp vectors to numpy arrays for unit testing
        vc=12*[None]
        vp=12*[None]
        for i in range(12):
            vc[i]=copy_vector_to_numpy_array(cpp_variants[i])
            vp[i]=copy_vector_to_numpy_array(python_variants[i])
        self.vc=vc
        self.vp=vp
        self.cpp_variants=cpp_variants
        self.python_variants=python_variants

    def test_compare_python_variants_to_cpp_variants(self):
        # Sanity check to ensure both implementations of the variants are identical
        # Relative tolerance
        rtol = 1e-15
        # Absolute tolerance
        atol = 1e-15
        # Use numpy.testing.assert_allclose to compare the arrays
        for i in range(12):
            # print("vc =%s"%(self.vc[i]))
            # print("vp =%s"%(self.vp[i]))
            numpy.testing.assert_allclose(self.vc[i], self.vp[i], rtol=rtol, atol=atol)
        pass

    def test_variants(self):
        # Symmetries
        hcp=cppq.hcp.get_symmetries()
        # Variants stored as cpp vectors
        vc=self.cpp_variants
        vp=self.python_variants
        # Disorientations relative to the first variant
        refc=vc[0]
        refp=vp[0]
        # Expected disorientation angles based upon choosing first variant as reference
        expected=[0.0,10.53,90.00,90.00,60.83,63.26,60.00,60.83,60.83,60.00,63.26,60.83]
        di=numpy.fromiter(expected,dtype=numpy.float64)
        dic=numpy.ndarray(shape=(12,),dtype=numpy.float64)
        dip=numpy.ndarray(shape=(12,),dtype=numpy.float64)

        for i in range(12):
            dic[i]=cppq.disorientation.compute_disorientation(hcp,refc,vc[i])
            dip[i]=cppq.disorientation.compute_disorientation(hcp,refp,vp[i])
            # print("Disorientation cpp: i=%d, dic[%d]=%15.12f"%(i,i,dic[i]))
            # print("Disorientation python: i=%d, dic[%d]=%15.12f"%(i,i,dip[i]))

        # Sanity check to ensure both implementations of the variants are identical
        # Relative tolerance
        rtol = 1e-15
        # Absolute tolerance
        atol = 1e-15
        # Use numpy.testing.assert_allclose to compare the arrays
        numpy.testing.assert_allclose(dic, dip, rtol=rtol, atol=atol)
        numpy.testing.assert_allclose(di, dip, rtol=.01, atol=.01)
        numpy.testing.assert_allclose(di, dic, rtol=.01, atol=.01)
        pass

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
