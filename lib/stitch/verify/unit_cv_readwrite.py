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

import unittest
import numpy
from math import fabs,pi,cos,sin,log
from numpy import array, ndarray, dot, transpose
from stitch.libstitch import libstitch

def pretty_string_block(block):
    """
    Use this to print string representation of block.

    input bb: type=array(shape=(2,3), dtype=numpy.int32)
              represents bounding box for input to spparks for use
              in the 'region' command.

    Following strings are concatenated for output.
    output: string = 'Number of sites' %d
    output: string='x0,x1,y0,y1,z0,z1 = %d, %d, %d, %d, %d, %d\n'
    """
    x_str="x0,x1,y0,y1,z0,z1={:d},{:d},{:d},{:d},{:d},{:d}\n"
    t=get_block_tuple(block)
    block_str=x_str.format(*t)
    return block_str

def get_block_tuple(block):
    """
    input block: type numpy.ndarray, shape=(2,3), order='Fortran'
    output x0,x1,y0,y1,z0,z1: type=tuple

    Example
    -------
    input:
        box=numpy.fromiter([0,127, 86,344,0,30],dtype=numpy.int32).reshape(2,3,order='F')
    output:
        (0,127,86,344,0,30)
    """
    shape=block.shape
    if 2!=shape[0] or 3!=shape[1]:
        raise ValueError("Input 'block' must be ndarray shape=2,3")
    if not numpy.dtype(numpy.int32) is block.dtype:
        raise TypeError("Input block must be ndarray dtype=numpy.32")

    t=(block[0,0],block[1,0],block[0,1],block[1,1],block[0,2],block[1,2])
    return t

def get_block_width(block):
    """
    input block: type numpy.ndarray, shape=(2,3), order='Fortran'
    output nx,ny,nz: type tuple, length=3, values=nx,ny,nz

    Example input:
        block=numpy.fromiter([0,127, 86,344,0,30],dtype=numpy.int32).reshape(2,3,order='F')
    """
    shape=block.shape
    if 2!=shape[0] or 3!=shape[1]:
        raise ValueError("Input 'block' must be ndarray shape=2,3")
    if not numpy.dtype(numpy.int32) is block.dtype:
        raise TypeError("Input block must be ndarray dtype=numpy.32")

    nx=block[1,0]-block[0,0];
    ny=block[1,1]-block[0,1];
    nz=block[1,2]-block[0,2];
    return nx,ny,nz

def get_Q(box):
    nx,ny,nz=get_block_width(box)
    return nx*ny*nz


class unit_fixed_size_cv_write_read(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_fixed_size_cv_write_read'
        # 'suffix' 'st' for 'stitch file'
        self._fname=name+'.st'

        self._state_dtype=numpy.int32
        # x1, x2, y1, y2, z1, z2
        b1=numpy.fromiter([0,127, 86,344,0,30],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b1=b1
        nx,ny,nz=get_block_width(b1)
        self._nx=nx;
        self._ny=ny;
        self._nz=nz;
        self._Q=nx*ny*nz
        (rc,self._file) = libstitch.open (self._fname);
        (rc,self._field_id) = libstitch.create_field (self._file, 'spin', 1, 1, -1)
        self._absolute_tolerance = 1.0e-9;
        self._relative_tolerance = 1.0e-15;
        self._no_value_present = -1;
        rc = libstitch.set_parameters (self._file, self._absolute_tolerance, self._relative_tolerance, self._no_value_present);
        pass

    def write(self, time, block, state):
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, new_time) = libstitch.write_block (self._file, field_id, time, block, state)
        return new_time

    def read(self, time, block):
        # Read data
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, state, new_time) = libstitch.read_block (self._file, field_id, time, block)
        return (state, new_time)

    def test_write_read(self):
        self.assertTrue(self._Q==982980)
        nx,ny,nz=get_block_width(self._b1)
        state=numpy.ndarray(shape=(nx,ny,nz),dtype=self._state_dtype,order='F')
        t0=0.0
        self._t0 = t0
        # Assign state
        for k in range(nz):
            # Assign values according to 'z' elevation
            v=k;
            state[:,:,k]=v

        b1_str=pretty_string_block(self._b1)
        self.write(self._t0, self._b1, state)

        # read
        (trial_state, new_time)=self.read(self._t0, self._b1)
        # print ('t0 ', self._t0)
        # print ('trial_state')
        # print (trial_state)
        # print ('state')
        # print (state)
        # print ('get the parameters')
        # (rc, absolute_tolerance, relative_tolerance, first_t, last_t) = libstitch.get_parameters (self._file);
        # print ('absolute_tolerance ', absolute_tolerance, ' relative_tolerance ', relative_tolerance, ' first_time ', first_t, ' last_time ', last_t);

        # Assert that trial was correctly read
        self.assertTrue((trial_state == state).all())
        pass

    def tearDown(self):
        # close the file
        libstitch.close (self._file);
        pass

class unit_variable_size_cv_write_read(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        name='unit_variable_size_cv_write_read'

        # 'suffix' 'st' for 'stitch file'
        self._fname=name+'.st'
        (rc,self._file) = libstitch.open (self._fname);
        (rc,field_id) = libstitch.create_field (self._file, 'spin', 1, 1, -1)
        self._absolute_tolerance = 1.0e-9;
        self._relative_tolerance = 1.0e-15;
        self._no_value_present = -1;
        rc = libstitch.set_parameters (self._file, self._absolute_tolerance, self._relative_tolerance, self._no_value_present);

        self._state_dtype=numpy.int32
        # Define 3 distinct blocks and then a union block
        # x1, x2, y1, y2, z1, z2
        # Blocks are 'non-overlapping'
        self._b1=numpy.fromiter([0,2,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b2=numpy.fromiter([2,4,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        self._b3=numpy.fromiter([4,6,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        # This is a union of all 3 blocks b1 b2 b3
        self._bu=numpy.fromiter([0,6,0,2,0,2],dtype=self._state_dtype).reshape(2,3,order='F')
        # Define state which will be written to each block
        self._s1=  numpy.ones(shape=get_block_width(self._b1),dtype=self._state_dtype,order='F')
        self._s2=2*numpy.ones(shape=get_block_width(self._b2),dtype=self._state_dtype,order='F')
        self._s3=3*numpy.ones(shape=get_block_width(self._b3),dtype=self._state_dtype,order='F')

        self._state_dtype_i32=numpy.int32
        self._state_dtype_i64=numpy.int64
        self._state_dtype_f64=numpy.float64
        # multi-field test area
        self._b4=numpy.fromiter([6,8,0,2,0,2],dtype=self._state_dtype_i32).reshape(2,3,order='F')
        # multi-field test blocks
        self._s4=4*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_i32,order='F')
        self._s5=5*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_i64,order='F')
        self._s6=6.0*numpy.ones(shape=get_block_width(self._b4),dtype=self._state_dtype_f64,order='F')

        # retrieve or define a set of times; 
        # (rc,times) = libstitch.get_times (self._file) // don't do
        self._t0=0.0
        self._t1=1.0
        self._t2=2.0
        self._t3=3.0
        pass

    def tearDown(self):
        # close the file
        libstitch.close (self._file);
        pass

    # NOTE: for these tests to work correctly, they must be run in order of occurrence 
    #   as they are in this file; otherwise the latter test could fail.
    # One way to apparently achieve this is to add alphabetical orders to 
    #   at beginning of each test name.
    def test_a_write_read_b1_at_t1(self):
        # Write/read b1 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, new_time) = libstitch.write_block (self._file, field_id, self._t1, self._b1, self._s1);
        # print ('write time ', self._t1)
        (rc, trial_s1, new_time) = libstitch.read_block (self._file, field_id, self._t1, self._b1)
        self.assertTrue((trial_s1 == self._s1).all())
        # add test for the bounding box once that is implemented
        pass

    def test_b_read_flag_tests(self):
        # read the t1/b1 block and check the new_time flag
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_s2, new_time) = libstitch.read_block (self._file, field_id, self._t1, self._b1)
        #print (new_time)
        self.assertTrue (new_time == 0) # a true would be a non-zero value
        pass

    def test_c_read_on_block_b1_at_nonexistent_prior_time(self):
        # This time 't0' is prior to time 't0' previously written;
        #    So read should get all '-1' values
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_m1, new_time) = libstitch.read_block (self._file, field_id, self._t0, self._b1)
        m1=-1*numpy.ones(shape=get_block_width(self._b1),dtype=self._state_dtype,order='F')
        # print ('t0 ', self._t0)
        # print ('b1 ', self._b1)
        # print ('m1')
        # print (m1)
        # print ('trial_m1')
        # print (trial_m1)
        self.assertTrue((trial_m1 == m1).all())
        pass

    def test_d_write_read_b2_at_t2(self):
        # Write/read b2 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, new_time) = libstitch.write_block (self._file, field_id, self._t2, self._b2, self._s2);
        (rc, trial_s2, new_time) = libstitch.read_block (self._file, field_id, self._t2, self._b2)
        self.assertTrue((trial_s2 == self._s2).all())
        pass

    def test_e_write_read_b3_at_t3(self):
        # Write/read b3 and 'assert'
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, new_time) = libstitch.write_block (self._file, field_id, self._t3, self._b3, self._s3);
        (rc, trial_s3, new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b3)
        self.assertTrue((trial_s3 == self._s3).all())
        pass

    def test_f_read_write_on_union(self):
        # Create su on block union 'bu'
        # It is already correct for s1 @ t3
        # Assign for s2 @ t3
        # Assign for s3 @ t3
        su=numpy.ones(shape=get_block_width(self._bu),dtype=self._state_dtype,order='F')
        su[2:4,:,:]=2
        su[4:6,:,:]=3
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_su, new_time) = libstitch.read_block (self._file, field_id, self._t3, self._bu)
        # FAILS
        # print ('self._t3 ', self._t3)
        # print ('bb')
        # print (self._bu)
        # print ('su')
        # print (su)
        # print ('trial_su')
        # print (trial_su)
        self.assertTrue((trial_su == su).all())
        # Now write union at later time
        t4=4.0
        su[:,:,:]=4
        (rc, new_time) = libstitch.write_block (self._file, field_id, t4, self._bu, su);
        # Assert 
        self.assertIsNotNone(self._file)
        self.assertIsNotNone(self._bu)
        self.assertIsNotNone(t4)
        self.assertIsNotNone(libstitch.read_block)
        ret = libstitch.read_block (self._file, field_id, t4, self._bu)
        # (rc, trial_s4, new_time) = libstitch.read_block (self._file, field_id, t4, self._bu)
        self.assertIsNotNone (ret)
        self.assertIsNotNone (ret [0])
        self.assertIsNotNone (ret [1])
        self.assertIsNotNone (ret [2])
        rc = ret [0]
        trial_s4 = ret [1]
        new_time = ret [2]
        self.assertTrue((trial_s4 == su).all())
        pass

    def test_g_read_blocks_at_later_time(self):
        # re-run existing test
        self.test_c_read_on_block_b1_at_nonexistent_prior_time()

        # Re-read original blocks at t3 and assert values
        # b1
        # 1 = STITCH_TYPE_INT32. Need to use an enum
        (rc, field_id) = libstitch.query_field (self._file, 'spin')
        (rc, trial_s1, new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b1)
        self.assertTrue((trial_s1 == self._s1).all())
        # b2
        (rc, trial_s2, new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b2)
        self.assertTrue((trial_s2 == self._s2).all())
        # b3
        (rc, trial_s3, new_time) = libstitch.read_block (self._file, field_id, self._t3, self._b3)
        self.assertTrue((trial_s3 == self._s3).all())
        pass

    def test_h_multifield_write_read (self):
        # create three fields (one of each type), write each, read each
        #print ('create 3 different type fields')
        (rc, field_id_i32) = libstitch.create_field (self._file, 'field_int32', 1, 1, -1)
        (rc, field_id_i64) = libstitch.create_field (self._file, 'field_int64', 2, 1, -3000000000)
        (rc, field_id_f64) = libstitch.create_field (self._file, 'field_float64', 3, 1, -1.0)

        (rc, field_ids, labels, t, lengths, no_value_presents) = libstitch.get_fields (self._file)
        #for i in range (len (field_ids)):
        #    print (field_ids [i], ' ', labels [i], ' ', t [i], ' ', lengths [i], ' ', no_value_presents [i])

        # setup for getting the block where data is new?
        libstitch.write_block (self._file, field_id_i32, self._t1, self._b4, self._s4)
        libstitch.write_block (self._file, field_id_i64, self._t1, self._b4, self._s5)
        libstitch.write_block (self._file, field_id_f64, self._t1, self._b4, self._s6)

        (rc, trial_su1_i32, is_new_time) = libstitch.read_block (self._file, field_id_i32, self._t1, self._b4)
        (rc, trial_su1_i64, is_new_time) = libstitch.read_block (self._file, field_id_i64, self._t1, self._b4)
        (rc, trial_su1_f64, is_new_time) = libstitch.read_block (self._file, field_id_f64, self._t1, self._b4)

        self.assertTrue ((trial_su1_i32 == self._s4).all ())
        self.assertTrue ((trial_su1_i64 == self._s5).all ())
        self.assertTrue ((trial_su1_f64 == self._s6).all ())
        pass

if __name__ == "__main__":
    unittest.main()
