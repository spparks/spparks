# cpp_quaternion.py
# Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
# Date: March 8 2024

import cppyy
import numpy
from cppyy.gbl.std import vector

cppyy.include("quaternion.h")
cppyy.include("disorientation.h")
cppyy.include("cubic_symmetries.h")
cppyy.include("hcp_symmetries.h")

def to_numpy(vector,reshape=None):
    """
    Input 'vector' must be std::vector<double>
    Returns numpy array view of input 'vector'

    if reshape==None: returns array same shape/size
    if type(reshape)==tuple -- then returns reshaped array
       per specified reshape=(nrows,mcols)

    NOTE on numpy views of std::vector
    If in a python function a temporary std::vector is created and from that vector
    a numpy array is created and returned from said function, then the numpy array
    becomes invalid outside the function because the temporary local std::vector
    that numpy array was created from goes out of scope and its underlying memory
    is deleted -- this is the memory that the numpy array is pointing to -- hence
    the returned numpy array is invalid outside the function.  For this reason,
    if a python wrapper function calls a cpp function, and the return value from
    the cpp function is a std::vector, then its best to return the std::vector directly
    from the python function and then call 'to_numpy' if its desirable to also
    have a numpy array view of the std::vector.  Then the numpy array will remain
    valid as long as the std::vector is in scope.
    """
    n=vector.size()
    # Following cppyy documentation
    d=vector.data()
    d.reshape((n,))
    # Convert to numpy array
    nv=numpy.frombuffer(d,dtype=numpy.float64,count=n)
    # Reshape 
    if None is not reshape and type(reshape) is tuple:
        c=reshape[1]
        if 0!=n%c:
            raise ValueError("reshape specified invalid; require 0=n%reshape[1] where n is len(vector)")
        nv=nv.reshape(reshape)
    return nv

def get_identity_quaternion():
    uvw=vector['double']((1.0,0.0,0.0,0.0))
    return uvw

def generate_random_unit_quaternions(n):
    """
    returns uqv 
      where uqv is Lowlevel c++ vector<double> len=n*4
    """
    # Returns vector<double> uqv(n*4)
    uqv=cppyy.gbl.SPPARKS_NS.quaternion.generate_random_unit_quaternions(n)
    return uqv

def estimate_disorientation(uqv,symmetry='CUBIC'):
    """
    Inputs
    uqv: type Lowlevel c++ vector<double>, len=n*4 where n
      is number of quaternions as returned from above
      'generate_random_unit_quaternions
    symmetry: type 'string' label of symmetry type -- expecting
      'CUBIC' or 'HCP'
    """
    s=None
    di=None
    if 'CUBIC'==symmetry:
        s=cppyy.gbl.SPPARKS_NS.CUBIC.get_symmetries()
    elif 'HCP'==symmetry:
        s=cppyy.gbl.SPPARKS_NS.HCP.get_symmetries()
    else:
        print("Warning: did not match symmetry of 'CUBIC' or 'HCP'; returning 'None'\n")
        return None

    di=cppyy.gbl.SPPARKS_NS.disorientation.estimate_disorientation_distribution(uqv,s)
    return di

def estimate_min100(uqv):
    """
    Input uqv 
      where uqv is Lowlevel c++ vector<double> len=n*4
      -- as returned from above 'generate_random_unit_quaternions
    """
    min100=cppyy.gbl.SPPARKS_NS.CUBIC.estimate_min100(uqv)
    return min100

def estimate_min110(uqv):
    """
    Input uqv
      where uqv is Lowlevel c++ vector<double> len=n*4
      -- as returned from above 'generate_random_unit_quaternions
    """
    min110=cppyy.gbl.SPPARKS_NS.CUBIC.estimate_min110(uqv)
    return min110


def estimate_min123(uqv):
    """
    Input uqv 
      where uqv is Lowlevel c++ vector<double> len=n*4
      -- as returned from above 'generate_random_unit_quaternions
    """
    uvw=vector['double']((1.0,2.0,3.0))
    min123=cppyy.gbl.SPPARKS_NS.CUBIC.estimate_minuvw(uvw,uqv)
    return min123

def estimate_min112(uqv):
    """
    Input uqv 
      where uqv is Lowlevel c++ vector<double> len=n*4
      -- as returned from above 'generate_random_unit_quaternions
    """
    uvw=vector['double']((1.0,1.0,2.0))
    min112=cppyy.gbl.SPPARKS_NS.CUBIC.estimate_minuvw(uvw,uqv)
    return min112



