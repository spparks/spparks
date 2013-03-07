# ----------------------------------------------------------------------
#   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
#   http://www.cs.sandia.gov/~sjplimp/spparks.html
#   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
#   Copyright (2008) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level SPPARKS directory.
# -------------------------------------------------------------------------

# Python wrapper on SPPARKS library via ctypes

import sys,traceback,types
from ctypes import *

class spparks:
  def __init__(self,name="",cmdargs=None):

    # load libspparks.so by default
    # if name = "g++", load libspparks_g++.so
    
    try:
      if not name: self.lib = CDLL("libspparks.so",RTLD_GLOBAL)
      else: self.lib = CDLL("libspparks_%s.so" % name,RTLD_GLOBAL)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise OSError,"Could not load SPPARKS dynamic library"

    # create an instance of SPPARKS
    # don't know how to pass an MPI communicator from PyPar
    # no_mpi call lets SPPARKS use MPI_COMM_WORLD
    # cargs = array of C strings from args
    
    if cmdargs:
      cmdargs.insert(0,"spparks.py")
      narg = len(cmdargs)
      cargs = (c_char_p*narg)(*cmdargs)
      self.spk = c_void_p()
      self.lib.spparks_open_no_mpi(narg,cargs,byref(self.spk))
    else:
      self.spk = c_void_p()
      self.lib.spparks_open_no_mpi(0,None,byref(self.spk))
      # could use just this if SPPARKS lib interface supported it
      # self.spk = self.lib.spparks_open_no_mpi(0,None)

  def __del__(self):
    if self.spk: self.lib.spparks_close(self.spk)

  def close(self):
    self.lib.spparks_close(self.spk)
    self.spk = None

  def file(self,file):
    self.lib.spparks_file(self.spk,file)

  def command(self,cmd):
    self.lib.spparks_command(self.spk,cmd)

  def extract(self,name,type):
    if type == 0:
      self.lib.spparks_extract.restype = POINTER(c_int)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr[0]
    if type == 1:
      self.lib.spparks_extract.restype = POINTER(c_int)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 2:
      self.lib.spparks_extract.restype = POINTER(POINTER(c_int))
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 3:
      self.lib.spparks_extract.restype = POINTER(c_double)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr[0]
    if type == 4:
      self.lib.spparks_extract.restype = POINTER(c_double)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 5:
      self.lib.spparks_extract.restype = POINTER(POINTER(c_double))
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    return None

  def energy(self):
    self.lib.spparks_energy.restype = c_double
    return self.lib.spparks_energy(self.spk)
