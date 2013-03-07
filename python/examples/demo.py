#!/usr/local/bin/python -i
# preceeding line should have path for Python on your machine

# demo.py
# Purpose: illustrate use of many library interface commands
# Syntax:  demo.py
#          uses in.ising as SPPARKS input script

import sys

# parse command line

argv = sys.argv
if len(argv) != 1:
  print "Syntax: demo.py"
  sys.exit()

me = 0
# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()

from spparks import spparks
spk = spparks()

# test out various library functions after running in.demo

spk.file("in.ising")

if me == 0: print "\nPython output:"

nglobal = spk.extract("nglobal",0)
nlocal = spk.extract("nlocal",0)
print "Nglobal, nlocal =",nglobal,nlocal

xyz = spk.extract("xyz",5)
print "Y coord of 100th lattice site =",xyz[99][1]

eng = spk.energy()
print "Energy of system =",eng

site = spk.extract("site",1)
print "Site values for 1,10,100 =",site[0],site[9],site[99]

# uncomment if running in parallel via Pypar
#print "Proc %d out of %d procs has" % (me,nprocs), spk
#pypar.finalize()
