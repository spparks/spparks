#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Stitch library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-stitch args="-b"
Syntax from lib dir: python Install.py -b

specify option

  -b = build the Stitch library and create links to it

Example:

make lib-stitch args="-b"   # build in lib/stitch and set links
"""

# settings

version = "stitch-1.0"
url = "URL for Stitch library tarball when avaiable for download %s" % version

# print error message or help

def error(str=None):
  if not str: print(help)
  else: print("ERROR",str)
  sys.exit()

# expand to full path name
# process leading '~' or relative path

def fullpath(path):
  return os.path.abspath(os.path.expanduser(path))

def which(program):
  def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      path = path.strip('"')
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file

  return None

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

homepath = "."

buildflag = True
linkflag = True

version = "libstitch"
homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,version)

# build Stitch

if buildflag:
  print("Building Stitch ...")
  cmd = "cd %s; make" % homedir
  txt = subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  print(txt.decode('UTF-8'))

# create 2 links in lib/stitch to Stitch src dir

if linkflag:
  print("Creating links to Stitch include and lib files")
  if os.path.isfile("includelink") or os.path.islink("includelink"):
    os.remove("includelink")
  if os.path.isfile("liblink") or os.path.islink("liblink"):
    os.remove("liblink")
  cmd = 'ln -s "%s" includelink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
  cmd = 'ln -s "%s" liblink' % homedir
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
