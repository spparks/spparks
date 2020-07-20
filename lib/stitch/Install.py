#!/usr/bin/env python

# Install.py tool to download, unpack, build, and link to the Stitch library
# used to automate the steps described in the README file in this dir

from __future__ import print_function
import sys,os,re,subprocess

# help message

help = """
Syntax from src dir: make lib-stitch args="-b"
                 or: make lib-stitch args="-p /usr/local/stitch-1.0"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/stitch-1.0

specify one or more options, order does not matter

  -b = download and build the Stitch library
  -p = specify folder of existing Stitch installation

Example:

make lib-stitch args="-b"   # download/build in lib/stitch/stitch-1.0
make lib-stitch args="-p $HOME/stitch-1.0" # use existing Stitch installation in $HOME/stitch-1.0
"""

# settings

version = "stitch-1.0"
url = "URL for Stitch library tarball when avaiable for download" % version

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

def geturl(url,fname):
  success = False

  if which('curl') != None:
    cmd = 'curl -L -o "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling curl failed with: %s" % e.output.decode('UTF-8'))

  if not success and which('wget') != None:
    cmd = 'wget -O "%s" %s' % (fname,url)
    try:
      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
      success = True
    except subprocess.CalledProcessError as e:
      print("Calling wget failed with: %s" % e.output.decode('UTF-8'))

  if not success:
    error("Failed to download source code with 'curl' or 'wget'")
  return

# parse args

args = sys.argv[1:]
nargs = len(args)
if nargs == 0: error()

homepath = "."
homedir = version

buildflag = False
pathflag = False
linkflag = True

iarg = 0
while iarg < nargs:
  if args[iarg] == "-v":
    if iarg+2 > nargs: error()
    version = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-p":
    if iarg+2 > nargs: error()
    stitchpath = fullpath(args[iarg+1])
    pathflag = True
    iarg += 2
  elif args[iarg] == "-b":
    buildflag = True
    iarg += 1
  else: error()

homepath = fullpath(homepath)
homedir = "%s/%s" % (homepath,version)

if (pathflag):
    if not os.path.isdir(stitchpath): error("Stitch path does not exist")
    homedir = stitchpath

if (buildflag and pathflag):
    error("Cannot use -b and -p flag at the same time")

if (not buildflag and not pathflag):
    error("Have to use either -b or -p flag")

# download and unpack Stitch tarball

#if buildflag:
#  print("Downloading Stitch ...")
#  geturl(url,"%s/%s.tar.gz" % (homepath,version))
#
#  print("Unpacking Stitch tarball ...")
#  if os.path.exists("%s/%s" % (homepath,version)):
#    cmd = 'rm -rf "%s/%s"' % (homepath,version)
#    subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
#  cmd = 'cd "%s"; tar -xzvf %s.tar.gz' % (homepath,version)
#  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
#  os.remove("%s/%s.tar.gz" % (homepath,version))
#  if os.path.basename(homedir) != version:
#    if os.path.exists(homedir):
#      cmd = 'rm -rf "%s"' % homedir
#      subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
#    os.rename("%s/%s" % (homepath,version),homedir)

# build Stitch

if buildflag:
  print("Building Stitch ...")
  cmd = "cd %s; make -f Makefile.steve stitch.lib" % homedir
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
