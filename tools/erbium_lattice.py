#!/usr/local/bin/python

# Python script to create a 3-fold lattice for Erbium model
# creates a lattice file that SPPARKS can read in as part of app_style command

# Syntax: erbium_lattice.py Nx Ny Nz latfile
#   Nx,Ny,Nz = cubic unit cells in each dimension
#   latfile = SPPARKS lattice file to create

# Lattice geometry = 3 conjoined lattices
# F = fcc site
# O = octahedral site
# T = tetrahedral site
#
# F sites are on fcc lattice
#   4 sites per cubic unit cell
# O sites are edge centers of cubic unit cell + 1 in center (bcc site)
#   4 sites per cubic unit cell
# T sites are diamond sites 1/4 towards center from each of 8 corners
#   8 sites per cubic unit cell
#   form a simple cubic lattice with spacing 1/2 of fcc cubic unit cell
#
# neighbors of each F site: 12 F, 6 O, 8 T
# neighbors of each O site: 6 F, 12 O, 8 T
# neighbors of each T site: 4 F, 4 O, 6 T
#
# for SPPARKS to assign F,O,T to site IDs, this script IDs them as follows:
# within each unit cell, 4 F sites are listed first, then 4 O, then 8 T

import sys
from math import sqrt

if len(sys.argv) != 5:
  raise StandardError, "Syntax: erbium_lattice.py Nx Ny Nz latfile"

nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])
latfile = sys.argv[4]

FCC = 0
OCTA = 1
TETRA = 2

nbasis = 4 + 4 + 8
maxconn = 12 + 6 + 8

# nearest neighbor distances for all combinations of lattice sites

epsilon = 0.01
ffdist = sqrt(2)/2 + epsilon
fodist = 0.5 + epsilon
ftdist = sqrt(3)/4 + epsilon
oodist = sqrt(2)/2 + epsilon
otdist = sqrt(3)/4 + epsilon
ttdist = 0.5 + epsilon

# generate vertices, one unit cell at a time

xyz = []

for k in xrange(nz):
  for j in xrange(ny):
    for i in xrange(nx):
      xyz.append([i+0.0,j+0.0,k+0.0])       # 4 FCC sites
      xyz.append([i+0.5,j+0.5,k+0.0]) 
      xyz.append([i+0.5,j+0.0,k+0.5]) 
      xyz.append([i+0.0,j+0.5,k+0.5]) 

      xyz.append([i+0.5,j+0.0,k+0.0])       # 4 OCTA sites
      xyz.append([i+0.0,j+0.5,k+0.0])
      xyz.append([i+0.0,j+0.0,k+0.5])
      xyz.append([i+0.5,j+0.5,k+0.5])
      
      xyz.append([i+0.25,j+0.25,k+0.25])    # 8 TETRA sites
      xyz.append([i+0.75,j+0.25,k+0.25])
      xyz.append([i+0.25,j+0.75,k+0.25])
      xyz.append([i+0.75,j+0.75,k+0.25])
      xyz.append([i+0.25,j+0.25,k+0.75])
      xyz.append([i+0.75,j+0.25,k+0.75])
      xyz.append([i+0.25,j+0.75,k+0.75])
      xyz.append([i+0.75,j+0.75,k+0.75])


# create 3x3x3 array of unit cells of atoms
# for each atom, store x,y,z,type,idelta,jdelta,kdelta,basis
# type = FCC,OCTA,TETRA
# idelta,jdelta,kdelta = -1,0,+1 for which unit cell atom is in
# basis = 0 to nbasis-1 for which basis atom it is

x = 27*nbasis*[0]
y = 27*nbasis*[0]
z = 27*nbasis*[0]
type = 27*nbasis*[0]
idelta = 27*nbasis*[0]
jdelta = 27*nbasis*[0]
kdelta = 27*nbasis*[0]
basis = 27*nbasis*[0]

n = 0
for kk in [-1,0,1]:
  for jj in [-1,0,1]:
    for ii in [-1,0,1]:
      for m in xrange(nbasis):
        x[n] = xyz[m][0] + ii
        y[n] = xyz[m][1] + jj
        z[n] = xyz[m][2] + kk
        if m < 4: type[n] = FCC
        elif m < 8: type[n] = OCTA
        else: type[n] = TETRA
        idelta[n] = ii
        jdelta[n] = jj
        kdelta[n] = kk
        basis[n] = m
        n += 1

# create cmap connectivity array
# tells how to find each neighbor of each basis atom
# dimensions of cmap = nbasis by nneigh by 4
# [0],[1],[2] = offset of neighbor in I,J,K cubic unit cells
# [3] = which basis atom the neighbor is, in new unit cell

# double loop over atoms in central cell by loop over all atoms
# find classes of nearest F,O,T neighbors and append them to cmap array

cmap = []
nstart = 13*nbasis
m = 0

for i in xrange(nstart,nstart+nbasis,1):
  cmap.append([])
  
  for j in xrange(27*nbasis):
    if i == j: continue

    dx = x[i] - x[j]
    dy = y[i] - y[j]
    dz = z[i] - z[j]
    r = sqrt(dx*dx + dy*dy + dz*dz)

    itype = type[i]
    jtype = type[j]
    if itype == FCC:
      if jtype == FCC and r < ffdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == OCTA and r < fodist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == TETRA and r < ftdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
    elif itype == OCTA:
      if jtype == FCC and r < fodist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == OCTA and r < oodist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == TETRA and r < otdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
    elif itype == TETRA:
      if jtype == FCC and r < ftdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == OCTA and r < otdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
      elif jtype == TETRA and r < ttdist:
        cmap[m].append([idelta[j],jdelta[j],kdelta[j],basis[j]])
  m += 1

# sanity check on cmap counts of each kind of neighbor

for m in xrange(nbasis):
  if m < 4 and len(cmap[m]) != 26:
    print "ERROR: Bad cmap array FCC length",m,len(cmap[m])
    sys.exit(1)
  if m >= 4 and m < 8 and len(cmap[m]) != 26:
    print "ERROR: Bad cmap array OCTA length",m,len(cmap[m])
    sys.exit(1)
  if m >= 8 and len(cmap[m]) != 14:
    print "ERROR: Bad cmap array TETRA length",m,len(cmap[m])
    sys.exit(1)

  fcount = ocount = tcount = 0
  for i in xrange(len(cmap[m])):
    if cmap[m][i][3] < 4: fcount += 1
    elif cmap[m][i][3] < 8: ocount += 1
    else: tcount += 1
  if m < 4:
    if fcount != 12 or ocount != 6 or tcount != 8:
      print "ERROR: Bad cmap array FCC counts"
      sys.exit(1)
  elif m < 8:
    if fcount != 6 or ocount != 12 or tcount != 8:
      print "ERROR: Bad cmap array OCTA counts"
      sys.exit(1)
  else:
    if fcount != 4 or ocount != 4 or tcount != 6:
      print "ERROR: Bad cmap array TETRA counts"
      sys.exit(1)
    
# generate edges

edges = []

for k in xrange(nz):
  for j in xrange(ny):
    for i in xrange(nx):
      for m in xrange(nbasis):
        neighs = []

        for n in xrange(len(cmap[m])):      # neighbors per site
          ii = i + cmap[m][n][0]            # ii,jj,kk = new cubic unit cell
          jj = j + cmap[m][n][1]
          kk = k + cmap[m][n][2]
          mm = cmap[m][n][3]                # new basis atom in new cell
          if ii < 0: ii += nx               # apply PBC
          if ii >= nx: ii -= nx
          if jj < 0: jj += ny
          if jj >= ny: jj -= ny
          if kk < 0: kk += nz
          if kk >= nz: kk -= nz
          id = kk*nx*ny*nbasis + jj*nx*nbasis + ii*nbasis + mm + 1
          neighs.append(id)
          
        edges.append(neighs)
      
# write SPPARKS lattice file

fp = open(latfile,"w")
print >>fp,\
      "Erbium 3-fold lattice for %d by %d by %d cubic unit cells" % (nx,ny,nz)
print >>fp
print >>fp,3,"dimension"
print >>fp,len(xyz),"vertices"
print >>fp,maxconn,"max connectivity"
print >>fp,0,nx,"xlo xhi"
print >>fp,0,ny,"ylo yhi"
print >>fp,0,nz,"zlo zhi"

print >>fp
print >>fp,"Vertices"
print >>fp

for i in xrange(len(xyz)):
  print >>fp,i+1,xyz[i][0],xyz[i][1],xyz[i][2]
  
print >>fp
print >>fp,"Edges"
print >>fp

for i in xrange(len(xyz)):
  print >>fp,i+1,
  for j in xrange(len(edges[i])): print >>fp,edges[i][j],
  print >>fp

fp.close()
