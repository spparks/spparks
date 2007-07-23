# Script:  lattice.py
# Purpose: view screen dump of a SPPARKS lattice snapshot
# Syntax:  lattice.py dumpfile N nx ny (nz)
#          N = timestamp on snapshot
#          nx,ny,nz = lattice size (nz not used for 2d problem)

def compare(a,b):
  if a[0] < b[0]: return -1
  elif a[0] > b[0]: return 1
  else: return 0
  
if len(argv) != 5 and len(argv) != 6:
  raise StandardError, "Syntax: lattice.py dumpfile N nx ny (nz)"

dumpfile = argv[1]
nstep = int(argv[2])
nx = int(argv[3])
ny = int(argv[4])
nz = 0
if len(argv) == 6: nz = int(argv[5])

m = mdump(dumpfile)
isnap = m.findtime(nstep)
evalues = m.snaps[isnap].evalues

spins = []
for value in evalues:
  spins.append([int(value[0]),int(value[1])])

spins.sort(compare)

if nz == 0:
  if nx*ny != len(spins):
    raise StandardError, "Nx*Ny != # of spins"
  print "\n** 2d lattice for snapshot %d **" % nstep
  for j in range(ny-1,-1,-1):
    for i in range(nx):
      n = j*nx + i
      print "%2i" % spins[n][1],
    print

if nz:
  if nx*ny*nz != len(spins):
    raise StandardError, "Nx*Ny*Nz != # of spins"
  print "\n** 3d lattice for snapshot %d **" % nstep
  for k in range(nz):
    print "\n2d slice for k = %d" % (k+1)
    for j in range(ny-1,-1,-1):
      for i in range(nx):
        n = k*ny*nx + j*nx + i
        print "%2d" % spins[n][1],
      print
  
sys.exit()

  

