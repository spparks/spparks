import argparse
import os, sys
import numpy
from numpy import min,max
import random
from matplotlib import pyplot as plt
from stitch.libstitch import libstitch

def main(args):
    prefix=args.prefix
    # Stitch filename
    stitch_filename=prefix+".st"
    # open stitch file
    rc,fid=libstitch.open(stitch_filename)

    # Command line arguments
    bx=args.block[0:2]
    by=args.block[2:4]
    bz=args.block[4:6]
    nx=bx[1]-bx[0]
    ny=by[1]-by[0]
    nz=bz[1]-bz[0]
    num_sites=nx*ny*nz
    bb=bx+by+bz
    state=numpy.ndarray(shape=(num_sites,),dtype=numpy.float64)
    state.fill(args.value)
    state=state.reshape((nx,ny,nz),order='F')
    block=numpy.fromiter(bb,dtype=numpy.int32).reshape(2,3,order='F')
    site=numpy.ndarray(shape=(num_sites,),dtype=numpy.int32)
    print("Randomizing initial site values with range=(1,%d)\n"%(num_sites,))
    random.seed()
    for i in range(num_sites):
        site[i]=random.randint(1,num_sites)
    site=site.reshape((nx,ny,nz),order='F')
    print("Randomizing site values min=%d, max=%d\n"%(min(site),max(site)))

    # 64 bit double field_type
    field_type=3
    num_vals=1
    default_val=-1.0
    (rc, float_id) = libstitch.create_field (fid, args.field, field_type, num_vals, default_val)

    # 32 bit integer field_type
    field_type=1
    num_vals=1
    default_val=-1
    (rc, site_id) = libstitch.create_field (fid, 'site', field_type, num_vals, default_val)

    # write 
    libstitch.write_block (fid, float_id, args.time, block, state);
    libstitch.write_block (fid, site_id,  args.time, block, site);

    # Close stitch file
    libstitch.close(fid)


# python init.py init --field='d1' 0 200 0 200 0 1
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="write stitch file with uniform value.")
    parser.add_argument('prefix',default="",type=str, \
                       help="stitch output file name prefix without .st suffix")
    s="dimensions of stitch file;\n"
    s+="block=xlo xhi ylo yhi zlo zhi;\n"
    parser.add_argument('block',nargs=6, type=int, help=s)
    parser.add_argument('--field',default='d1',type=str, \
                       help="stitch field created")
    parser.add_argument('--value',default=0.0,type=float, \
                       help="field value written to file")
    parser.add_argument('--time',default=0.0,type=float, \
                       help="field written to stitch file for this time")
    args=parser.parse_args()
    print(args)
    main(args)
