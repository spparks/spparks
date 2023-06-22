import argparse
import os, sys
import numpy
from numpy import min,max
from matplotlib import pyplot as plt
from stitch.libstitch import libstitch
#import libstitch

def main(args):
    prefix=args.prefix
    # Stitch filename
    stitch_filename=prefix+".st"
    # open stitch file
    rc,fid=libstitch.open(stitch_filename)
    # Get all times stored in file
    rc,times=libstitch.get_times(fid)

    # Command line arguments
    cut_axis=args.cut_axis
    locus=args.locus
    w=args.window

    if 'x'==cut_axis:
        bx=[locus,locus+1]
        by=[w[0],w[1]]
        bz=[w[2],w[3]]
    elif 'y'==cut_axis:
        bx=[w[0],w[1]]
        by=[locus,locus+1]
        bz=[w[2],w[3]]
    elif 'z'==cut_axis:
        bx=[w[0],w[1]]
        by=[w[2],w[3]]
        bz=[locus,locus+1]
    else:
        raise ValueError("Invalid 'cut_axis'")

    # Times in stitch file
    #print("Times read=\n%s"%(times,))

    bb=bx+by+bz
    block=numpy.fromiter(bb,dtype=numpy.int32).reshape(2,3,order='F')
    layer=0

    field=args.field
    rc,site_id=libstitch.query_field(fid,field)
    if -1==site_id:
        raise RuntimeError("Field "+ field+" NOT FOUND")

    # Read step specified
    t=times[args.time_step]
    n=args.time_step
    if -1==args.time_step:
        n=len(times)-1
    print("Number of time steps in file = %d\n"%(len(times),))
    print("\tfirst time value t=%7.1f; last time VALUE t=%7.1f \n"%(times[0],t,))
    (rc, state, read_flag) = libstitch.read_block (fid, site_id, t, block)
    print("\t@t=%7.1f, field '%s' min value = %d, max value = %d\n"%(t,field,min(state),max(state),))

    # Close stitch file
    libstitch.close(fid)

    # 'state' is a numpy array
    print("numpy array state.flags=%s"%(state.flags,))
    print("numpy array state.shape=%s"%(state.shape,))
    print("numpy array state.dtype=%s"%(state.dtype,))

    outfile=prefix+'_'+cut_axis+"_"+str(locus).zfill(4)+"_"+"step_"+str(n).zfill(2)+".png"
    print("output image file name = %s"%(outfile,))
    fig=plt.figure()
    ax=fig.subplots(1)
    cmap=plt.get_cmap("gist_rainbow")
    ## NEXT LINE is xcut
    if 'x'==cut_axis:
        #ax.imshow(state[layer,:,:].T,origin='lower',interpolation=None,cmap=cmap)
        ax.imshow(state[layer,:,:].T,origin='lower',interpolation='none',cmap=cmap)
    ## NEXT LINE is ycut
    elif 'y'==cut_axis:
        #ax.imshow(state[:,layer,:].T,origin='lower',interpolation=None,cmap=cmap)
        ax.imshow(state[:,layer,:].T,origin='lower',interpolation='none',cmap=cmap)
    ## NEXT LINE is zcut
    elif 'z'==cut_axis:
        #ax.imshow(state[:,:,layer].T,origin='lower',interpolation=None,cmap=cmap)
        ax.imshow(state[:,:,layer].T,origin='lower',interpolation='none',cmap=cmap)
    else:
        raise ValueError("Invalid 'cut_axis'")
    #ax.axhline(y=70,xmin=0,xmax=628,color='k')
    #ax.axhline(y=70+64,xmin=0,xmax=628,color='k')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.savefig(outfile,bbox_inches='tight')
    plt.close(fig)

# python plot_stitch_cut.py stitch_rectangle_zcut --field=site 0 6000 0 1250
# python plot_stitch_cut.py equiaxed --time_step=20 --field=site --cut_axis='x' --locus=100 0 200 0 200
# python plot_stitch_cut.py equiaxed --time_step=0 --field=d1 --cut_axis='z' --locus=0 0 500 0 500
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="create 2d images from stitch file.")
    parser.add_argument('prefix',default="",type=str, \
                       help="stitch file name prefix without .st suffix")
    s="dimensions of cut/image plane;\n"
    s+="must make sense relative to input stitch file;\n"
    s+="if --cut_axis='z'; then window=xlo xhi ylo yhi;\n"
    parser.add_argument('window',nargs=4,
                       help=s)
    parser.add_argument('--cut_axis',default="z",type=str, choices=['x','y','z'], \
                       help="define axis cut plane ")
    parser.add_argument('--locus',default=0,type=int, \
                       help="define locus of cut plane; must be an integer into stitch file ")
    parser.add_argument('--field',default='i1',type=str, \
                       help="stitch field to plot")
    parser.add_argument('--time_step',default='-1',type=int, \
                       help="integer time step # in file to plot")
    
    args=parser.parse_args()
    print(args)
    main(args)
