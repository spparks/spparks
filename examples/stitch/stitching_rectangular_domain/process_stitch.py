import getopt
import os, sys
import numpy
from matplotlib import pyplot as plt
from stitch.libstitch import libstitch

def main(prefix):
    # Stitch filename
    stitch_filename=prefix+".st"
    # open stitch file
    rc,fid=libstitch.open(stitch_filename)
    # Get all times stored in file
    rc,times=libstitch.get_times(fid)

    # Times in stitch file
    print("Times read=\n%s"%(times,))

    # Below are z-cuts; but x or y cuts may also be done
    bx=[49,50]
    bx=[0,100]
    by=[0,280]
    bz=[0,48] 
    bz=[0,1] # zcut bottom
    bz=[2,3] # zcut top of 1st layer
    bb=bx+by+bz
    block=numpy.fromiter(bb,dtype=numpy.int32).reshape(2,3,order='F')
    layer=0

    # Reads 'double value' stored by spparks as 'd1'
    rc,site_id=libstitch.query_field(fid,'site')
    if -1==site_id:
        raise RuntimeError("Field 'site' NOT FOUND")
    # rc,d1_id=libstitch.query_field(fid,'d1')
    #if -1==d1_id:
    #    raise RuntimeError("Field 'd1' NOT FOUND")

    for n,t in enumerate(times):
        # (rc, state, read_flag) = libstitch.read_block (fid, d1_id, t, block)
        (rc, state, read_flag) = libstitch.read_block (fid, site_id, t, block)
        print("time t=%7.1f, read_flag=%d, \n\tblock=\n%s"%(t,read_flag,block,))
        # NOTE 'zfill VERY necessary so that files are naturally order
        #    on linux/bash command line
        outfile=prefix+'_'+str(n).zfill(2)+".png"
        #print("output image file name = %s"%(outfile,))
        fig=plt.figure()
        ax=fig.subplots(1)
        cmap=plt.get_cmap("gist_rainbow")
        ## NEXT LINE is xcut
        #ax.imshow(state[layer,:,:].T,origin='lower',interpolation=None,cmap=cmap)
        ## NEXT LINE is ycut
        #ax.imshow(state[:,layer,:].T,origin='lower',interpolation=None,cmap=cmap)
        ## NEXT LINE is zcut
        ax.imshow(state[:,:,layer].T,origin='lower',interpolation=None,cmap=cmap)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        fig.savefig(outfile,bbox_inches='tight')
        plt.close(fig)

if __name__=='__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:",longopts=[])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)
    s=""
    prefix=None
    for o, a in opts:
        if o=="-p":
            prefix=a
            w="specified stitch file prefix = {0:s}\n".format(prefix)
            s+=w
        else:
            print("ERROR: invalid option(s)")
            usage()
            sys.exit(2)
    print(s)
    if None is prefix:
        raise ValueError("Invalid or missing stitch file prefix")
    main(prefix)
