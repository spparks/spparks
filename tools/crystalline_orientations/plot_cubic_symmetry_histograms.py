# plot_cubic_symmetry_histograms.py
# Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
# Date: March 8 2024

import matplotlib.pyplot as plt
import numpy
from numpy import average, std, max, linspace, pi, cos, sin, sqrt

def get_axes():
    fig=plt.figure()
    ax=fig.subplots()
    return ax,fig

dbins={"min100":[0,10,20,30,40,50],
       "min110":[0,4,8,12,16,20,24,28,32],
       "min112":[0,2,4,6,8,10,12,14],
       "min123":[0,2,4,6,8,10],
       "disorientation":linspace(0,65,27)}

bandwidth={"min100":5, "min110":2, "min112":1, "min123":1, "disorientation":5}

dxy={"min100":[[0,.03],.0025],
     "min110":[[22,.0525],0.005],
     "min112":[[10,.10],0.00625],
     "min123":[[6,.2],0.025],
     "disorientation":[[0,.0275],0.0025]}

def last_interval(d):
    # TODO
    pass

def exact_distribution_handscomb(d):
    # if d<=45:
    #     return (2/15)*(1.0-cos(d*pi/180))
    # elif d<=60:
    #     return (2/15)*(3*(sqrt(2)-1)*sin(d)-2*(1-cos(d)))
    # elif d<60.6:
    #     return (2/15)*((3*(sqrt(2)-1)+4/sqrt(3))*sin(d)-6*(1-cos(d)))
    # Example usage
    # d = np.linspace(0, 60.5999, 100)
    # h = exact_distribution_handscomb(h)
    # Last interval is a mess! (d>=60.6) & (d<62.8)
    condlist = [d <= 45, (d>45) & (d<=60), (d>60) & (d<60.6)]
    flist=[
        lambda d: (2/15)*(1.0-cos(d*pi/180)),
        lambda d: (2/15)*(3*(sqrt(2)-1)*sin(d*pi/180)-2*(1-cos(d*pi/180))),
        lambda d: (2/15)*((3*(sqrt(2)-1)+4/sqrt(3))*sin(d*pi/180)-6*(1-cos(d*pi/180)))
            ]
    return numpy.piecewise(d, condlist, flist)

def plot_distribution(d,kind='min100'):
    """
    N=1000000
    import cpp_quaternion as cppq
    import plot_cubic_symmetry_histograms as pcsh

    Example(s):
    uqv=cppq.generate_random_unit_quaternions(N)
    min100v=cppq.estimate_min100(uqv)
    min100=cppq.to_numpy(min100v)
    min123v=cppq.estimate_min123(uqv)
    min123=cppq.to_numpy(min123v)
    min112v=cppq.estimate_min112(uqv)
    min112=cppq.to_numpy(min112v)
    div=cppq.estimate_disorientation(uqv)
    di=cppq.to_numpy(div)
    // This one is slow but it will finish
    min110v=cppq.estimate_min110(uqv)
    min110=cppq.to_numpy(min110v)
    pcsh.plot_distribution(min100,'min100')
    pcsh.plot_distribution(min110,'min110')
    pcsh.plot_distribution(min123,'min123')
    pcsh.plot_distribution(min112,'min112')
    pcsh.plot_distribution(di,'disorientation')
    """
    symmetry='CUBIC'
    bins=dbins[kind]
    ax,fig=get_axes()
    ### Add curve for kernel density estimation
    ax.hist(d,bins=bins,density=True,label=symmetry+" "+kind+' estimated')
    data=d[:,numpy.newaxis]
    ### Plot exact distribution (from Handscomb)
    x=linspace(0,60.5999,100)
    ex=exact_distribution_handscomb(x)
    ax.plot(x,ex,label="Exact (Handscomb)")
    ###
    # ax.hist(d,bins=bins,density=True,label=kind)
    ax.legend()
    ax.set_xlabel("Angle degrees")
    ax.set_ylabel("Probability Density")
    mean=average(d)
    sigma=std(d)
    mx=max(d)
    s1="Mean $\\mu={:4.1f}$".format(mean)
    s2="Std $\\sigma={:4.1f}$".format(sigma)
    s3="Max = ${:4.1f}$".format(mx)
    s4="N samples = ${:7d}$".format(len(d))
    xy=dxy[kind][0]
    delt=dxy[kind][1]
    ax.annotate(s1,xy)
    # adjust annotation y-component down
    y=xy[1]-delt
    ax.annotate(s2,[xy[0],y])
    y=y-delt
    ax.annotate(s3,[xy[0],y])
    y=y-delt
    ax.annotate(s4,[xy[0],y])
    fig.show()
    fig.savefig(symmetry+"_"+kind+'.pdf')
    return data,bins
