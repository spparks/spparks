# plot_hcp_symmetry_histograms.py
# Author: John Mitchell, jamitch@sandia.gov, Sandia National Laboratories
# Date: March 8 2024
# Updated Sept 2025 John Mitchell and Svetoslav Nikolov

import matplotlib.pyplot as plt
import numpy
from numpy import average, std, max, linspace
# from sklearn.neighbors import KernelDensity

def get_axes():
    fig=plt.figure()
    ax=fig.subplots()
    return ax,fig

dbins={"disorientation":linspace(0,100,21)}

bandwidth={"disorientation":2.5}

dxy={"disorientation":[[0,.015],0.00125]}

def plot_distribution(d,kind='disorientation'):
    """
    N=1000000
    import cpp_quaternion as cppq
    import plot_hcp_symmetry_histograms as phcp

    Example(s):
    uqv=cppq.generate_random_unit_quaternions(N)
    div=cppq.estimate_disorientation(uqv,symmetry='HCP')
    di=cppq.to_numpy(div)
    phcp.plot_distribution(di,'disorientation')
    """
    symmetry='HCP'
    filename = "mtex-hcp-analytical"
    x1, y1 = numpy.loadtxt(filename, skiprows=0, usecols=(0,1), unpack=True)
    bins=dbins[kind]
    ax,fig=get_axes()
    ### Add curve for kernel density estimation
    ax.hist(d,bins=bins,density=True,label="SPPARKS "+symmetry+" "+kind)
    #n, binz = numpy.histogram(d,bins=bins,density=False)
    #ax.hist(mi,bins=bins,density=True,alpha=0.5,label='MTEX')
    data=d[:,numpy.newaxis]
    #db=bandwidth[kind]
    #kde=KernelDensity(kernel="gaussian",bandwidth=db).fit(data)
    # Create equally spaced points across orientation space from bins
    #domain=linspace(bins[0],bins[-1],50)[:,numpy.newaxis]
    #log_dens = kde.score_samples(domain)
    #ax.plot(domain[:,0],numpy.exp(log_dens),label="KDE, Gaussian")
    ax.plot(x1,y1,label='MTEX-Analytical')
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
    # Makes a tight bounding box
    plt.subplots_adjust(left=.135, right=.986, bottom=.094, top=.974, wspace=None, hspace=None)
    fig.show()
    fig.savefig(symmetry+"_"+kind+'.pdf')
    return data,bins
