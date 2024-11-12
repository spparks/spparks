import argparse
import json
import os
import shutil
from math import pow
import numpy
from numpy import sqrt, fabs
import matplotlib.pyplot as plt

def compute_coefficients(control_points):
    """
    Function which makes for simple and fast evaluation of bezier 
    curve and its evaluation as a power series.
    
    Allows bezier polynomial to be evaluated as a power series.  Coefficients 
    computed below correspond with 4th order bezier polynomial with control 
    points.
        
    input: array len 5 of control points for 4th order bezier polynomial.
    output: array len 5 of coefficients (c0,c1,c2,c3,c4) as defined below.
    double value = c0 + c1 * u + c2 * pow(u,2) + c3 * pow(u,3) + c4 * pow(u,4);
    """
    # alias to for simplicity
    x=control_points
    x0 = x[0];
    x1 = x[1];
    x2 = x[2];
    x3 = x[3];
    x4 = x[4];

    c=numpy.zeros(shape=(5,),dtype=numpy.float64)
    c[0] = x0;
    c[1] = (-4 * x0 + 4 * x1);
    c[2] = (6 * x0 - 12 * x1 + 6 * x2);
    c[3] = (-4 * x0 + 12 * x1 - 12 * x2 + 4 * x3);
    c[4] = (x0 - 4 * x1 + 6 * x2 - 4 * x3 + x4);

    return c

def get_4th_order_curve(xpts,ypts):
    """
    Compute power series polynomial.

    Input: array len=5, double values for x and y control points

    Output: tuple, len=2, power series operators for xpts and ypts

    Power series curve
    p numpy.polynomial.Polynomial
    p = c0 + c1 * u + c2 * u2     + c3 * u3 +        c4 * u4;
    """
    cx=compute_coefficients(xpts)    
    cy=compute_coefficients(ypts)    
    x=numpy.polynomial.Polynomial(cx,domain=(0,1),window=(0,1))
    y=numpy.polynomial.Polynomial(cy,domain=(0,1),window=(0,1))
    return x,y

def get_normal_to_curve(g):
    x=g[0]
    y=g[1]
    dx=x.deriv(1)
    dy=y.deriv(1)
    nx=dy
    ny=-dx
    return nx,ny

def compute_umax_at_maximum_width(gy):
    """
    input: y-component operator for top bezier curve
    """
    dgy=gy.deriv(1)
    ddgy=gy.deriv(2)
    u0=0.5
    tolerance=1.0e-9
    max_iter=10
    u=u0
    f=dgy(u)
    j=ddgy(u)
    resid=fabs(f)
    iter=0
    while(resid>tolerance and iter<=max_iter):
        # Compute increment and update newton iteration
        # u+du = inv(Jacobian) * residual
        u += -f / j;
        # Evaluate residual and Jacobian at new iterate
        f=dgy(u)
        j=ddgy(u)
        resid = fabs(f);
        iter += 1;

    return u


# HOWTO RUN this script.
# Input json file: variable_haz_demo.json
# Run this script WITHOUT '.json' file extension
# python variable_haz_plot.py vhaz
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot utility for variable haz.")
    parser.add_argument('prefix',default=None,type=str, \
                       help="'json' filename prefix (no 'json' file extension).")
    args=parser.parse_args()
    prefix=args.prefix
    with open(prefix+".json",'r') as f:
        params=json.load(f)
    f.close()
    cp=params['control_points']
    haz=params['variable_haz']
    box_limits=params['box_limits']
    xpts=cp['X']
    ypts=cp['Y']

    # Pool length
    pool_length=xpts[-1]-xpts[0]

    # Insert '0' at start and end of y control_points.
    # SPPARKS implicitly assumes pool is symmetric about 
    #   centerline -- first and last y-component
    #   of control points are on centerline and hence y=0
    ypts.insert(0,0)
    ypts.append(0)
    g=get_4th_order_curve(xpts,ypts)

    # Parametric coordinate array [0.0,1.0]
    u = numpy.linspace(0.0, 1.0, 1000)

    # Evaluate curve
    x=g[0](u)
    y=g[1](u)
    # Lower curve
    plt.plot(x,y,color='red',label='melt pool')
    # Upper curve
    plt.plot(x,-y,color='red')

    # Parametric coordinate where pool is maximum width
    umax=compute_umax_at_maximum_width(g[1])

    # Compute pool width
    pool_width=2*fabs(g[1](umax))

    # Index location of value that is strictly less than or equal to umax
    # imax=numpy.argmax(u>umax)-1
    # print("umax=%10.8f\n"%(umax,))
    # vals=(imax,u[imax-1],u[imax],u[imax+1])
    # print("imax=%d, u[imax-1]=%10.8f, u[imax]=%10.8f, u[imax+1]=%10.8f\n"%vals)

    # Create normal to curve, normalize normal, and compute variable haz
    h0=haz['h0']
    ht=haz['ht']
    n=haz['n']
    h1=(ht-h0)/pow(umax,n);
    gn=get_normal_to_curve(g)
    gnx=gn[0](u)
    gny=gn[1](u)
    hx=numpy.zeros(shape=len(u),dtype=numpy.float64)
    hy=numpy.zeros(shape=len(u),dtype=numpy.float64)
    for i in range(len(u)):
        nx=gnx[i]
        ny=gny[i]
        norm=sqrt(nx*nx+ny*ny)
        vx=nx/norm
        vy=ny/norm
        delt=u[i]-umax
        if delt<0:
            # Region of variable haz
            delt=fabs(u[i]-umax)
            hx[i]=x[i]+h0*vx+h1*vx*pow(delt,n)
            hy[i]=y[i]+h0*vy+h1*vy*pow(delt,n)
        else:
            # Foward of max pool width -- fixed haz if used at all
            hx[i]=x[i]+h0*vx
            hy[i]=y[i]+h0*vy

    ## Plot variable haz
    # Lower curve
    plt.plot(hx,hy,color='blue',label='haz',linewidth=2.5)
    # Upper curve
    plt.plot(hx,-hy,color='blue',linewidth=2.5)

    # Add label to plot for input parameters
    label="pool: width={w:5.1f}, length={l:5.1f}\n".format(w=pool_width,l=pool_length)
    label+="$h_0={h0}$, $h_t={ht}$, $n={n:4.2f}$".format(h0=h0,ht=ht,n=n)
    # Text position may require adjustment that depends upon control points
    # This is a good start, but vertical positioning may need adjustment
    px=xpts[0]+pool_length/10
    py=-10
    plt.text(px,py, label, fontsize=16)
    plt.axis('equal')
    plt.axis('off')
    xlim=box_limits['xlim']
    ylim=box_limits['ylim']
    plt.xlim(xlim[0],xlim[1])
    plt.ylim(ylim[0],ylim[1])
    plt.tight_layout()
    image_name=params['image_name']
    plt.savefig(image_name,bbox_inches='tight', pad_inches=0,transparent=True)
