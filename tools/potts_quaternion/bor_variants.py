import argparse
import numpy
from numpy import sqrt, sin, cos, pi, arccos, arcsin
import cppyy
from cppyy.gbl.std import vector
cppyy.include("quaternion.h")
quaternion=cppyy.gbl.SPPARKS_NS.quaternion
cppyy.include("burgers_orientation_relation_bcc_to_hcp.h")
bor=cppyy.gbl.SPPARKS_NS.BurgersOrienationRelation_BCC_to_HCP

# Allowable variants ['V1', 'V2', ..., 'V12']
VALID_VARIANTS = [f'V{i}' for i in range(1, 13)]  

def cpp_get_variants(q0=vector['double']((1,0,0,0))):
    """
    Calls c++ header 'burgers_orientation_relation_bcc_to_hcp' to get variants.
    """
    v=bor.get_variants(q0);
    return v

def cpp_get_variants_dictionary(q0):
    """
    Calls c++ header 'burgers_orientation_relation_bcc_to_hcp' to get variants.
    Adds a dictionary for use with 'orient_hcp.py' plotting utility.
    """
    v=bor.get_variants(q0);
    variants={}
    for i,c in enumerate(VALID_VARIANTS):
        variants[c]=v[i]
    return variants;

## HACK to tweak variants calculation and still plot
def get_variants_dictionary(q0):
    """
    Calls python function below to get variants.
    Adds a dictionary for use with 'orient_hcp.py' plotting utility.
    """
    variants=get_variants(q0)
    variants_dic={}
    for i,v in enumerate(variants):
        variants_dic[VALID_VARIANTS[i]]=v
    return variants_dic

# Utility lambda functions for calculating quaternions
rot_x=lambda x: vector['double']((cos(x*pi/180/2),sin(x*pi/180/2),0,0))
rot_y=lambda x: vector['double']((cos(x*pi/180/2),0,sin(x*pi/180/2),0))
rot_z=lambda x: vector['double']((cos(x*pi/180/2),0,0,sin(x*pi/180/2)))

def get_variant(q0,q1,q2,q3,q4):
    """
    Utility function for computing active variants.  Called by python 
    function below 'get_variants'
    """
    r1=quaternion.composition(q0,q1)
    r2=quaternion.composition(r1,q2)
    r3=quaternion.composition(r2,q3)
    r4=quaternion.composition(r3,q4)
    return r4

def get_variants(q0=vector['double']((1,0,0,0))):
    """

    Returns so-called active orientations.

    SPPARKS stores the active quaternion orientations on
    each lattice site.

    Active or passive orientations can be used to compute disorientations
    between variants -- provided that the disorientation calculation
    is tweaked to take into account 'active' or 'passive' variants.

    See compute_disorientation in disorientation.h

    For plotting, the output of this function
    will orient/rotate the variant and depicts
    the variant in its actual orientation.

    """
    # HEXAGONAL prism
    # Orientation of hexagonal prism following bor variants
    sqrt2=sqrt(2)
    sqrt3=sqrt(3)
    phi=180*arcsin(1.0/sqrt3)/pi

    # Variants are ordered following 
    # Farabi 2020, Table 2
    variants=12*[None]

    # xy edge
    qxy=rot_x(0)
    # xz edge
    qxz=rot_x(90)
    # yz edge
    qyz=quaternion.composition(rot_z(90), rot_x(90))

    #
    # Table 2 Karthikeyan (110)
    variants[0]=get_variant(q0,qxy,rot_z(-45),rot_x(90),rot_z(phi))
    variants[1]=get_variant(q0,qxy,rot_z(-45),rot_x(90),rot_z(60-phi))
    #
    # Table 2 Karthikeyan (1b10) or (11b0)
    variants[2]=get_variant(q0,qxy,rot_z(45),rot_x(90),rot_z(60-phi))
    variants[3]=get_variant(q0,qxy,rot_z(45),rot_x(90),rot_z(phi))
    #
    # Table 2 Karthikeyan (101)
    variants[4]=get_variant(q0,qxz,rot_z(-45),rot_x(90),rot_z(phi))
    variants[5]=get_variant(q0,qxz,rot_z(-45),rot_x(90),rot_z(60-phi))
    #
    # Table 2 Karthikeyan (1b01)
    variants[6]=get_variant(q0,qxz,rot_z(45),rot_x(90),rot_z(phi))
    variants[7]=get_variant(q0,qxz,rot_z(45),rot_x(90),rot_z(60-phi))
    #
    # Table 2 Karthikeyan (011)
    variants[8]=get_variant(q0,qyz,rot_z(-45),rot_x(90),rot_z(phi))
    variants[9]=get_variant(q0,qyz,rot_z(-45),rot_x(90),rot_z(60-phi))
    #
    # Table 2 Karthikeyan (01b1)
    variants[10]=get_variant(q0,qyz,rot_z(45),rot_x(90),rot_z(phi))
    variants[11]=get_variant(q0,qyz,rot_z(45),rot_x(90),rot_z(60-phi))
    return variants
