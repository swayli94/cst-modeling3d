
import numpy as np

from cst_modeling.basic import output_curves_igs, output_plot3d
from cst_modeling.section import cst_foil


if __name__ == "__main__":
    
    
    #* Build an airfoil
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    x0, yu, yl, _, _ = cst_foil(101, cst_u, cst_l, x=None, t=None, tail=0.0)
    
    x = [x0[None,:], x0[None,:]]
    y = [yu[None,:], yl[None,:]]
    z = [np.zeros_like(x[0]), np.zeros_like(x[0])]

    output_plot3d(x, y, z, fname='airfoil.grd')
    
    x = np.concatenate(x, axis=0)*1000
    y = np.concatenate(y, axis=0)*1000
    z = np.zeros_like(x)
    
    #output_curves_igs(x0, yu, np.zeros_like(x0), fname='airfoil-u.igs', n_degree=1, is_planar=True)
    #output_curves_igs(x0, yl, np.zeros_like(x0), fname='airfoil-l.igs', n_degree=0, is_planar=True)
    
    output_curves_igs(x, y, z, fname='airfoil.igs', n_degree=3, is_planar=True)
    