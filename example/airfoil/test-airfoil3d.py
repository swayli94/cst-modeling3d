'''
This is an example of generating 2D geometry in 3D form for ICEM CFD and CFL3D
'''

import numpy as np
from cst_modeling.foil import cst_foil
from cst_modeling.basic import BasicSection, BasicSurface, plot3d_to_igs


if __name__ == "__main__":
    
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    airfoil = BasicSection(thick=None, chord=1, twist=0.0, lTwistAroundLE=True)
    
    airfoil.xx, airfoil.yu, airfoil.yl, t0, R0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)


    geo3d = BasicSurface(n_sec=0, name='Wing', nn=airfoil.xx.shape[0], ns=11)
    
    geo3d.secs = [airfoil]
    
    geo3d.geo()
    
    geo3d.output_tecplot(fname='Wing.dat')
    
    geo3d.output_plot3d(fname='Wing.grd')
    
    plot3d_to_igs(fname='Wing')


