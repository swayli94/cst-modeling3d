import os
import sys
sys.path.append('.')

import numpy as np
from matplotlib import pyplot as plt

from cst_modeling.section import cst_foil
from cst_modeling.basic import BasicSection, BasicSurface
from cst_modeling.io import plot3d_to_igs


if __name__ == "__main__":

    np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=200)
    path = os.path.dirname(sys.argv[0])


    #* Build an airfoil
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])


    x1, yu1, yl1, tmax1, rLE1 = cst_foil(201, cst_u, cst_l, x=None, t=None, tail=0.0)
    x2, yu2, yl2, tmax2, rLE2 = cst_foil(201, cst_u, cst_l, x=None, t=None, tail=0.004)
    x3, yu3, yl3, tmax3, rLE3 = cst_foil(201, cst_u, cst_l, x=None, t=0.11, tail=0.004)
    
    print('> maximum thickness:   %.3f  %.3f  %.3f'%(tmax1, tmax2, tmax3))
    print('> leading edge radius: %.3f  %.3f  %.3f'%(rLE1, rLE2, rLE3))
    print()
    
    
    plt.figure()
    plt.plot(x1, yu1, 'b', label='airfoil (tail=0.000)')
    plt.plot(x1, yl1, 'b', label='_nolegend_')
    
    plt.plot(x2, yu2, 'r--', lw=1, label='airfoil (tail=0.004)')
    plt.plot(x2, yl2, 'r--', lw=1, label='_nolegend_')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.07, 0.07))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title(r'Adding $t_\text{TE}$ without keeping $t_\text{max}$')
    plt.savefig(os.path.join(path, 'cst-airfoil-add-tail.png'), dpi=300)


    plt.figure()
    plt.plot(x1, yu1, 'b', label='airfoil (tail=0.000)')
    plt.plot(x1, yl1, 'b', label='_nolegend_')
    
    plt.plot(x3, yu3, 'r--', lw=1, label='airfoil (tail=0.004)')
    plt.plot(x3, yl3, 'r--', lw=1, label='_nolegend_')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.07, 0.07))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title(r'Adding $t_\text{TE}$ while keeping $t_\text{max}$')
    plt.savefig(os.path.join(path, 'cst-airfoil-add-tail-same-tmax.png'), dpi=300)


    #* 3D airfoil geometry
    airfoil = BasicSection(thick=None, chord=1, twist=0.0, lTwistAroundLE=True)
    
    airfoil.xx = x1
    airfoil.yu = yu1
    airfoil.yl = yl1

    geo3d = BasicSurface(n_sec=0, name='wing', nn=airfoil.xx.shape[0], ns=5)
    
    geo3d.secs = [airfoil]
    
    geo3d.geo()
    
    geo3d.output_tecplot(fname=os.path.join(path, 'wing.dat'))
    
    geo3d.output_plot3d(fname=os.path.join(path, 'wing.xyz'))
    
    #* Format transformation
    plot3d_to_igs(fname=os.path.join(path, 'wing'))

