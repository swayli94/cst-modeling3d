'''
Fit an airfoil (with raw data points) using CST method and output the geometry to other formats.
'''
import os
import sys
sys.path.append('.')

import numpy as np
from matplotlib import pyplot as plt

from cst_modeling.section import cst_foil, cst_foil_fit
from cst_modeling.io import output_curves_igs, output_plot3d


if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])
    
    #* Load airfoil raw data
    with open(os.path.join(path, 'naca1408.dat'), 'r') as f:
        lines = f.readlines()
    
        x = []
        y = []
        
        for line in lines[1:]:
            x_, y_ = line.split()
            x.append(float(x_))
            y.append(float(y_))
    
    xu0 = x[:18]; xu0.reverse(); xu0 = np.array(xu0)
    yu0 = y[:18]; yu0.reverse(); yu0 = np.array(yu0)
    xl0 = np.array(x[17:])
    yl0 = np.array(y[17:])
    

    #* Fit the airfoil
    cst_u, cst_l = cst_foil_fit(xu0, yu0, xl0, yl0, n_cst=10)

    xx, yu, yl, _, _ = cst_foil(101, cst_u, cst_l, x=None, t=None, tail=0.0)
    
    with open(os.path.join(path, 'cst-coefficients.txt'), 'w') as f:
        for i, c in enumerate(cst_u):
            f.write(f'cst_u-{i} {c}\n')
        for i, c in enumerate(cst_l):
            f.write(f'cst_l-{i} {c}\n')
    
    plt.figure()
    plt.plot(xu0, yu0, 'bo', label='Raw data')
    plt.plot(xl0, yl0, 'bo', label='_nolegend_')
    plt.plot(xx, yu, 'r', label='Fitting curve')
    plt.plot(xx, yl, 'r', label='_nolegend_')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'CST fitting of an airfoil')
    plt.savefig(os.path.join(path, 'cst-fitting-airfoil.png'), dpi=300)
    
    
    #* Output the geometry to other formats
    
    x = [xx[None,:], xx[None,:]]
    y = [yu[None,:], yl[None,:]]
    z = [np.zeros_like(x[0]), np.zeros_like(x[0])]

    output_plot3d(x, y, z, fname=os.path.join(path, 'airfoil-2d.grd'))
    
    x = np.concatenate(x, axis=0)*1000
    y = np.concatenate(y, axis=0)*1000
    z = np.zeros_like(x)
    
    output_curves_igs(x, y, z, fname=os.path.join(path, 'airfoil-2d.igs'), n_degree=3, is_planar=True)
    