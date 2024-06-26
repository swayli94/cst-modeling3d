
import numpy as np
from cst_modeling.basic import output_foil
from cst_modeling.section import cst_foil
from cst_modeling.tools.xfoil import read_xfoil_polar, run_xfoil, read_xfoil_dump, foil_for_XFoil
import os
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=200)
    
    
    #* Build an airfoil
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    x, yu, yl, t0, R0 = cst_foil(201, cst_u, cst_l, x=None, t=0.10, tail=0.0)
    
    foil_for_XFoil(x, yu, yl, fname='airfoil.dat')
    output_foil(x, yu, yl, fname='foil.dat')
    
    #* Run XFoil
    run_xfoil(AoAs=[-2.0, 5.0, 1.0], Cls=None,
    Minf=0.1, Re=1e6, nNodes=161, iterVis=10, 
    fname_airfoil='airfoil.dat', delete_temp=True,
    fname_cp=None, fname_raw='raw-aoas.bin', fname_polar='polar-aoas.dat')

    result1 = read_xfoil_dump('raw-aoas.bin')
    result2 = read_xfoil_polar('polar-aoas.dat')

    numCase = result1['numCase']
    AoAs= result1['AoAs']
    Xs  = result1['X']
    Cps = result1['Cp']
    Cfs = result1['Cf']
    
    with open('Cp-aoas.dat', 'w') as f:
        
        f.write('Variables= X Cp Cf \n')
        
        for k in range(numCase):
            f.write('zone T= "%f" \n'%(AoAs[k]))

            xx = Xs[k][1]
            nn = xx.shape[0]
            for i in range(nn):
                f.write('% 20.10f  %20.10f  %20.10f \n'%(
                    xx[nn-i-1], Cps[k][1][nn-i-1], Cfs[k][1][nn-i-1]))
                
            xx = Xs[k][0]
            nn = xx.shape[0]
            for i in range(nn-1):
                f.write('% 20.10f  %20.10f  %20.10f \n'%(
                    xx[i+1], Cps[k][0][i+1], Cfs[k][0][i+1]))
            
            f.write('\n')
    
    os.remove('airfoil.dat')
    os.remove('raw-aoas.bin')

    plt.figure(figsize=(12,8))
    plt.subplot(2,2,1)
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.subplot(2,2,2)
    plt.plot(result2['CLs'], result2['CDs'],  'b*-')
    plt.plot(result2['CLs'], result2['CDps'], 'g*-')
    plt.xlabel('CL')
    plt.ylabel('Cd')

    plt.subplot(2,2,3)
    for k in range(numCase):
        plt.plot(Xs[k][0], -Cps[k][0], 'b')
        plt.plot(Xs[k][1], -Cps[k][1], 'b')
    plt.xlabel('X')
    plt.ylabel('-Cp')
    
    plt.subplot(2,2,4)
    for k in range(numCase):
        plt.plot(Xs[k][0], Cfs[k][0], 'b')
        plt.plot(Xs[k][1], Cfs[k][1], 'b')
    plt.xlabel('X')
    plt.ylabel('Cf')
    
    plt.savefig('xfoil-aoas.png', dpi=300)
    
    