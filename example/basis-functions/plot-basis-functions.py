'''
Plot class shape functions and Hicks-Henne bump functions
'''
import os
import sys
sys.path.append('.')

import numpy as np

from cst_modeling.section import cst_curve, bump_function, RoundTipSection

from matplotlib import pyplot as plt


if __name__ == "__main__":

    n_cst = 7
    xx = np.linspace(0,1,1001,endpoint=True)
    
    path = os.path.dirname(sys.argv[0])
    
    #* CST shape functions (Kulfan, 2008)
    
    plt.figure()

    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst)
        plt.plot(x, y)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'CST shape functions $x_{n1}$=0.5, $x_{n2}$=1.0')
    plt.savefig(os.path.join(path, 'basis-functions-xn1-0.5-xn2-1.0.png'), dpi=300)

    plt.figure()

    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst, xn1=0.5, xn2=0.5)
        plt.plot(x, y)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'CST shape functions $x_{n1}$=0.5, $x_{n2}$=0.5')
    plt.savefig(os.path.join(path, 'basis-functions-xn1-0.5-xn2-0.5.png'), dpi=300)


    #* Round tip shape functions
    
    plt.figure()
    
    x, y = RoundTipSection.base_shape(xx, 0.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, i_split=None)
    
    dy = RoundTipSection.base_camber(x, a_LE=30, a_TE=-10)

    plt.plot(x,  y, 'b')
    plt.plot(x, -y, 'b', label='_nolegend_')
    
    plt.plot(x,  y+dy, 'r--', lw=0.5)
    plt.plot(x, -y+dy, 'r--', lw=0.5, label='_nolegend_')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.legend(['Round tip basis shape', 'Round tip shape with camber'])
    plt.savefig(os.path.join(path, 'round-tip-section-basis.png'), dpi=300)
    

    #* Hicks-Henne bump functions

    plt.figure(figsize=(10,10))
    
    plt.subplot(2,2,1)
    yy = bump_function(xx, xc=0.01, h=0.01, s=0.1, kind='H'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.01, h=0.01, s=0.5, kind='H'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.01, h=0.01, s=1.0, kind='H'); plt.plot(xx, yy, 'g')
    yy = bump_function(xx, xc=0.99, h=0.01, s=0.1, kind='H'); plt.plot(xx, yy, 'k--')
    yy = bump_function(xx, xc=0.99, h=0.01, s=0.5, kind='H'); plt.plot(xx, yy, 'b--')
    yy = bump_function(xx, xc=0.99, h=0.01, s=1.0, kind='H'); plt.plot(xx, yy, 'g--')
    
    plt.subplot(2,2,2)
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.4, kind='H'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.5, h=0.01, s=0.4, kind='H'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.8, h=0.01, s=0.4, kind='H'); plt.plot(xx, yy, 'g')
    
    plt.subplot(2,2,3)
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.2, kind='H'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.8, kind='H'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.2, h=0.01, s=1.0, kind='H'); plt.plot(xx, yy, 'g')
    
    plt.subplot(2,2,4)
    yy = bump_function(xx, xc=0.5, h=-0.002, s=0.6, kind='H'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.5, h= 0.005, s=0.6, kind='H'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.5, h= 0.010, s=0.6, kind='H'); plt.plot(xx, yy, 'g')
    
    plt.savefig(os.path.join(path, 'hicks-henne-bumps.png'), dpi=300)
    
    
    #* Hicks-Henne bump functions
    xx = np.linspace(0,1,201,endpoint=True)
    
    plt.figure(figsize=(10,10))
    
    plt.subplot(2,2,1)
    yy = bump_function(xx, xc=0.01, h=0.01, s=0.1, kind='G'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.01, h=0.01, s=0.5, kind='G'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.01, h=0.01, s=1.0, kind='G'); plt.plot(xx, yy, 'g')
    yy = bump_function(xx, xc=0.99, h=0.01, s=0.1, kind='G'); plt.plot(xx, yy, 'k--')
    yy = bump_function(xx, xc=0.99, h=0.01, s=0.5, kind='G'); plt.plot(xx, yy, 'b--')
    yy = bump_function(xx, xc=0.99, h=0.01, s=1.0, kind='G'); plt.plot(xx, yy, 'g--')
    
    plt.subplot(2,2,2)
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.4, kind='G'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.5, h=0.01, s=0.4, kind='G'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.8, h=0.01, s=0.4, kind='G'); plt.plot(xx, yy, 'g')
    
    plt.subplot(2,2,3)
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.2, kind='G'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.2, h=0.01, s=0.8, kind='G'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.2, h=0.01, s=1.0, kind='G'); plt.plot(xx, yy, 'g')
    
    plt.subplot(2,2,4)
    yy = bump_function(xx, xc=0.5, h=-0.002, s=0.6, kind='G'); plt.plot(xx, yy, 'k')
    yy = bump_function(xx, xc=0.5, h= 0.005, s=0.6, kind='G'); plt.plot(xx, yy, 'b')
    yy = bump_function(xx, xc=0.5, h= 0.010, s=0.6, kind='G'); plt.plot(xx, yy, 'g')
    
    plt.savefig(os.path.join(path, 'gaussian-bumps.png'), dpi=300)
    
    