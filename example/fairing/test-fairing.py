

import numpy as np
from cst_modeling.foil import cst_curve
from cst_modeling.surface import BasicSurface

import matplotlib.pyplot as plt


if __name__ == "__main__":


    fairing = BasicSurface(n_sec=8, name='Fairing', nn=51, ns=51, project=False)
    fairing.read_setting('Fairing.txt')

    #* Fairing end curves
    with open('..\\fuselage\\tube-section.dat', 'r') as f:
        Y_end = []
        Z_end = []
        lines = f.readlines()
        for i in range(len(lines)-2):
            if i >= 102 and i<= 152:
                line = lines[i].split()
                Y_end.append(float(line[0]))
                Z_end.append(float(line[1]))

        Y_end = np.array(Y_end) - Y_end[0]
        Z_end = np.array(Z_end) - Z_end[0]

    for i in range(fairing.n_sec):
        fairing.secs[i].xx = Y_end.copy()
        fairing.secs[i].yy = Z_end.copy()

    #* Inner curves

    # Increment curve of yy
    n_cst = 5
    nx  = Y_end.shape[0]
    cc  = Y_end[-1]
    xx  = Y_end.copy()/cc
    
    cst = np.zeros(n_cst)
    cst[0] = 0.0
    cst[1] = 0.8
    cst[2] = 1.0
    cst[3] = 1.0

    _, y_i = cst_curve(nx, cst, x=xx)
    fairing.secs[2].yy += np.flip(y_i)*np.abs(cc)
    fairing.secs[3].yy += np.flip(y_i)*np.abs(cc)
    fairing.secs[4].yy += np.flip(y_i)*np.abs(cc)
    fairing.secs[5].yy += np.flip(y_i)*np.abs(cc)

    '''
    plt.plot(Y_end, Z_end, 'b')
    plt.plot(Y_end, fairing.secs[2].yy, 'r')
    plt.axis('equal')
    plt.show()
    '''

    #* Build geometry
    fairing.geo()

    fairing.smooth(1, 3, smooth0=True, smooth1=True)
    fairing.smooth(5, 7, smooth0=True, smooth1=True)

    fairing.flip(axis='+Z +Y')

    fairing.translate(dX=15.0, dY=0.0, dZ=2.6)

    fairing.output_tecplot(fname='fairing.dat')

