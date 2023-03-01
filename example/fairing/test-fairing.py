

import numpy as np
from cst_modeling.surface import BasicSurface

from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


if __name__ == "__main__":

    fairing = BasicSurface(n_sec=8, name='Fairing', nn=51, ns=101, projection=False)
    fairing.read_setting('Fairing.txt')


    # Read raw data
    with open('..\\fuselage\\tube-curve.dat', 'r') as f:
        X_end = []
        Y_end = []
        lines = f.readlines()
        for i in range(len(lines)-2):
            line = lines[i+2].split()
            X_end.append(float(line[0]))
            Y_end.append(float(line[1]))

    X_end = np.array(X_end)
    Y_end = np.array(Y_end)

    for i in range(fairing.n_sec):
        fairing.secs[i].xx = X_end.copy()
        fairing.secs[i].yy = Y_end.copy()


    #* Control section curve
    X0 = [X_end[-1], -2.30, -1.80, -0.80, X_end[0]]
    Y0 = [0.0,        0.40,  0.55,  0.10, 0.0]
    curve = CubicSpline(X0, Y0, bc_type=((1,0.0), (1,0.0)))

    y_increment = curve(X_end)

    if True:
        #! Plot curves on screen
        plt.plot(X_end, Y_end, 'r')
        plt.plot(X0, Y0, 'go')
        plt.plot(X_end, y_increment, 'g')
        plt.plot(X_end, Y_end+y_increment, 'b')
        plt.axis('equal')
        plt.legend(['baseline curve', 'control points', 'incremental curve', 'final curve'])
        plt.show()

    fairing.secs[2].yy += y_increment
    fairing.secs[3].yy += y_increment
    fairing.secs[4].yy += y_increment*1.2
    fairing.secs[5].yy += y_increment*1.2


    #* Build geometry

    fairing.geo()
    
    fairing.smooth(1, 3, smooth0=True, smooth1=True)
    fairing.smooth(5, 7, smooth0=True, smooth1=True)

    fairing.flip(axis='+Z +Y')

    fairing.output_tecplot(fname='fairing.dat')

