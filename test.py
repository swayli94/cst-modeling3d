
import os
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

from cst_modeling.surface import Surface


if __name__ == "__main__":

    os.system('cls')
    print()




    wing = Surface(n_sec=9, n_cst=7, name='Wing',nn=101, ns=51)
    wing.read_setting('Wing.txt')
    # wing.add_sec(location=[8, 12.0])
    wing.geo(split=True, showfoil=True)
    wing.bend(isec0=7, isec1=8, leader=[[22.9, 1.2, 27.1, 0.75]], kx=[0.4, 1.0], ky=[0.0, 1.6])
    # toCylinder

    wing.output_tecplot(fname='Wing.dat', one_piece=False)
    # wing.output_plot3d(fname='Wing.grd')

    # wing.flip(axis='+X +Z')
    # wing.plot()

    
    print()
    exit()


    u = np.zeros(3)      # independent variable list, e.g., z
    v = np.zeros(3)     # dependent variable list  , e.g., x, y
    slope0 = (0.0)
    slope1 = (0.0)

    u[0] = 0.0
    u[1] = 0.5
    u[2] = 1.0
    v[0] = 0.0
    v[1] = 0.5
    v[2] = 1.0

    leader_curve = CubicSpline(u, v, bc_type=((1,slope0), (1,slope1)))

    xx = np.linspace(0,1,100)
    fig, ax = plt.subplots()

    yy = leader_curve(xx)

    print(leader_curve(0.5))

    ax.plot(xx, yy)
    plt.show()
    exit()


