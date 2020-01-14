
import os
import numpy as np
import copy
import matplotlib.pyplot as plt

# Run this in the directory where the folder cst_modeling is
from cst_modeling.surface import Surface
from cst_modeling.foil import cst_foil_fit, rotate

def sortX(loc):
    return loc[0]

if __name__ == "__main__":

    # For turbomachinery, the sections are a 3D curve on cylinders
    # The origin is (x,y,z)=(0,0,0), or (r,t,z)=(0,0,0). (t~theta, rad)
    
    # This .py transfer the coordinates of curves on cylinders
    # to 2D sections on planes

    #TODO: ==================================================
    print()
    print('This is a example for transfer cylinder to plane')
    print()

    # Leading edge
    x = [15.180, 18.083, 18.222, 18.006, 17.168, 15.544, 20.749]
    y = [13.292, 17.263, 21.259, 23.995, 28.183, 32.693, 31.469]
    z = [15.057, 14.767, 16.453, 17.622, 19.434, 21.420, 22.103]
    XLE, YLE, ZLE = Surface.fromCylinder(x, y, z, flip=True)

    # Trailing edge upper
    x = [20.926, 25.423, 28.316, 30.028, 32.158, 33.512, 33.857]
    y = [-3.837, -0.504,  2.591,  5.010,  9.223, 13.689, 17.154]
    z = [-3.425, -2.192, -1.199, -0.539,  0.446,  1.353,  9.067]
    XTE1, YTE1, ZTE1 = Surface.fromCylinder(x, y, z, flip=True)

    # Trailing edge lower
    x = [20.996, 25.411, 28.208, 29.859, 31.902, 33.181, 33.931]
    y = [-2.382,  0.356,  3.467,  5.865, 10.032, 14.473, 17.118]
    z = [-3.425, -2.215, -1.219, -0.565,  0.427,  1.328,  6.731]
    XTE2, YTE2, ZTE2 = Surface.fromCylinder(x, y, z, flip=True)

    XTE = copy.deepcopy(XLE)
    YTE = copy.deepcopy(YLE)
    ZTE = copy.deepcopy(ZLE)
    chord = copy.deepcopy(XLE)
    twist = copy.deepcopy(XLE)
    for i in range(len(x)):
        XTE[i] = 0.5*(XTE1[i]+XTE2[i])
        YTE[i] = 0.5*(YTE1[i]+YTE2[i])
        ZTE[i] = 0.5*(ZTE1[i]+ZTE2[i])

        chord[i] = np.linalg.norm([XLE[i]-XTE[i], YLE[i]-YTE[i], ZLE[i]-ZTE[i]])
        twist[i] = np.arctan((YTE[i]-YLE[i])/(XTE[i]-XLE[i]))*180/np.pi

        print('%.3f  %.3f  %.3f  %.3f  %.3f'%(XLE[i], YLE[i], ZLE[i], chord[i], twist[i]))

    #TODO: ==============================================
    print()
    print('This is a example for getting CST parameters')
    print()

    n_sec = 6
    with open('ori-sections.dat', 'r') as f:
        line = f.readline()    # Variables= X Y Z

        for i in range(n_sec):

            fig = plt.figure(i)

            # Upper Surface
            line = f.readline().split()
            n_point = int(line[2])
            xu = []
            yu = []
            zu = []
            for j in range(n_point):
                line = f.readline().split()
                xu.append(float(line[0]))
                yu.append(float(line[1]))
                zu.append(float(line[2]))
            line = f.readline() # Empty line
            xu, yu, _ = Surface.fromCylinder(xu, yu, zu, flip=True)

            # Lower Surface
            line = f.readline().split()
            n_point = int(line[2])
            xl = []
            yl = []
            zl = []
            for j in range(n_point):
                line = f.readline().split()
                xl.append(float(line[0]))
                yl.append(float(line[1]))
                zl.append(float(line[2]))
            line = f.readline() # Empty line
            xl, yl, _ = Surface.fromCylinder(xl, yl, zl, flip=True)

            # Translation, scale, rotation
            for j in range(n_point):
                xu[j] = (xu[j]-XLE[i])/chord[i]
                xl[j] = (xl[j]-XLE[i])/chord[i]
                yu[j] = (yu[j]-YLE[i])/chord[i]
                yl[j] = (yl[j]-YLE[i])/chord[i]

            zz = None
            xu, yu, _ = rotate(xu, yu, zz, angle=-twist[i], axis='Z')
            xl, yl, _ = rotate(xl, yl, zz, angle=-twist[i], axis='Z')

            # Sort and fit
            point_u = []
            point_l = []
            for j in range(n_point):
                point_u.append([xu[j], yu[j]])
                point_l.append([xl[j], yl[j]])

            point_u.sort(key=sortX)
            point_l.sort(key=sortX)

            for j in range(n_point):
                xu[j] = point_u[j][0]
                yu[j] = point_u[j][1]
                xl[j] = point_l[j][0]
                yl[j] = point_l[j][1]

            coef_upp, coef_low = cst_foil_fit(xu, yu, xl, yl, n_order=7)

            print('Section %d ---------------'%(i+1))
            for ii in range(len(coef_upp)):
                print(' %.6f '%(coef_upp[ii]), end='')
            print()
            for ii in range(len(coef_low)):
                print(' %.6f '%(coef_low[ii]), end='')
            print()

            plt.plot(xu, yu, 'b')
            plt.plot(xl, yl, 'b')
            plt.show()


    exit()
