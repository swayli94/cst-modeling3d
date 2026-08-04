
import numpy as np

from cst_modeling.io import output_curve
from cst_modeling.basic import BasicSurface

from scipy.interpolate import CubicSpline


if __name__ == "__main__":

    #* ============================================
    #* Nose
    #* ============================================

    nose = BasicSurface(n_section=7, name='Nose', nn=201, ns=51, projection=False)
    nose.read_setting('./files/jet_fuselage.txt')

    #* Top curve
    xc = np.array([0.00, 2.00, 3.50, 4.50, 6.00])
    yc = np.array([0.00, 0.90, 2.00, 2.60, 3.00])

    curv = CubicSpline(xc, yc, bc_type=((1, 0.0), (1, 0.0)))

    xx = np.linspace(xc[0], xc[-1], num=nose.nn)
    yy = curv(xx)

    an = 1.2
    bn = 0.8
    for i in range(xx.shape[0]):
        if xx[i] <= an:
            yy[i] += np.sqrt(1-(xx[i]/an-1)**2)*bn
        else:
            yy[i] += bn

    nose.sections[0].xx = xx/xc[-1]
    nose.sections[0].yy = yy/xc[-1]
    nose.sections[-1].xx = xx/xc[-1]
    nose.sections[-1].yy = yy/xc[-1]

    #* Side top curves
    xc = np.array([0.00,  2.00,  6.00])
    yc = np.array([0.00, -0.70, -2.20])

    curv = CubicSpline(xc, yc, bc_type=((1, 0.0), (1, 0.0)))

    xx = np.linspace(xc[0], xc[-1], num=nose.nn)
    yy = curv(xx)

    an = 1.2
    bn = 0.8
    for i in range(xx.shape[0]):
        if xx[i] <= an:
            yy[i] -= np.sqrt(1-(xx[i]/an-1)**2)*bn
        else:
            yy[i] -= bn

    nose.sections[1].xx = xx/xc[-1]
    nose.sections[1].yy = yy/xc[-1]

    nose.sections[-2].xx = nose.sections[1].xx
    nose.sections[-2].yy = nose.sections[1].yy

    #* Side curve and bottom curve
    xc = np.array([0.00,  2.00,  6.00])
    yc = np.array([0.00, -0.60, -1.40])

    curv = CubicSpline(xc, yc, bc_type=((1, 0.0), (1, 0.0)))

    xx = np.linspace(xc[0], xc[-1], num=nose.nn)
    yy = curv(xx)

    an = 2.0
    bn = 1.2
    for i in range(xx.shape[0]):
        if xx[i] <= an:
            yy[i] -= np.sqrt(1-(xx[i]/an-1)**2)*bn
        else:
            yy[i] -= bn

    nose.sections[2].xx = xx/xc[-1]
    nose.sections[2].yy = yy/xc[-1]

    nose.sections[-3].xx = nose.sections[2].xx
    nose.sections[-3].yy = nose.sections[2].yy

    nose.sections[3].xx = xx/xc[-1]
    nose.sections[3].yy = yy/xc[-1]

    phi = [0.0, 30.0, 90.0, 180.0, 270.0, 330.0, 360.0]

    nose.geo_axisymmetric(phi)

    nose.smooth_axisymmetric(0, 6, phi, linear_TEx=True)

    nose.output_tecplot(fname='./dump/nose.dat')

    #* Nose end curve
    Y = nose.surfaces[0][1][:,-1]
    Z = nose.surfaces[0][2][:,-1]
    for i in range(nose.n_section-2):
        Y = np.concatenate((Y, nose.surfaces[i+1][1][1:,-1]), axis=0)
        Z = np.concatenate((Z, nose.surfaces[i+1][2][1:,-1]), axis=0)

    scale = 5.6
    Y = (Y - Y[0]) / scale
    Z = Z / scale

    output_curve(Y, Z, fname='./dump/tube-section.dat', zone_id=0)

    #* ============================================
    #* Tube and aft body
    #* ============================================
    if True:

        tube = BasicSurface(n_section=3, name='Tube', nn=Y.shape[0], ns=51, projection=False)
        tube.read_setting('./files/jet_fuselage.txt')

        tube.sections[0].xx = Y
        tube.sections[0].yy = Z

        tube.sections[1].xx = tube.sections[0].xx
        tube.sections[1].yy = tube.sections[0].yy

        aa = np.linspace(0.0, 2*np.pi, num=tube.nn)
        tube.sections[2].xx = (np.cos(aa)-1)*0.5
        tube.sections[2].yy = np.sin(aa)*0.5

        tube.geo()

        tube.smooth(1, 2, smooth0=True)

        tube.flip(axis='+Z +Y')

        tube.translate(dx=6.0, dy=3.8)

        tube.output_tecplot(fname='./dump/tube.dat')
