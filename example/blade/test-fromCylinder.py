

import numpy as np
from cst_modeling.foil import fit_curve_with_twist, fromCylinder, find_circle_3p
from cst_modeling.surface import OpenSurface


if __name__ == "__main__":

    np.set_printoptions(formatter={'float': '{: 10.6f}'.format}, linewidth=200)

    # For turbomachinery, the sections are a 3D curve on cylinders
    # The origin is (x,y,z)=(0,0,0), or (r,t,z)=(0,0,0). (t~theta, rad)
    
    # This .py transfer the coordinates of curves on cylinders
    # to 2D sections on planes

    n_sec = 12

    blade = OpenSurface(n_sec=n_sec, name='Blade', nn=101, ns=101, project=False)

    CST = []
    origins = []

    print('Layout')
    with open('ori-sections.dat', 'r') as f:
        line = f.readline()    # Variables= X Y Z

        for i in range(n_sec):

            #* Original curve on cylinder
            line = f.readline().split()
            n_point = int(line[2])
            x = np.zeros(n_point)
            y = np.zeros(n_point)
            z = np.zeros(n_point)
            for j in range(n_point):
                line = f.readline().split()
                x[j] = float(line[0])
                y[j] = float(line[1])
                z[j] = float(line[2])
            line = f.readline() # Empty line

            if x[0]>x[-1]:
                x = np.flip(x)
                y = np.flip(y)
                z = np.flip(z)

            #* Locate Origin of cylinder
            ii = int(0.5*n_point)
            _, origin = find_circle_3p([x[0], y[0]], [x[ii], y[ii]], [x[-1], y[-1]])
            origins.append(origin)
            
            #* Convert to plane curve
            xx, yy, zz = fromCylinder(x, y, z, flip=True, origin=origin)

            XLE = xx[0]
            YLE = yy[0]
            ZLE = zz[0]

            #* CST coefficients
            cst, chord, twist, thick = fit_curve_with_twist(xx, yy, n_order=7)
            CST.append(cst)

            print(np.array([XLE, YLE, ZLE, chord, twist, thick]))

            blade.secs[i].xLE   = XLE
            blade.secs[i].yLE   = YLE
            blade.secs[i].zLE   = ZLE
            blade.secs[i].chord = chord
            blade.secs[i].twist = twist
            blade.secs[i].thick = thick
            blade.secs[i].cst   = cst.copy()

    print()
    print('CST Parameters')
    for cst in CST:
        print(cst)
        
    print()
    print('Cylinder Origins')
    for origin in origins:
        print(origin)

    blade.geo_secs()
    blade.Surf2Cylinder(flip=True, origin=origins)
    blade.output_section(TwoD=False)

