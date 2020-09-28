
import numpy as np
from scipy.interpolate import CubicSpline

# Run this in the directory where the folder cst_modeling is
from cst_modeling.surface import Surface


if __name__ == "__main__":

    print('This is a example for constructing a fan')
    # Define a [Surface] object, which is called fan
    # It has n_sec control sections, and its name is 'Fan'
    # Each section's upper/lower surface are constructed by CST method (has n_cst parameters)
    # The chord-wise has nn points in the upper/lower surface
    # The span-wise has ns points in each surface between two adjcent sections
    fan = Surface(n_sec=7, name='Fan',nn=51, ns=51, project=False)

    # Read settings from file 'Fan.txt'
    # The settings of [fan] object is under its name 'Fan'
    # First part is the parameters fro the layout
    # Second part is the CST parameters of upper/lower surface of each section
    fan.read_setting('Fan.txt', tail=0.86)

    #!===================================================
    if True:
        #! This part modifies the layout parameters of fan
        #! These coordinates are the real value of the fan geometry
        #! The XLE, XTE, chord, etc. are values in the converted space
        # Leading edge
        x = np.array([15.180, 18.083, 18.222, 18.006, 17.168, 15.544, 20.749])
        y = np.array([13.292, 17.263, 21.259, 23.995, 28.183, 32.693, 31.469])
        z = np.array([15.057, 14.767, 16.453, 17.622, 19.434, 21.420, 22.103])
        XLE, YLE, ZLE = Surface.fromCylinder(x, y, z, flip=True)

        # Trailing edge center
        x = np.array([20.961, 25.417, 28.262, 29.943, 32.030, 33.346, 33.894])
        y = np.array([-3.109, -0.074,  3.029,  5.437,  9.627, 14.081, 17.136])
        z = np.array([-3.425, -2.203, -1.209, -0.552,  0.436,  1.340,  7.899])
        XTE, YTE, ZTE = Surface.fromCylinder(x, y, z, flip=True)

        for i in range(fan.n_sec):
            chord = np.linalg.norm([XLE[i]-XTE[i], YLE[i]-YTE[i], ZLE[i]-ZTE[i]])
            twist = np.arctan((YTE[i]-YLE[i])/(XTE[i]-XLE[i]))*180/np.pi

            fan.secs[i].xLE = XLE[i]
            fan.secs[i].yLE = YLE[i]
            fan.secs[i].zLE = ZLE[i]
            fan.secs[i].chord = chord
            fan.secs[i].twist = twist
    #!===================================================

    # This constructs the surfaces between sections
    # split and showfoil are options
    fan.geo(split=True, showfoil=False)

    # Smooth
    fan.smooth(isec0=0, isec1=5)

    # Convert back to cylinder
    fan.Surf2Cylinder(flip=True)

    # This outputs the surface to fname in tecplot format
    # one_piece is an option for combining all surfaces in different sections into one piece
    fan.output_tecplot(fname='Fan.dat', one_piece=False)

    # fan.output_plot3d(fname='Fan.grd')

    print()

    exit()



