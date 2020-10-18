
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
    fan = Surface(n_sec=7, name='Blade-basic',nn=51, ns=51, project=False)

    # Read settings from file 'Fan.txt'
    # The settings of [fan] object is under its name 'Fan'
    # First part is the parameters fro the layout
    # Second part is the CST parameters of upper/lower surface of each section
    fan.read_setting('Fan.txt', tail=0.86)

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




