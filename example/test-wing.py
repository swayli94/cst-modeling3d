
import os
import numpy as np
import matplotlib.pyplot as plt

# Run this in the directory where the folder cst_modeling is
from cst_modeling.surface import Surface


if __name__ == "__main__":

    os.system('cls')

    print('This is a example for constructing a wing')
    # Define a [Surface] object, which is called wing
    # It has n_sec control sections, and its name is 'Wing'
    # Each section's upper/lower surface are constructed by CST method
    # The chord-wise has nn points in the upper/lower surface
    # The span-wise has ns points in each surface between two adjcent sections
    wing = Surface(n_sec=9, name='Wing',nn=101, ns=51)

    # Read settings from file 'Wing.txt'
    # The settings of [wing] object is under its name 'Wing'
    # First part is the parameters fro the layout
    # Second part is the CST parameters of upper/lower surface of each section
    wing.read_setting('Wing.txt')

    # This function can add additional sections to specified locations
    # The new sections are linearly interploted from original control sections
    wing.add_sec(location=[8, 12.0], axis='Z')

    # This constructs the surfaces between sections
    # split and showfoil are options
    wing.geo(split=True, showfoil=True)

    # This is for constructing surfaces with curved leading edge lines
    # e.g., strut of the SBW aircraft, wing-let
    wing.bend(isec0=9, isec1=10, leader=[[22.9, 1.2, 27.1, 0.75]], kx=[0.4, 1.0], ky=[0.0, 1.6], rot_x=True)

    # This outputs the surface to fname in tecplot format
    # one_piece is an option for combining all surfaces in different sections into one piece
    wing.output_tecplot(fname='Wing.dat', one_piece=False)

    # This outputs the surface to fname in plot3d format
    wing.output_plot3d(fname='Wing.grd')

    # This is flipping the surface, it can do turing 90deg or mirror
    wing.flip(axis='+X +Z')

    # Plot the surface by python
    wing.plot()

    print()

    exit()



