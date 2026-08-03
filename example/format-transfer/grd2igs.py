'''
Format transfer: plot3d (*.xyz) to IGES (*.igs).

A demo surface is generated in plot3d format, then converted to IGES.
To convert your own surface, write it to `wing.xyz` in plot3d format
(or change `fname` below) and run this script.
'''
import numpy as np

from cst_modeling.math import surf_axisymmetric
from cst_modeling.io import output_plot3d, plot3d_to_igs


if __name__ == "__main__":

    #* Generate a demo surface (axisymmetric body) in plot3d format
    xx = np.linspace(0.0, 4.0, 41)
    rr = 1.0 - (np.linspace(0.0, 1.0, 41)*2.0 - 1.0)**2
    surf = surf_axisymmetric(xx, rr, phi0=0.0, phi1=360.0, ns=61)
    output_plot3d([surf[0]], [surf[1]], [surf[2]], fname='wing.xyz')

    #* Transfer wing.xyz to wing.igs
    plot3d_to_igs(fname='wing')
