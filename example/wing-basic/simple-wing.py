import os
import sys
sys.path.append('.')

from cst_modeling.surface2 import Surface
from cst_modeling.io import plot3d_to_igs


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])
    
    #* Define a [Surface] object, which is called wing
    #  It has n_sec control sections, and its name is 'Wing'
    #  Each section's upper/lower surface are constructed by CST method
    #  The chord-wise has nn points in the upper/lower surface
    #  The span-wise has ns points in each surface between two adjacent sections
    wing = Surface(n_sec=3, name='Wing-basic', nn=101, ns=51,
                    smooth_surface=False, smooth_sections=None)

    #* Read settings from file 'Wing.txt'
    #  The settings of [wing] object is under its name 'Wing'
    #  First part is the parameters for the layout
    #  Second part is the CST parameters of upper/lower surface of each section
    wing.read_setting(os.path.join(path, 'Wing.txt'))

    wing.prepare()
    
    wing.geo()


    #* Output tecplot format to fname
    # one_piece is an option for combining all surfaces in different sections into one piece
    # split is an option for splitting to upper and lower surfaces
    wing.output_tecplot(fname=os.path.join(path, 'Wing-basic.dat'), one_piece=False, split=True)

    wing.output_plot3d(fname=os.path.join(path, 'Wing-basic.xyz'), split=True)

    plot3d_to_igs(fname=os.path.join(path, 'Wing-basic'))
    
    wing.output_guide_curve(os.path.join(path, 'Wing-basic-guide-curve.dat'))


    #* Plot the wing surface on screen
    '''
    #  This is flipping the surface, it can do turing 90deg or mirror
    wing.flip(axis='+X +Z')
    wing.plot()
    '''
