import numpy as np
import matplotlib.pyplot as plt


file_dump = './dump/'


def airfoil3d():
    '''
    3D surface for airfoil meshing in ICEM CFD
    '''
    from cst_modeling.basic import plot3d_to_igs
    from cst_modeling.surface import Surface
    
    wing = Surface(n_sec=0, name='airfoil3d', nn=501, ns=11, projection=True)
    wing.read_setting('./files/swept-wing.txt', tail=0.01)
    wing.geo()
    
    wing.output_tecplot(fname=file_dump+'airfoil3d.dat', split=False)
    wing.output_plot3d(fname=file_dump+'airfoil3d.grd', split=False)
    plot3d_to_igs(fname=file_dump+'airfoil3d')
    
    ax = wing.plot(type='surface', show=False)
    
    ax.view_init(elev=120, azim=-90)
        
    plt.savefig('figures/airfoil3d.jpg', dpi=300)
    plt.close()
    
def transonic_wing():
    '''
    Swept wing for transonic jet 
    '''
    from cst_modeling.surface import Surface
    
    #* Define a [Surface] object, which is called 'transonic-wing'
    #  It has n_sec control sections.
    #  Each section's upper/lower surface are constructed by CST method.
    #  The chord-wise has nn points in the upper/lower surface.
    #  The span-wise has ns points in each surface between two adjacent sections.
    wing = Surface(n_sec=6, name='transonic-wing', nn=101, ns=52)

    #* Read settings from file
    #  The settings of `wing` object is under its name 'transonic-wing'
    #  The first part is the parameters for the layout
    #  The second part is the CST coefficients of the upper/lower surfaces of each section
    wing.read_setting('./files/swept-wing.txt', tail=0.01)

    #* Add an auxiliary section
    wing.add_sec([5.0], axis='Z')

    #* Constructs the surfaces between sections
    wing.geo()

    #* Output tecplot format to fname
    # one_piece is an option for combining all surfaces in different sections into one piece
    # split is an option for splitting the upper/lower surfaces
    wing.output_tecplot(fname=file_dump+'transonic_wing.dat', one_piece=False, split=True)
    wing.output_plot3d(fname=file_dump+'transonic_wing.grd', split=True)

    #* Output control sections
    #  3D curves by default, 2D curves when `TwoD=True`.
    wing.output_section(fname=file_dump+'transonic_wing-sections.dat', TwoD=False)

    #* Plot figure
    ax = wing.plot(type='surface', show=False)
    ax.view_init(elev=150, azim=-90)
        
    plt.savefig('figures/transonic_wing.jpg', dpi=300)
    plt.close()
    
def transonic_wing_winglet():
    '''
    Swept wing for transonic jet (with wing tip)
    '''
    from cst_modeling.surface import Surface
    
    wing = Surface(n_sec=6, name='transonic-wing', nn=101, ns=52)
    
    wing.read_setting('./files/swept-wing.txt', tail=[0.05, 0.05, 0.05, 0.05, 0.05, 0.01])

    wing.geo()

    wing.bend(4, 5, leader=[[21.0, 2.1, 30.0, 1.6]], 
              kx=[0.6983, 4.0], ky=[0.1043, 1.10], rot_x=True)

    wing.smooth(0,2)
    wing.smooth(2,4)

    wing.output_tecplot(fname=file_dump+'transonic_wing_winglet.dat', 
                        one_piece=True, split=True)

    ax = wing.plot(type='surface', show=False)
    ax.view_init(elev=150, azim=-90)
    
    plt.savefig('figures/transonic_wing_winglet.jpg', dpi=300)
    plt.close()


if __name__ == '__main__':
    
    airfoil3d()
    
    transonic_wing()
    
    transonic_wing_winglet()
    