
import numpy as np
from cst_modeling.surface import OpenSurface


if __name__ == "__main__":

    #* ==========================================
    #* Suction side
    blade1 = OpenSurface(n_sec=6, name='Blade-suction',nn=101, ns=51, project=False)

    blade1.read_setting('Fan.txt')

    origins = blade1.read_cylinder_origins('Fan.txt')

    blade1.geo()

    # Smooth
    blade1.smooth(isec0=0, isec1=5)

    # Convert back to cylinder
    blade1.Surf2Cylinder(flip=True, origin=origins)

    # This outputs the surface to fname in tecplot format
    blade1.output_tecplot(fname='Blade-suction.dat')

    #* ==========================================
    #* Pressure side
    blade2 = OpenSurface(n_sec=6, name='Blade-pressure',nn=101, ns=51, project=False)

    blade2.read_setting('Fan.txt')

    origins = blade2.read_cylinder_origins('Fan.txt')

    blade2.geo()

    blade2.smooth(isec0=0, isec1=5)

    blade2.Surf2Cylinder(flip=True, origin=origins)

    blade2.output_tecplot(fname='Blade-pressure.dat')


