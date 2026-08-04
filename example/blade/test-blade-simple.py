
from cst_modeling.basic import Surface


if __name__ == "__main__":

    blade = Surface(n_section=6, name='Blade-simple',nn=101, ns=51, projection=False)

    blade.read_setting('Fan.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.1, 0.05])

    blade.geo()

    blade.smooth(i_sec0=0, i_sec1=4)
    blade.smooth(i_sec0=4, i_sec1=5, smooth0=True)

    blade.surf_to_cylinder(flip=True)

    blade.output_tecplot(fname='Blade-simple.dat')

