
from cst_modeling.surface import Surface


if __name__ == "__main__":

    blade = Surface(n_sec=6, name='Blade-simple',nn=101, ns=51, project=False)

    blade.read_setting('Fan.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.1, 0.05])

    blade.geo()

    blade.smooth(isec0=0, isec1=4)
    blade.smooth(isec0=4, isec1=5, smooth0=True)

    blade.Surf2Cylinder(flip=True)

    blade.output_tecplot(fname='Blade-simple.dat')

