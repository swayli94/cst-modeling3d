
from cst_modeling.surface import OpenSurface


if __name__ == "__main__":

    fairing = OpenSurface(n_sec=3, name='Fairing-simple', nn=51, ns=51, projection=False)

    fairing.read_setting('Fairing.txt')

    phi = [0.0, 90.0, 180.0]

    fairing.geo_axisymmetric(phi)

    fairing.smooth_axisymmetric(0, 2, phi, linear_TEx=True)

    fairing.output_tecplot(fname='Fairing-simple.dat')


