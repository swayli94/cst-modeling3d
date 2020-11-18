
from cst_modeling.surface import OpenSurface


if __name__ == "__main__":

    patch = OpenSurface(n_sec=3, name='Patch', nn=101, ns=101, project=False)

    patch.read_setting('Fan.txt')

    patch.geo(flip_x=False, update_sec=True)

    patch.smooth(0,2)

    patch.output_tecplot(fname='Patch.dat', one_piece=False)
