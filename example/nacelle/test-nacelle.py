
from cst_modeling.surface import Surface


if __name__ == "__main__":

    if True:

        nacelle = Surface(n_sec=5, name='Nacelle-simple', nn=201, ns=101)

        nacelle.read_setting('Nacelle.txt', tail=0.1)

        phi = [0.0, 90.0, 180.0, 270.0, 360.0]

        nacelle.geo_axisymmetric(phi)

        nacelle.smooth_axisymmetric(0, 4, phi, linear_TEx=True)

        nacelle.output_tecplot(fname='Nacelle-simple.dat', one_piece=False, split=False)

    if True:
        
        def func_trans(tx, x_cri=0.6):
            if tx < x_cri:
                ratio = 0.0
            else:
                ratio = ((tx-x_cri)/(1-x_cri))**20
            return ratio

        nacelle = Surface(n_sec=7, name='Nacelle', nn=51, ns=51)

        nacelle.read_setting('Nacelle.txt', tail=0.02)

        phi = [0.0, 90.0, 135.0, 180.0, 225.0, 270.0, 360.0]

        nacelle.geo_axisymmetric(phi)

        nacelle.smooth_axisymmetric(0, 6, phi, linear_TEx=True, RTE=0.8, RTE_=0.78, func_trans=func_trans)

        nacelle.output_tecplot(fname='Nacelle.dat', one_piece=False, split=False)

