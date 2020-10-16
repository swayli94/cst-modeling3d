
import numpy as np

from cst_modeling.foil import cst_foil, cst_curve, cst_foil_fit

from matplotlib import pyplot as plt


if __name__ == "__main__":

    np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=200)

    #* CST shape functions (Kulfan, 2008)
    n_cst = 7
    plt.figure()

    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst)
        plt.plot(x, y)

    plt.show()


    #* Build an airfoil
    cst_u = np.array([ 0.699328, 0.701191, 0.918286, 0.806252, 1.233955, 0.874498, 1.141528])
    cst_l = np.array([-0.681143,-0.791295,-0.643583,-1.493055,-0.072058,-0.698530, 0.377973])

    x, yu, yl, t0, R0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)

    plt.figure()
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.show()


    #* Fit an airfoil

    cst_u, cst_l = cst_foil_fit(x, yu, x, yl, n_order=n_cst)

    print(cst_u)
    print(cst_l)
    print()

    yu = yu * 1.2
    cst_u, cst_l = cst_foil_fit(x, yu, x, yl, n_order=n_cst)

    print(cst_u)
    print(cst_l)
    print()

