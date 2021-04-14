
import numpy as np

from cst_modeling.foil import cst_foil, cst_curve, cst_foil_fit, foil_bump_modify, foil_increment

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
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    x, yu, yl, t0, R0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)
    print(t0)

    plt.figure()
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.07, 0.07))
    plt.show()


    #* Fit an airfoil
    cst_u_, cst_l_ = cst_foil_fit(x, yu, x, yl, n_cst=n_cst)

    print(cst_u_)
    print(cst_l_)
    print()

    yu_ = yu * 1.2
    cst_u_, cst_l_ = cst_foil_fit(x, yu_, x, yl, n_cst=n_cst)

    print(cst_u_)
    print(cst_l_)
    print()


    #* Airfoil bump modification
    xc = 0.7
    h  = 0.02
    s  = 0.6
    yu_, yl_ = foil_bump_modify(x, yu, yl, xc, h, s, side=1)

    plt.figure()
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.plot(x, yu_, 'r--')
    plt.plot(x, yl_, 'r--')
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.07, 0.07))
    plt.show()


    #* Airfoil increment
    cst_u_ = np.zeros(16)
    cst_l_ = np.zeros(16)
    cst_u_[12] = 0.05
    yu_, yl_ = foil_increment(x, yu, yl, cst_u_, cst_l_, t=t0)

    x1, y1 = cst_curve(101, cst_u_)

    plt.figure()
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.plot(x, yu_, 'r--')
    plt.plot(x, yl_, 'r--')
    plt.plot(x1, y1, 'g--')
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.07, 0.07))
    plt.show()

