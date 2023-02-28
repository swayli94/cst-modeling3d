'''

'''
import numpy as np
import matplotlib.pyplot as plt
from cst_modeling.section import cst_curve


def cst_base_function():
    '''
    Plot CST base function
    '''
    n_cst = 7
    plt.figure()

    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst, xn1=0.5, xn2=1.0)
        plt.plot(x, y)

    plt.text(0.5, 0.14, 'xn1=0.5, xn2=1.0', fontsize=16)
    plt.savefig('figures/cst_base_function-1.jpg', dpi=300)
    plt.close()


    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst, xn1=0.1, xn2=1.0)
        plt.plot(x, y)

    plt.text(0.5, 0.5, 'xn1=0.1, xn2=1.0', fontsize=16)
    plt.savefig('figures/cst_base_function-2.jpg', dpi=300)
    plt.close()







if __name__ == '__main__':
    
    cst_base_function()



