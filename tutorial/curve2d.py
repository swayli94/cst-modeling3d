'''

'''
import numpy as np
import matplotlib.pyplot as plt



def cst_base_function():
    '''
    Plot CST base function
    '''
    from cst_modeling.section import cst_curve
    
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

def cosine_distribution():
    '''
    Point distribution on x-axis
    '''
    from cst_modeling.section import dist_clustcos
    from sklearn.neighbors import KernelDensity

    x1 = dist_clustcos(101, a0=0.0079, a1=0.96, beta=1.0)
    x2 = dist_clustcos(101, a0=0.1000, a1=0.60, beta=1.0)
    x3 = dist_clustcos(101, a0=0.0079, a1=0.96, beta=1.2)
    
    kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(x1[:,None])
    log_density = kde.score_samples(x1[:,None])
    density = np.exp(log_density)

    plt.plot(x1, np.ones_like(x1)*2, 'k*')
    plt.plot(x1, density, 'k')
    plt.text(0.1, 2.2, 'a0=0.0079, a1=0.96, beta=1.0', fontsize=16, color='k')
    
    
    kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(x2[:,None])
    log_density = kde.score_samples(x2[:,None])
    density = np.exp(log_density)
    
    plt.plot(x2, np.ones_like(x1), 'g*')
    plt.plot(x2, density, 'g')
    plt.text(0.1, 1.2, 'a0=0.1000, a1=0.60, beta=1.0', fontsize=16, color='g')
    
    
    kde = KernelDensity(kernel='gaussian', bandwidth=0.01).fit(x3[:,None])
    log_density = kde.score_samples(x3[:,None])
    density = np.exp(log_density)
    
    plt.plot(x3, np.ones_like(x1)*3, 'r*')
    plt.plot(x3, density, 'r')
    plt.text(0.1, 3.2, 'a0=0.0079, a1=0.96, beta=1.2', fontsize=16, color='r')
    
    
    plt.savefig('figures/cosine_distribution.jpg', dpi=300)
    plt.close()



if __name__ == '__main__':
    
    cst_base_function()

    cosine_distribution()


