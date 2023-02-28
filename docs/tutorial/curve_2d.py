
import numpy as np
import matplotlib.pyplot as plt


def cst_base_function():
    '''
    Plot CST base function
    '''
    from cst_modeling.section import cst_curve
    
    n_cst = 7
    plt.figure(figsize=(9, 4))
    plt.subplot(121)

    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst, xn1=0.5, xn2=1.0)
        plt.plot(x, y)

    plt.text(0.2, 0.14, 'xn1=0.5, xn2=1.0', fontsize=16)
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.subplot(122)
    for i in range(n_cst):
        cst = np.zeros(n_cst)
        cst[i] = 1.0

        x, y = cst_curve(101, cst, xn1=0.1, xn2=1.0)
        plt.plot(x, y)

    plt.text(0.2, 0.5, 'xn1=0.1, xn2=1.0', fontsize=16)
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/cst_base_function.jpg', dpi=300)
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
    
    plt.xlabel('X')

    plt.savefig('figures/cosine_distribution.jpg', dpi=300)
    plt.close()

def cst_airfoil():
    '''
    Build an airfoil
    '''
    from cst_modeling.section import cst_foil
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    #* The influence of t and tail
    
    plt.figure()
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)
    plt.plot(x, yu, 'k')
    plt.plot(x, yl, 'k')
    plt.text(0.2, 0.10, 'tmax= %.3f, rLE= %.4f, tail=0.00'%(t0, r0), fontsize=12, color='k')
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.02)
    plt.plot(x, yu, 'r--')
    plt.plot(x, yl, 'r--')
    plt.text(0.2, 0.09, 'tmax= %.3f, rLE= %.4f, tail=0.02'%(t0, r0), fontsize=12, color='r')
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=0.05, tail=0.0)
    plt.plot(x, yu, 'b')
    plt.plot(x, yl, 'b')
    plt.text(0.2, 0.08, 'tmax= %.3f, rLE= %.4f, tail=0.00'%(t0, r0), fontsize=12, color='b')
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=0.05, tail=0.01)
    plt.plot(x, yu, 'g--')
    plt.plot(x, yl, 'g--')
    plt.text(0.2, 0.07, 'tmax= %.3f, rLE= %.4f, tail=0.01'%(t0, r0), fontsize=12, color='g')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.08, 0.12))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/cst_airfoil-t-tail.jpg', dpi=300)
    plt.close()
    
    
    #* The influence of CST parameters xn1 and xn2
    
    plt.figure()
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
    plt.plot(x, yu, 'k')
    plt.plot(x, yl, 'k')
    plt.text(0.1, 0.10, 'tmax= %.3f, rLE= %.4f, xn1=%.2f, xn2=%.2f'%(t0, r0, 0.5, 1.0), fontsize=12, color='k')
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0, xn1=0.1, xn2=1.0)
    plt.plot(x, yu, 'r')
    plt.plot(x, yl, 'r')
    plt.text(0.1, 0.09, 'tmax= %.3f, rLE= %.4f, xn1=%.2f, xn2=%.2f'%(t0, r0, 0.1, 1.0), fontsize=12, color='r')
    
    x, yu, yl, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0, xn1=0.5, xn2=0.5)
    plt.plot(x, yu, 'g--')
    plt.plot(x, yl, 'g--')
    plt.text(0.1, 0.08, 'tmax= %.3f, rLE= %.4f, xn1=%.2f, xn2=%.2f'%(t0, r0, 0.5, 0.5), fontsize=12, color='g')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.08, 0.12))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/cst_airfoil-xn1-xn2.jpg', dpi=300)
    plt.close()

def curve_curvature():
    '''
    Calculate the curvature distribution of a curve
    '''
    from cst_modeling.section import cst_curve, curve_curvature

    cst = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    x, y = cst_curve(501, cst, xn1=0.5, xn2=1.0)
    curvature = curve_curvature(x, y)

    plt.figure()
    plt.plot(x, y*20, 'k')
    plt.plot(x, -curvature, 'g--')
    plt.legend(['20*curve', '- curvature'])
    plt.xlim((-0.05, 1.05))
    plt.ylim((-1, 10))
    plt.xlabel('X')
    plt.savefig('figures/curve_curvature.jpg', dpi=300)
    plt.close()

def fitting_airfoil():
    '''
    Fitting an airfoil with CST
    '''

    from cst_modeling.section import cst_foil, cst_foil_fit
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    x0, yu0, yl0, _, _ = cst_foil(1001, cst_u, cst_l)

    cst_u, cst_l = cst_foil_fit(x0, yu0, x0, yl0, n_cst=10, xn1=0.3, xn2=1.0)
    
    _, yu1, yl1, _, _ = cst_foil(x0.shape[0], cst_u, cst_l, x=x0, xn1=0.3, xn2=1.0)
    
    plt.figure()
    plt.plot(x0, yu0, 'k')
    plt.plot(x0, yu1, 'g--')
    plt.plot(x0, 100*np.abs(yu1-yu0), 'r')

    plt.legend(['original', 'fitting', '100*error'])
    
    plt.plot(x0, yl0, 'k')
    plt.plot(x0, yl1, 'g--')
    plt.plot(x0, 100*np.abs(yl1-yl0), 'r')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.25, 0.25))
    plt.xlabel('X')
    plt.savefig('figures/fitting_airfoil.jpg', dpi=300)
    plt.close()

def fitting_curve():
    '''
    Fitting a curve with CST
    '''
    from cst_modeling.basic import transform
    from cst_modeling.section import cst_curve, fit_curve_with_twist, fit_curve_partial
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    
    x, y = cst_curve(1001, cst_u)
    
    x_new, _, y_new, _ = transform(x, x, y, y, scale=1.5, rot=20, dx=-0.5, dy=0.2)

    coef, chord, twist, _ = fit_curve_with_twist(x_new, y_new, n_cst=10, xn1=0.3, xn2=1.0)
    
    _, y_ = cst_curve(x.shape[0], coef, x=x, xn1=0.3, xn2=1.0)
    
    plt.figure()
    plt.plot(x_new, y_new, 'b')
    plt.plot(x, y_, 'g')
    plt.plot(x, 100*np.abs(y-y_), 'r')
    plt.legend(['original', 'fitting', '100*error'])
    plt.text(-0.2, 0.3, 'scale= %.2f, rotation= %.2f deg'%(chord, twist), fontsize=12, color='k')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('figures/fitting_curve.jpg', dpi=300)
    plt.close()
    

    coef = fit_curve_partial(x, y, ip0=200, ip1=501, ic0=2, ic1=11, n_cst=20, xn1=0.3, xn2=1.0)
    
    _, y_ = cst_curve(x.shape[0], coef, x=x, xn1=0.3, xn2=1.0)
    
    plt.figure()
    plt.plot(x, y, 'k--')
    plt.plot(x[200:501], y[200:501], 'k')
    plt.plot(x[200:501], y_[200:501], 'g--')
    plt.plot(x[200:501], 100*np.abs(y[200:501]-y_[200:501]), 'r')
    plt.legend(['original complete', 'original', 'fitting', '100*error'])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('figures/fitting_curve_partial.jpg', dpi=300)
    plt.close()
    

if __name__ == '__main__':
    
    cst_base_function()

    cosine_distribution()

    cst_airfoil()

    curve_curvature()
    
    fitting_airfoil()
    
    fitting_curve()
    