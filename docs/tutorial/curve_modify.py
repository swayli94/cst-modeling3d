
import numpy as np
import matplotlib.pyplot as plt

def bumps():
    '''
    Bump function
    '''
    from cst_modeling.section import bump_function

    xcs = [0.02, 0.05, 0.5, 0.5, 0.98]
    ss =  [0.1,  0.5,  0.1, 1.0, 0.1]
    
    xx = np.linspace(0, 1, 501, endpoint=True)

    plt.figure(figsize=(10, 4))
    plt.subplot(121)

    for i in range(len(xcs)):

        plt.plot(xx, bump_function(xx, xcs[i], 0.5, ss[i], kind='G'))

    plt.title('Gaussian bump')
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.subplot(122)
    for i in range(len(xcs)):

        plt.plot(xx, bump_function(xx, xcs[i], 0.5, ss[i], kind='H'))

    plt.title('Hicks-Henne bump')
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/bumps.jpg', dpi=300)
    plt.close()
    
def foil_increment_curve():
    '''
    Foil with incremental curves
    '''
    
    from cst_modeling.section import cst_foil, bump_function, foil_increment_curve
    
    tmax = 0.11
    h = 0.01
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    plt.figure()
    
    x, yu, yl, _, _ = cst_foil(1001, cst_u, cst_l, x=None, t=tmax, tail=0.01)

    yu_i = bump_function(x, 0.01, h, 0.2, kind='H')
    yl_i = bump_function(x, 0.90, h, 0.8, kind='H')
    yu1, yl1 = foil_increment_curve(x, yu, yl, yu_i, yl_i, t=tmax)
    
    yu_i = bump_function(x, 0.30, h, 0.4, kind='H')
    yu2, yl2 = foil_increment_curve(x, yu, yl, yu_i, yl_i=None, t=tmax)
    
    
    plt.plot(x, yu, 'k')
    plt.plot(x, yu1, 'r')
    plt.plot(x, yu2, 'g--')
    
    plt.legend(['original', 'add bumps', 'add a bump near tmax location'])
    
    plt.plot(x, yl, 'k')
    plt.plot(x, yl1, 'r')
    plt.plot(x, yl2, 'g--')

    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.08, 0.12))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/foil_increment_curve.jpg', dpi=300)
    plt.close()

def foil_bump_modification():
    '''
    Foil with incremental curves
    '''
    
    from cst_modeling.section import cst_foil, foil_bump_modify
    
    tmax = 0.11
    h = 0.01
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    x, yu, yl, _, _ = cst_foil(1001, cst_u, cst_l, x=None, t=tmax, tail=0.01)

    yu1, yl1 = foil_bump_modify(x, yu,  yl,  0.01, h/tmax, 0.2,  1, n_cst=10, keep_tmax=True)
    yu1, yl1 = foil_bump_modify(x, yu1, yl1, 0.90, h/tmax, 0.8, -1, n_cst=10, keep_tmax=True)
    
    yu2, yl2 = foil_bump_modify(x, yu, yl, 0.30, h/tmax, 0.4, 1, n_cst=10, keep_tmax=True)


    plt.figure()
    plt.plot(x, yu, 'k')
    plt.plot(x, yu1, 'r')
    plt.plot(x, yu2, 'g--')
    
    plt.legend(['original', 'add bumps', 'add a bump near tmax location'])

    plt.plot(x, yl, 'k')
    plt.plot(x, yl1, 'r')
    plt.plot(x, yl2, 'g--')

    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.08, 0.12))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/foil_bump_modification.jpg', dpi=300)
    plt.close()


    
if __name__ == '__main__':
    
    bumps()
    
    foil_increment_curve()
    
    foil_bump_modification()
