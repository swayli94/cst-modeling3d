
import numpy as np
import matplotlib.pyplot as plt

def base_shape():
    '''
    Base shape of RoundTipSection
    '''

    from cst_modeling.section import RoundTipSection
    
    xx = np.linspace(0, 1, 501)
    
    x_, y_ = RoundTipSection.base_shape(xx, x_LE=0, x_TE=1, l_LE=0.1, l_TE=0.1,
                                        r_LE=0.02, r_TE=0.02, h=0.03, i_split=None)
    dy_    = RoundTipSection.base_camber(x_, a_LE=10, a_TE=-10)
    
    
    plt.figure()
    plt.plot(x_, y_, 'k')
    plt.plot(x_, y_+dy_, 'b--')
    
    plt.legend(['base shape', 'base shape with camber'])
    
    plt.plot(x_, -y_, 'k')
    plt.plot(x_, -y_+dy_, 'b--')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.20, 0.20))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/base_shape.jpg', dpi=300)
    plt.close()

def fitting_blade():
    '''
    Fitting blade (DFVLR L030-4)
    '''
    
    from scipy.interpolate import interp1d
    from cst_modeling.section import dist_clustcos, RoundTipSection, cst_foil_fit, cst_foil

    with open('data/blade-ref.dat', 'r') as f:
        
        lines = f.readlines()

        xu = []
        yu = []
        xl = []
        yl = []
        
        i = 0
        for line in lines[2:]:
            
            line = line.split()
            
            if len(line)==0:
                continue
            elif len(line)>2:
                i += 1
                continue
                    
            if i%2==1:
                xu.append(float(line[0]))
                yu.append(float(line[1]))
            else:
                xl.append(float(line[0]))
                yl.append(float(line[1]))

        xx = dist_clustcos(501)
        fu = interp1d(xu, yu, kind='linear')
        fl = interp1d(xl, yl, kind='linear')
        yu = fu(xx)
        yl = fl(xx)

    #* ==================================================================
    #* Manually try out base shape coefficients

    BaseLERatio     = 0.4
    BaseTERatio     = 0.4
    BaseAbsThick    = 0.04
    BaseRelRadiusLE = 0.02
    BaseRelRadiusTE = 0.03
    
    BaseAlphaLE = 6
    BaseAlphaTE = -10

    x_, y_ = RoundTipSection.base_shape(xx, 0.0, 1.0, BaseLERatio, BaseTERatio, 
                                        BaseRelRadiusLE, BaseRelRadiusTE, BaseAbsThick/2.0, i_split=None)
    
    dy_ = RoundTipSection.base_camber(xx, a_LE=BaseAlphaLE, a_TE=BaseAlphaTE)


    #* ==================================================================
    #* Fit with CST
    xn1 = 0.1
    xn2 = 0.1

    yu_ = yu - y_ + dy_
    yl_ = yl + y_ + dy_
    
    cst_u, cst_l = cst_foil_fit(xx, yu_, xx, yl_, n_cst=7, xn1=xn1, xn2=xn2)
    
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, linewidth=200)
    print(cst_u)
    print(cst_l)
    print()
    
    _, yuc, ylc, _, _ = cst_foil(xx.shape[0], cst_u, cst_l, xx, t=None, xn1=xn1, xn2=xn2)

    plt.figure(figsize=(10,4))
    plt.plot(xx, yu, 'k')
    plt.plot(xx, y_+dy_, 'b--')
    plt.plot(xx, yuc + y_ - dy_, 'r--')
    plt.plot(xx, 10*np.abs(yu_-yuc), 'g')
    
    plt.legend(['reference', 'base shape', 'fitting', '10*error'])
    
    plt.plot(xx, yl, 'k')
    plt.plot(xx,-y_+dy_, 'b--')
    plt.plot(xx, ylc - y_ - dy_, 'r--')
    plt.plot(xx, 10*np.abs(ylc-ylc), 'g')
    
    plt.xlim((-0.05, 1.05))
    plt.ylim((-0.01, 0.06))
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/fitting_blade.jpg', dpi=300)
    plt.close()



if __name__ == '__main__':

    base_shape()
    
    fitting_blade()
    