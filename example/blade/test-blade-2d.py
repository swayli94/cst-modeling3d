
'''
Fit a reference blade geometry (DFVLR L030-4) with RoundTipSection

'''

import numpy as np
from scipy.interpolate import interp1d

from cst_modeling.foil import dist_clustcos, RoundTipSection, cst_foil_fit, cst_foil

import matplotlib.pyplot as plt

    
if __name__ == '__main__':


    with open('blade-ref.dat', 'r') as f:
        
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

    yu_ = yu - y_ + dy_
    yl_ = yl + y_ + dy_
    
    cst_u, cst_l = cst_foil_fit(xx, yu_, xx, yl_, n_cst=10)
    
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format}, linewidth=200)
    print(cst_u)
    print(cst_l)
    print()
    
    x_, yuc, ylc, t0, R0 = cst_foil(xx.shape[0], cst_u, cst_l, xx, t=None)

    plt.plot(xx, yu, 'k')
    plt.plot(xx, y_+dy_, 'b--')
    plt.plot(xx, yuc + y_ - dy_, 'r--')
    plt.plot(xx, 10*(yu_-yuc), 'g')
    
    plt.legend(['reference', 'base shape', 'fitting', '10*error'])
    
    plt.plot(xx, yl, 'k')
    plt.plot(xx,-y_+dy_, 'b--')
    plt.plot(xx, ylc - y_ - dy_, 'r--')
    plt.plot(xx, 10*(ylc-ylc), 'g')
    
    plt.show()
            
            