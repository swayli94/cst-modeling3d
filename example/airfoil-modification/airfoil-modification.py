'''
Modify airfoil geometry with bump functions
'''
import os
import sys
sys.path.append('.')

import numpy as np
from matplotlib import pyplot as plt

from cst_modeling.section import cst_foil
from cst_modeling.foil import FoilGeoFeatures, FoilModification


if __name__ == '__main__':

    path = os.path.dirname(sys.argv[0])

    #* Initialize an airfoil
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])


    #* Leading edge modification
    if True:

        for tail, rLE, width_bump in [(0.0, 0.015, 1.0), (0.004, 0.005, 0.8)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=tail)
            
            geo_old = FoilGeoFeatures(x, yu, yl)

            modify = FoilModification(x, yu, yl)
            
            _, _, _, rLE_new = modify.set_leading_edge_radius(rLE_new=rLE, width_bump=width_bump)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)

            plt.figure(figsize=(16,8))
            
            plt.plot(x, yu, 'k')
            plt.plot(x, yl, 'k', label='_nolegend_')
            
            plt.plot(x, geo_old.get_feature('thickness'), 'k--', lw=0.5)
            plt.plot(x, geo_old.get_feature('camber'), 'b--', lw=0.5)
            
            plt.plot(modify.x, modify.yu, 'r', lw=1.0)
            plt.plot(modify.x, modify.yl, 'r', lw=1.0, label='_nolegend_')
            
            plt.plot(modify.x, geo_new.get_feature('thickness'), 'r--', lw=0.5)

            
            plt.xlim((-0.2, 1.2))
            plt.ylim((-0.2, 0.2))
            plt.axis('equal')
            plt.legend(['Airfoil (old)', 'Thickness (old)', 'Camber (old)', 'Airfoil (new)', 'Thickness (new)'])
            
            plt.title('Leading edge modification: %.3f -> %.3f'%(rLE_old, rLE_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-rLE-%.3f.png'%(rLE_new)), dpi=300)
            
            
    #* Trailing edge wedge angle modification
    if True:

        for tail, wedge_angle, width_bump in [(0.0, 15.0, 0.2), (0.004, 2.0, 0.6)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=tail)
            
            geo_old = FoilGeoFeatures(x, yu, yl)
            wedge_angle_old = geo_old.get_trailing_edge_wedge_angle()

            modify = FoilModification(x, yu, yl)
            
            modify.set_trailing_edge_wedge_angle(wedge_angle_new=wedge_angle, width_bump=width_bump)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            wedge_angle_new = geo_new.get_trailing_edge_wedge_angle()

            plt.figure(figsize=(16,8))
            
            plt.plot(x, yu, 'k')
            plt.plot(x, yl, 'k', label='_nolegend_')
            
            plt.plot(x, geo_old.get_feature('thickness'), 'k--', lw=0.5)
            plt.plot(x, geo_old.get_feature('camber'), 'b--', lw=0.5)
            
            plt.plot(modify.x, modify.yu, 'r', lw=1.0)
            plt.plot(modify.x, modify.yl, 'r', lw=1.0, label='_nolegend_')
            
            plt.plot(modify.x, geo_new.get_feature('thickness'), 'r--', lw=0.5)

            
            plt.xlim((-0.2, 1.2))
            plt.ylim((-0.2, 0.2))
            plt.axis('equal')
            plt.legend(['Airfoil (old)', 'Thickness (old)', 'Camber (old)', 'Airfoil (new)', 'Thickness (new)'])
            
            plt.title('Trailing edge wedge angle: %.1f -> %.1f'%(wedge_angle_old, wedge_angle_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-TE-wedge-angle-%.1f.png'%(wedge_angle_new)), dpi=300)
            
            
    #* Trailing edge slope angle modification
    if True:

        for tail, slope_angle, width_bump in [(0.0, 15.0, 0.8), (0.004, 1.0, 0.2)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=tail)
            
            geo_old = FoilGeoFeatures(x, yu, yl)
            slope_angle_old = geo_old.get_trailing_edge_slope_angle()

            modify = FoilModification(x, yu, yl)
            
            modify.set_trailing_edge_slope_angle(slope_angle_new=slope_angle, width_bump=width_bump)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            slope_angle_new = geo_new.get_trailing_edge_slope_angle()

            plt.figure(figsize=(16,8))
            
            plt.plot(x, yu, 'k')
            plt.plot(x, yl, 'k', label='_nolegend_')
            
            plt.plot(x, geo_old.get_feature('thickness'), 'k--', lw=0.5)
            plt.plot(x, geo_old.get_feature('camber'), 'b--', lw=0.5)
            
            plt.plot(modify.x, modify.yu, 'r', lw=1.0)
            plt.plot(modify.x, modify.yl, 'r', lw=1.0, label='_nolegend_')
            
            plt.plot(modify.x, geo_new.get_feature('camber'), 'r--', lw=0.5)

            
            plt.xlim((-0.2, 1.2))
            plt.ylim((-0.2, 0.2))
            plt.axis('equal')
            plt.legend(['Airfoil (old)', 'Thickness (old)', 'Camber (old)', 'Airfoil (new)', 'Camber (new)'])
            
            plt.title('Trailing edge slope angle: %.1f -> %.1f'%(slope_angle_old, slope_angle_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-TE-slope-angle-%.1f.png'%(slope_angle_new)), dpi=300)
            
    
    #* Change airfoil camber
    if True:
        
        for xc, h, w in [(0.3, 0.01, 1.0), (0.8, 0.005, 0.6)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)

            modify = FoilModification(x, yu, yl)
            
            modify.add_bump_to_camber(xc, h, w, kind='H', keep_tmax=True)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            
            plt.figure(figsize=(16,8))
            
            plt.plot(x, yu, 'k')
            plt.plot(x, yl, 'k', label='_nolegend_')
            
            plt.plot(x, geo_old.get_feature('thickness'), 'k--', lw=0.5)
            plt.plot(x, geo_old.get_feature('camber'), 'b--', lw=0.5)
            
            plt.plot(modify.x, modify.yu, 'r', lw=1.0)
            plt.plot(modify.x, modify.yl, 'r', lw=1.0, label='_nolegend_')
            
            plt.plot(modify.x, geo_new.get_feature('camber'), 'r--', lw=0.5)

            
            plt.xlim((-0.2, 1.2))
            plt.ylim((-0.2, 0.2))
            plt.axis('equal')
            plt.legend(['Airfoil (old)', 'Thickness (old)', 'Camber (old)', 'Airfoil (new)', 'Camber (new)'])

            plt.savefig(os.path.join(path, 'airfoil-modify-camber-x-%.1f.png'%(xc)), dpi=300)
    
      