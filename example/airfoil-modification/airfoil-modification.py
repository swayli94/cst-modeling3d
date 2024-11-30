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
    if False:

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
            
        
    #* Leading edge slope angle modification
    if True:

        for slope_angle, width_bump in [(3.0, 0.6), (6.0, 1.0)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)
            slope_angle_old = geo_old.get_leading_edge_slope_angle()

            modify = FoilModification(x, yu, yl)
            
            modify.set_leading_edge_slope_angle(slope_angle_new=slope_angle, width_bump=width_bump)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            slope_angle_new = geo_new.get_leading_edge_slope_angle()

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
            
            plt.title('Leading edge slope angle: %.1f -> %.1f'%(slope_angle_old, slope_angle_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-LE-slope-angle-%.1f.png'%(slope_angle_new)), dpi=300)
            
        
    #* Trailing edge wedge angle modification
    if False:

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
    if False:

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
            
    
    #* Change airfoil thickness at x=0.2, 0.7
    if False:
        
        for x_t, t_new in [(0.2, 0.085), (0.7, 0.055)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)
            t_old = geo_old.get_thickness_at(x_t)
            
            modify = FoilModification(x, yu, yl)
            
            modify.set_thickness_at(x_t, t_new, width_bump=0.6)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            t_new = geo_new.get_thickness_at(x_t)

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
            
            plt.title('Thickness at x=%.2f: %.3f -> %.3f'%(x_t, t_old, t_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-t-x%.2f.png'%(x_t)), dpi=300)

            plt.close()
            
    
    #* Change airfoil maximum thickness location
    if False:
        
        for x_t_new, slope1 in [(0.3, 1.0), (0.5, 1.0)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)
            t_old, x_t_old, _ = geo_old.get_maximum_thickness()

            modify = FoilModification(x, yu, yl)
            
            modify.set_maximum_thickness_location(x_t_new=x_t_new, slope0=1.0, slope1=slope1)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            t_new, x_t_new, _ = geo_new.get_maximum_thickness()

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
            
            plt.title('Location of t-max: %.1f -> %.1f'%(x_t_old, x_t_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-x_tmax-%.1f.png'%(x_t_new)), dpi=300)

            plt.close()
            
    
    #* Change airfoil camber
    if False:
        
        for xc, h, w in [(0.3, 0.01, 1.0), (0.8, 0.005, 0.6)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)

            modify = FoilModification(x, yu, yl)
            
            modify.add_bump_to_camber(xc, h, w, kind='H')
            
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
    
    
    #* Change airfoil camber at front and rear part
    if False:
        
        for side, c_new in [('front', 0.004), ('rear', 0.005)]:

            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)

            modify = FoilModification(x, yu, yl)
            
            if side == 'front':
                
                modify.set_camber_front(c_new)
                c_old = geo_old.get_average_camber_front_60p()
                
                geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
                c_new = geo_new.get_average_camber_front_60p()
                
            else:
                
                modify.set_camber_rear(c_new)
                c_old = geo_old.get_average_camber_rear_40p()
                
                geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
                c_new = geo_new.get_average_camber_rear_40p()

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

            if side == 'front':
                plt.title('Average camber of front part: %.3f -> %.3f'%(c_old, c_new))
            else:
                plt.title('Average camber of rear part: %.3f -> %.3f'%(c_old, c_new))

            plt.savefig(os.path.join(path, 'airfoil-modify-camber-%s.png'%(side)), dpi=300)
    
    
    #* Add bump to airfoil surfaces
    if False:
    
        for keep_tmax in [True, False]:
    
            x, yu, yl, t0, rLE_old = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.0)
            
            geo_old = FoilGeoFeatures(x, yu, yl)

            modify = FoilModification(x, yu, yl)
            
            modify.add_bump(bumps=[(0.3, 0.01, 1.0, 'upper', 'H'), (0.85, 0.008, 0.6, 'lower', 'H')], keep_tmax=keep_tmax)
            
            geo_new = FoilGeoFeatures(modify.x, modify.yu, modify.yl)
            
            plt.figure(figsize=(16,8))
            
            plt.plot(x, yu, 'k')
            plt.plot(x, yl, 'k', label='_nolegend_')
            
            plt.plot(x, geo_old.get_feature('thickness'), 'k--', lw=0.5)
            plt.plot(x, geo_old.get_feature('camber'), 'b--', lw=0.5)
            
            plt.plot(modify.x, modify.yu, 'r', lw=1.0)
            plt.plot(modify.x, modify.yl, 'r', lw=1.0, label='_nolegend_')
            
            plt.plot(modify.x, geo_new.get_feature('thickness'), 'r--', lw=0.5)
            plt.plot(modify.x, geo_new.get_feature('camber'), 'g--', lw=0.5)

            plt.xlim((-0.2, 1.2))
            plt.ylim((-0.2, 0.2))
            plt.axis('equal')
            plt.legend(['Airfoil (old)', 'Thickness (old)', 'Camber (old)', 'Airfoil (new)', 'Thickness (new)', 'Camber (new)'])

            if keep_tmax:
                plt.title('Add bump to airfoil surfaces (keep tmax)')
                plt.savefig(os.path.join(path, 'airfoil-modify-bump-keep-tmax.png'), dpi=300)
            else:
                plt.title('Add bump to airfoil surfaces')
                plt.savefig(os.path.join(path, 'airfoil-modify-bump.png'), dpi=300)
    
    