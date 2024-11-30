'''
Extract airfoil geometric features
'''
import os
import sys
sys.path.append('.')

import numpy as np
from matplotlib import pyplot as plt

from cst_modeling.section import cst_foil
from cst_modeling.foil import FoilGeoFeatures


if __name__ == '__main__':

    path = os.path.dirname(sys.argv[0])

    #* Initialize an airfoil
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    for tail in [0.0, 0.004]:

        x, yu, yl, t0, rLE = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=tail)

        geo = FoilGeoFeatures(x, yu, yl)
        
        geo.get_maximum_thickness()
        geo.get_maximum_camber()
        geo.get_thickness_at(0.2)
        geo.get_thickness_at(0.7)
        geo.get_volume()
        geo.get_average_camber()
        geo.get_average_camber_front_60p()
        geo.get_average_camber_rear_40p()
        geo.get_curvature()
        geo.get_leading_edge_radius()
        geo.get_leading_edge_slope_angle()
        geo.get_trailing_edge_wedge_angle()
        geo.get_trailing_edge_slope_angle()
        geo.get_upper_crest_point()
        geo.get_lower_crest_point()
        
        
        fig, ax = plt.subplots(figsize=(16, 8))
        plt.plot(x, yu, 'k')
        plt.plot(x, yl, 'k', label='_nolegend_')
        
        plt.plot(x, geo.get_feature('thickness'), 'b--', lw=0.5)
        plt.plot(x, geo.get_feature('camber'), 'r--', lw=0.5)
        
        plt.xlim((-0.2, 1.2))
        plt.ylim((-0.2, 0.2))
        plt.axis('equal')
        plt.legend(['Airfoil surface', 'Thickness', 'Camber'])
        
        def add_label(x, y, label, dx=0.01, dy=0.01, color='b'):
            ax.plot(x, y, color+'*')
            ax.text(x+dx, y+dy, label, color=color)
        
        add_label(  x=geo.get_feature('x_t'), y=0.0, 
                    label=r'$t_{max}$: '+'%.3f'%geo.get_feature('t_max'))
        
        add_label(  x=0.2, y=0.0, 
                    label=r'$t_{0.2}$: '+'%.3f'%geo.get_feature('t_20'))
        
        add_label(  x=0.7, y=0.0, 
                    label=r'$t_{0.7}$: '+'%.3f'%geo.get_feature('t_70'))
        
        add_label(  x=geo.get_feature('x_c'), y=0.0, 
                    label=r'$c_{max}$: '+'%.3f'%geo.get_feature('c_max'),
                    dx=-0.02, dy=-0.04, color='r')
        
        add_label(  x=geo.get_feature('x_uc'), y=geo.get_feature('y_uc'), 
                    label=r'$y_{uc}$: '+'%.3f'%geo.get_feature('y_uc'))
        
        add_label(  x=geo.get_feature('x_lc'), y=geo.get_feature('y_lc'), 
                    label=r'$y_{lc}$: '+'%.3f'%geo.get_feature('y_lc'))
        
        add_label(  x=0.8, y=0.20, 
                    label=r'Volume: '+'%.3f'%geo.get_feature('volume'), dy=0)
        
        add_label(  x=0.8, y=0.18, 
                    label=r'Mean curvature: '+'%.3f'%geo.get_feature('c_mean'), dy=0)
        
        add_label(  x=0.3, y=0.0, 
                    label=r'$c_{front 0.6}$: '+'%.3f'%geo.get_feature('c_f60'),
                    dx=-0.02, dy=-0.02, color='r')
        
        add_label(  x=0.8, y=0.0, 
                    label=r'$c_{rear 0.4}$: '+'%.3f'%geo.get_feature('c_r40'),
                    dx=-0.02, dy=-0.02, color='r')
        
        add_label(  x=0.8, y=0.14, 
                    label=r'Leading edge radius: '+'%.3f'%geo.get_feature('r_le'), dy=0)
        
        add_label(  x=0.8, y=0.12, 
                    label=r'Leading edge slope angle: '+'%.3f'%geo.get_feature('slope_angle_le'), dy=0)
        
        add_label(  x=0.8, y=0.10, 
                    label=r'Trailing edge wedge angle: '+'%.3f'%geo.get_feature('wedge_angle'), dy=0)
        
        add_label(  x=0.8, y=0.08, 
                    label=r'Trailing edge slope angle: '+'%.3f'%geo.get_feature('slope_angle_te'), dy=0)
        

        plt.savefig(os.path.join(path, 'airfoil-geometric-features-tail-%.3f.png'%(tail)), dpi=300)
        plt.show()

