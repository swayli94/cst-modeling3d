'''
Lofting a surface with multiple sections and 1 guide curve.
'''
import os
import sys
sys.path.append('.')

import numpy as np

from cst_modeling.operation import GuideCurve, Lofting
from cst_modeling.basic import BasicSection
from cst_modeling.io import output_surface
from cst_modeling.section import cst_foil


N_SPANWISE = 101
N_POINT = 201
N_SECTION = 5

if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])
    
    fname_save1 = os.path.join(path, 'guide-multi-sec-toroidal.dat')
    fname_save2 = os.path.join(path, 'surface-multi-sec-toroidal.dat')


    #* Guide curve
    guide = GuideCurve(N_SECTION, N_SPANWISE, section_s_loc=[0.00, 0.25, 0.50, 0.75, 1.00])

    control_points = {
        's':        [ 0.00, 0.25, 0.50, 0.75, 1.00],
        'x':        [ 0.0,  0.0,  0.2,  0.0,  0.0 ],
        'y':        [ 0.0,  0.5,  1.0,  0.5,  0.0 ],
        'z':        [ 0.0,  1.0,  0.0, -1.0,  0.0 ],
        'scale':    [ 1.0,  0.8,  0.5,  0.8,  1.0 ],
        'rot_x':    [ 0.0,  0.0,  0.0,  0.0,  0.0 ],
        'rot_y':    [ 0.0,  0.0,  0.0,  0.0,  0.0 ],
        'rot_z':    [ 0.0,  0.0, 20.0,  0.0,  0.0 ],
    }

    for key in ['x','y', 'z', 'scale', 'rot_x', 'rot_y', 'rot_z']:
        guide.generate_by_spline(control_points['s'], control_points[key], key=key, slope_s0=None, slope_s1=None, periodic=True)


    
    guide.generate_by_spline(control_points['s'], control_points['rot_z'], key='rot_z', slope_s0=None, slope_s1=None, periodic=True)
    
    guide.update_rotation_angle_with_tangent(key='rot_x')

    guide.output(fname_save1)


    #* Build airfoil sections 
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    thickness = [0.12, 0.10, 0.05, 0.10, 0.12]

    sections = []
    for i_sec in range(N_SECTION):
        
        sec = BasicSection()
        sec.xx, sec.yu, sec.yl, sec.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=thickness[i_sec], tail=0.00)
        
        sections.append(sec)


    #* Lofting
    
    loft = Lofting(sections, guide, is_guide_curve_at_LE=True)

    loft.sweep(interp_profile_kind='quadratic')

    for i_surf in range(loft.n_section-1):
        
        output_surface(loft.surfs[i_surf], fname_save2, ID=i_surf)

