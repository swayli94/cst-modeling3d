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
N_SECTION = 3

if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])
    
    fname_save1 = os.path.join(path, 'guide-multi-sec-1guide.dat')
    fname_save2 = os.path.join(path, 'surface-multi-sec-1guide.dat')


    #* Guide curve
    guide = GuideCurve(N_SECTION, N_SPANWISE, section_s_loc=[0.0, 0.5, 1.0])

    control_points = {
        's':        [ 0.0, 0.5,  1.0],
        'x':        [ 0.0, 0.1,  0.2],
        'y':        [ 0.0, 0.0,  0.5],
        'z':        [ 0.0, 0.5,  1.0],
        'scale':    [ 1.0, 0.75, 0.5],
        'rot_x':    [ 0.0, 0.0,  0.0],
        'rot_y':    [ 0.0, 0.0,  0.0],
        'rot_z':    [-2.0, 0.0,  5.0],
    }

    guide.generate_by_interp1d(control_points['s'], control_points['x'], key='x', kind='linear')

    guide.generate_by_spline(control_points['s'], control_points['y'], key='y', slope_s0=0.0, slope_s1=None, periodic=False)
    
    guide.generate_by_interp1d(control_points['s'], control_points['z'], key='z', kind='linear')
    
    guide.generate_by_interp1d(control_points['s'], control_points['scale'], key='scale', kind='linear')
    guide.generate_by_interp1d(control_points['s'], control_points['rot_x'], key='rot_x', kind='linear')
    guide.generate_by_interp1d(control_points['s'], control_points['rot_y'], key='rot_y', kind='linear')
    
    guide.generate_by_spline(control_points['s'], control_points['rot_z'], key='rot_z', slope_s0=0.0, slope_s1=None, periodic=False)
    
    guide.update_rotation_angle_with_tangent(key='rot_x')

    guide.output(fname_save1)


    #* Build airfoil sections 
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    sec0 = BasicSection(chord=1.0)
    sec1 = BasicSection(chord=0.8)
    sec2 = BasicSection(chord=0.6)

    sec0.xx, sec0.yu, sec0.yl, sec0.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.12, tail=0.00)
    sec1.xx, sec1.yu, sec1.yl, sec1.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.09, tail=0.00)
    sec2.xx, sec2.yu, sec2.yl, sec2.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.08, tail=0.00)


    sections = [sec0, sec1, sec2]


    #* Lofting
    
    loft = Lofting(sections, guide, is_guide_curve_at_LE=True)

    loft.sweep(interp_profile_kind='quadratic')

    for i_surf in range(loft.n_section-1):
        
        output_surface(loft.surfs[i_surf], fname_save2, ID=i_surf)

