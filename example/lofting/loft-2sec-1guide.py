'''
Lofting a surface with 2 sections and 1 guide curve.
'''
import os
import sys
sys.path.append('.')

import numpy as np

from cst_modeling.operation import GuideCurve, Lofting_2Profile
from cst_modeling.basic import BasicSection
from cst_modeling.io import output_surface
from cst_modeling.section import cst_foil


N_SPANWISE = 101
N_POINT = 201
N_SECTION = 2


if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])
    
    fname_save1 = os.path.join(path, 'guide-2sec-1guide.dat')
    fname_save2 = os.path.join(path, 'surface-2sec-1guide.dat')

    #* Guide curve
    guide = GuideCurve(N_SECTION, N_SPANWISE, section_s_loc=[0.0, 1.0])

    control_points = {
        's':        [ 0.0, 1.0],
        'x':        [ 0.0, 0.2],
        'y':        [ 0.0, 0.1],
        'z':        [ 0.0, 1.0],
        'scale':    [ 1.0, 0.5],
        'rot_x':    [ 0.0, 0.0],
        'rot_y':    [ 0.0, 0.0],
        'rot_z':    [-2.0, 5.0],
    }

    for key in ['x', 'y', 'z', 'scale', 'rot_x', 'rot_y', 'rot_z']:
        guide.generate_by_interp1d(control_points['s'], control_points[key], key=key, kind='linear')

    guide.output(fname_save1)
    

    #* Build airfoil sections 
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    sec0 = BasicSection()
    sec1 = BasicSection()
    
    sec0.xx, sec0.yu, sec0.yl, _, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.12, tail=0.00)
    sec1.xx, sec1.yu, sec1.yl, _, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.08, tail=0.00)

    profile_0 = sec0.get_profile()
    profile_1 = sec1.get_profile()


    #* Lofting
    
    loft = Lofting_2Profile(profile_0, profile_1, n_spanwise=101, is_guide_curve_at_LE=True)
    
    loft.update_guide_curve(**guide.global_guide_curve)

    surf_x, surf_y, surf_z = loft.sweep()

    output_surface([surf_x, surf_y, surf_z], fname_save2, ID=0, zone_name='rotate-about-LE')


    loft = Lofting_2Profile(profile_0, profile_1, n_spanwise=101, is_guide_curve_at_LE=False)
    
    loft.update_guide_curve(**guide.global_guide_curve)
    
    #* Update the guide curve to go through the trailing edge
    loft.guide_curve['x'] = loft.guide_curve['x'] + loft.guide_curve['scale']

    surf_x, surf_y, surf_z = loft.sweep()

    output_surface([surf_x, surf_y, surf_z], fname_save2, ID=1, zone_name='rotate-about-TE')


