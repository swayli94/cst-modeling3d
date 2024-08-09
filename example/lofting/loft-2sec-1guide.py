'''
Lofting a surface with 2 sections and 1 guide curve.
'''
import os
import sys
sys.path.append('.')

import numpy as np

from cst_modeling.operation import Lofting_2Profile
from cst_modeling.basic import BasicSection
from cst_modeling.io import output_surface
from cst_modeling.section import cst_foil


N_POINT = 201


if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])
    
    fname_save = os.path.join(path, 'surface-2sec-1guide.dat')


    #* Build airfoil sections 
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

    sec0 = BasicSection(chord=1.0)
    sec1 = BasicSection(chord=0.5)
    
    sec0.xx, sec0.yu, sec0.yl, sec0.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.12, tail=0.00)
    sec1.xx, sec1.yu, sec1.yl, sec1.thick, _ = cst_foil(N_POINT, cst_u, cst_l, t=0.08, tail=0.00)

    sec0.rot_z = -2.0

    sec1.xLE = 0.2
    sec1.yLE = 0.1
    sec1.zLE = 1.0

    sec1.rot_z = 5.0
    sec1.rot_x = 0.0
    sec1.rot_y = 0.0


    #* Lofting
    
    loft = Lofting_2Profile(sec0, sec1, n_spanwise=101, is_guide_curve_at_LE=True)

    surf_x, surf_y, surf_z = loft.sweep()

    output_surface([surf_x, surf_y, surf_z], fname_save, ID=0, zone_name='rotate-about-LE')


    loft = Lofting_2Profile(sec0, sec1, n_spanwise=101, is_guide_curve_at_LE=False)

    surf_x, surf_y, surf_z = loft.sweep()

    output_surface([surf_x, surf_y, surf_z], fname_save, ID=1, zone_name='rotate-about-TE')


