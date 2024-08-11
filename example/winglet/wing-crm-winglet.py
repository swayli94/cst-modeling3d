'''
NASA Common Research Model (CRM) Wing with a user-defined winglet.

https://commonresearchmodel.larc.nasa.gov/
'''

import os
import sys
sys.path.append('.')

from cst_modeling.surface2 import Surface


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])

    '''
    For winglets, rotating the wing section about x-axis is only necessary when the winglet is close to vertical.
    
    The rotation sections can be just winglet sections, because usually we do not want to rotate the wing sections.
    '''
    
    wing = Surface(n_sec=10, name='Wing-CRM-winglet', nn=201, ns=51, 
                    smooth_surface=True, smooth_sections=[(0, 2), (4, 7), (8, 9)],
                    rotate_x_section=True, rotation_sections=[(8, 9)])

    wing.read_setting(os.path.join(path, 'Wing.txt'))

    wing.prepare()
    
    wing.geo()
    
    wing.flip(plane='XY')

    wing.output_tecplot(fname=os.path.join(path, 'Wing-CRM-winglet.dat'), one_piece=False, split=True)

    wing.output_guide_curve(os.path.join(path, 'Wing-CRM-winglet-guide-curve.dat'))

