'''
Axisymmetric Flow Through Nacelle (FTN)
'''
import os
import sys
sys.path.append('.')

import numpy as np

from cst_modeling.io import output_surface
from cst_modeling.section import cst_foil
from cst_modeling.operation import Lofting_Revolution


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])


    #* Build an airfoil
    
    cst = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])

    x_, yu, yl, _, _ = cst_foil(201, cst, -cst, x=None, t=0.06, tail=0.004)

    xx = np.concatenate((np.flip(x_),x_[1:]), axis=0)
    yy = np.concatenate((np.flip(yl),yu[1:]), axis=0)


    #* Surface of revolution

    section_s_loc=[0.0, 0.25, 0.50, 0.75]

    profiles = []
    for i_profile in range(len(section_s_loc)):
        profiles.append([xx, yy])

    nacelle = Lofting_Revolution(
        profiles=profiles,
        section_s_loc=section_s_loc,
        section_x=0.0,
        section_radius=0.25,
        section_scale=1.0,
        n_spanwise=21,
    )

    surfs = nacelle.sweep(interp_profile_kind='linear')

    for i_surf, surf in enumerate(surfs):
        output_surface(surf, fname=os.path.join(path, 'nacelle-axisymmetric.dat'), ID=i_surf)

    nacelle.guide_curve.output(os.path.join(path, 'nacelle-axisymmetric-guide-curve.dat'))

