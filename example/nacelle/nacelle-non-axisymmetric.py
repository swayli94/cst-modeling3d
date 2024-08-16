'''
Non-axisymmetric Powered Engine Nacelle (PEN)
'''
import os
import sys
sys.path.append('.')

import numpy as np
import matplotlib.pyplot as plt

from cst_modeling.io import output_surface
from cst_modeling.operation import Lofting_Revolution
from cst_modeling.tools.nacelle import NacelleIntakeHighlight, PoweredNacelleProfile


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])
    
    circum_control_psi = [0.0, 90.0, 180.0, 270.0]


    #* Nacelle Intake Highlight
    
    nacelle_highlight = NacelleIntakeHighlight(
        l_intake=35.0, theta_droop=0.0, theta_scarf=0.0,
        n_circum=101, circum_control_psi=circum_control_psi,
        circum_control_r_highlight=39.0
    )
    
    curve, psi_curve = nacelle_highlight.calculate()
    
    with open(os.path.join(path, 'nacelle-highlight.dat'), 'w') as f:
        
        f.write('Variables= X Y Z psi \n')
        
        f.write('zone T="Nacelle Highlight" I=%d\n' % len(curve))
        for i in range(len(curve)):
            f.write('%f %f %f %f\n' % (curve[i][0], curve[i][1], curve[i][2], 0.0))
    
        f.write('zone T="Nacelle Highlight (psi)" I=%d\n' % len(curve))
        for i in range(len(curve)):
            f.write('%f %f %f %f\n' % (psi_curve[i][1], psi_curve[i][2], psi_curve[i][3], psi_curve[i][0]))
    
    
    #* Nacelle Profile
    
    profiles = []
    
    for i_sec in range(len(circum_control_psi)):
        
        psi = circum_control_psi[i_sec]
    
        nacelle_profile = PoweredNacelleProfile(psi=psi, n_point_segment=101)
        
        highlight = nacelle_highlight.get_coordinate_2d(nacelle_profile.psi)
        
        nacelle_profile.set_parameters(
            r_spinner=9.0, theta_spinner=36.0, r_fan=None, 
            highlight_x= highlight[0], highlight_y= highlight[1],
            intake_face_center= nacelle_highlight.intake_face_center, 
            l_nacelle=135.0, r_te=33.0,
            l_fan=35.0, r_bypass_outer=None, r_bypass_inner=19.0,
            x_core_cowl_0=95.0, y_core_cowl_0=22.0,
            x_core_cowl_1=125.0, y_core_cowl_1=14.0,
            x_core_duct=109.0, r_core_outer=15.0, r_core_inner=7.5,
            x_core_plug_0=123.0, y_core_plug_0=8.2, x_core_plug_1=150.0,
            cst_u=[ 0.10, 0.10, 0.10, 0.10, 0.10], 
            cst_l=[-0.10, 0.15,-0.10, 0.05, 0.05],
            bypass_inner_angle=0.0,
            bypass_inner_control_points=[],
            core_outer_control_points=[],
            core_inner_control_points=[],
        )
        
        profile_x, profile_y = nacelle_profile.get_profile()
        
        profiles.append([profile_x, profile_y])
        
        #nacelle_profile.plot(show=True)
    
    
    #* Nacelle surface of revolution
    
    nacelle = Lofting_Revolution(
        profiles=profiles,
        section_s_loc=[a/360.0 for a in circum_control_psi],
        section_x=0.0,
        section_radius=0.0,
        section_scale=1.0,
        n_spanwise=51,
    )

    surfs = nacelle.sweep(interp_profile_kind='periodic')

    for i_surf, surf in enumerate(surfs):
        output_surface(surf, fname=os.path.join(path, 'nacelle-non-axisymmetric.dat'), ID=i_surf)

