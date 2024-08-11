import os
import sys
sys.path.append('.')

from cst_modeling.surface2 import Surface


if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])


    propeller = Surface(n_sec=11, name='Propeller', nn=101, ns=51,
                    smooth_surface=True, smooth_sections=None,
                    rotate_x_section=False, rotation_sections=None)

    propeller.read_setting(os.path.join(path, 'Propeller.txt'), tail=0.0)

    propeller.prepare()
    
    
    #* Update 'zLE' distribution by spline interpolation
    propeller.lofting.guide_curve.generate_by_spline(
        propeller.section_s_loc, propeller.spanwise_locations, key='z')
    
    #* Update 'rot_axis' distribution to align the tangent of the guide curve
    propeller.lofting.guide_curve.generate_rotation_angle_with_tangent(key='rot_axis')
    
    
    propeller.geo()

    propeller.output_tecplot(fname=os.path.join(path, 'toroidal-propeller.dat'), one_piece=False, split=False)

    propeller.output_guide_curve(fname=os.path.join(path, 'toroidal-propeller-guide-curve.dat'))
    
