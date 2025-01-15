'''
NASA Common Research Model (CRM) Wing with a user-defined winglet.

https://commonresearchmodel.larc.nasa.gov/
'''

import os
import sys
sys.path.append('.')

import copy
from scipy.interpolate import interp1d
from cst_modeling.basic import BasicSection
from cst_modeling.surface2 import Surface


def extract_partial_section(sec: BasicSection, ratio: float, interp1d_kind='cubic') -> BasicSection:
    '''
    Extract part of the section, where xx in [0, ratio].
    '''
    sec = copy.deepcopy(sec)
    
    f_yu = interp1d(sec.xx, sec.yu, kind=interp1d_kind)
    f_yl = interp1d(sec.xx, sec.yl, kind=interp1d_kind)
    
    sec.yu = f_yu(sec.xx * ratio) / ratio
    sec.yl = f_yl(sec.xx * ratio) / ratio
    sec.scale = sec.scale * ratio
    
    return sec



if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])
    
    fname_setting = os.path.join(path, 'separate-winglet.txt')
    fname_tecplot = os.path.join(path, 'separate-winglet.dat')
    
    #* Main wing

    wing = Surface(n_sec=2, name='main_wing', nn=201, ns=51)
    wing.read_setting(fname_setting, tail=0.001)
    wing.geo()
    
    #* Wing let

    winglet = Surface(n_sec=2, name='winglet', nn=201, ns=51, rotate_x_section=True)
    winglet.read_setting(fname_setting, tail=0.001)
    winglet.geo()
    
    #* Transition from wing to winglet
    transition = Surface(n_sec=2, name='transition', nn=201, ns=501)
    
    # Set the two sections of the transition surface
    ratio = (winglet.sections[0].scale + winglet.sections[0].xLE - wing.sections[-1].xLE) / wing.sections[-1].scale
    transition.sections[0] = extract_partial_section(wing.sections[-1], ratio)
    transition.sections[1] = copy.deepcopy(winglet.sections[0])
    
    # Set the guide curve for the transition surface
    guide = transition.get_default_guide_curve()
    span  = (winglet.sections[0].zLE - wing.sections[-1].zLE)
    span_0 = (wing.sections[-1].zLE - wing.sections[0].zLE)
    span_1 = (winglet.sections[-1].zLE - winglet.sections[0].zLE)
    
    x0 = wing.sections[-1].xLE
    x1 = winglet.sections[0].xLE
    guide.generate_by_spline(global_control_s=[0.0, 1.0], key='x',
        global_values= [x0, x1],
        slope_s0= (x0 - wing.sections[0].xLE) / span_0 * span, 
        slope_s1= (winglet.sections[-1].xLE - x1) / span_1 * span, 
        )

    y0 = wing.sections[-1].yLE
    y1 = winglet.sections[0].yLE
    guide.generate_by_spline(global_control_s=[0.0, 1.0], key='y',
        global_values= [y0, y1],
        slope_s0= (y0 - wing.sections[0].yLE) / span_0 * span, 
        slope_s1= (winglet.sections[-1].yLE - y1) / span_1 * span,
        )
    
    for ii in range(transition.ns):
        value = winglet.sections[0].xLE + winglet.sections[0].scale - guide.get_value('x', ii)
        guide.set_value('scale', ii, value)
    
    guide.generate_rotation_angle_with_tangent('rot_x')
    
    
    # Generate the transition surface
    transition.prepare(guide=guide, update_section_profile=False)
    transition.geo()

    #* Output geometry
    
    wing.output_tecplot(fname=fname_tecplot, one_piece=False, split=True)
    winglet.output_tecplot(fname=fname_tecplot, one_piece=False, split=True, append=True)
    transition.output_tecplot(fname=fname_tecplot, one_piece=False, split=True, append=True)

    transition.output_guide_curve(os.path.join(path, 'guide-curve.dat'))
