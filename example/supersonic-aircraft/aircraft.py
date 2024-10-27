'''
Lockheed Martin 1021 test case (AIAA Sonic Boom Workshop)

Note: this is only a simplified version of the original geometry.
The planform is not accurate, and geometry is not a perfect match.

https://lbpw.larc.nasa.gov/
'''

import os
import sys
sys.path.append('.')

import numpy as np
from cst_modeling.operation import Lofting, GuideCurve
from cst_modeling.surface2 import Surface, BasicSurface


class Fuselage(BasicSurface):
    
    def __init__(self, n_sec=1, name='Surf', nn=1001, ns=101,
            smooth_surface=False, smooth_sections = None,
            rotate_x_section=False, rotation_sections = None, 
            is_guide_curve_at_LE=True):
        
        super().__init__(n_sec, name, nn, ns, smooth_surface, smooth_sections, 
                    rotate_x_section, rotation_sections, is_guide_curve_at_LE)

    @staticmethod
    def create_circle(n, r, x0, y0):
        '''
        Create a circle in x-y plane
        '''
        theta = np.linspace(0, 2*np.pi, n)
        x = r * np.cos(theta) + x0
        y = r * np.sin(theta) + y0

        return x, y

    def update_section(self) -> None:
        '''
        Update all sections' profile curves, i.e., a closed circle
        '''
        for i in range(self.n_section):
            
            theta = np.linspace(0, 2*np.pi, self.nn, endpoint=True)
            
            self.sections[i].xx = np.cos(theta)
            self.sections[i].yy = np.sin(theta)

            self.sections[i].x = np.ones(self.nn) * self.sections[i].xLE
            self.sections[i].y = np.cos(theta) * self.sections[i].chord / 2 + self.sections[i].yLE
            self.sections[i].z = np.sin(theta) * self.sections[i].chord / 2 + self.sections[i].zLE

    def read_setting(self, fname) -> None:
        '''
        Read in Surface layout and CST parameters from file

        Parameters
        ----------
        fname : str
            settings file name.
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1}

        found_surf = False
        found_key = 0
        with open(fname, 'r') as f:

            lines = f.readlines()
            i_line = 0

            while i_line<len(lines):

                line = lines[i_line].split()

                if len(line) < 1:
                    i_line += 1
                    continue
                
                if not found_surf and len(line) > 1:
                    if '[Surf]' in line[0] and self.name == line[1]:
                        found_surf = True

                elif found_surf and '[Surf]' in line[0]:
                    break

                elif found_surf and found_key == 0:
                    if line[0] in key_dict:
                        found_key = key_dict[line[0]]

                elif found_surf and found_key == 1:
                    for i in range(self.n_section):
                        i_line += 1
                        line = lines[i_line].split()
                        self.sections[i].xLE   = float(line[0])
                        self.sections[i].yLE   = float(line[1])
                        self.sections[i].zLE   = float(line[2])
                        self.sections[i].chord = float(line[3])

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                i_line += 1
        
        print('Read surface [%s] settings'%(self.name))

        # Locate layout center for plot
        self.update_layout_center()

    def get_default_guide_curve(self) -> GuideCurve:
        '''
        Initialize the default guide curve object.
        It has a piecewise linear distribution along the span, defined by the section parameters.
        
        Returns
        --------
        guide: GuideCurve
            default guide curve object.
        '''
        #* Calculate parametric coordinates for the sections
        
        section_s_loc = [0.0]
        
        for i in range(self.n_section-1):
            
            section_s_loc.append(section_s_loc[i] + np.sqrt((self.sections[i+1].xLE-self.sections[i].xLE)**2 + 
                                                            (self.sections[i+1].yLE-self.sections[i].yLE)**2 + 
                                                            (self.sections[i+1].zLE-self.sections[i].zLE)**2) )

        ds = section_s_loc[-1]
        
        for i in range(self.n_section):
            section_s_loc[i] = (section_s_loc[i])/ds
            
        self.section_s_loc = section_s_loc
        
        #* Setup the control points
        
        control_points = {
            'x':        [sec.xLE    for sec in self.sections],
            'y':        [sec.yLE    for sec in self.sections],
            'z':        [sec.zLE    for sec in self.sections],
            'scale':    [sec.chord  for sec in self.sections],
            'rot_x':    [0.0        for _   in self.sections],
            'rot_y':    [90.0       for _   in self.sections],
            'rot_z':    [0.0        for _   in self.sections],
            'rot_axis': [0.0        for _   in self.sections],
        }
                   

        #* Generate the default (piecewise linear) guide curve
        
        guide = GuideCurve(self.n_section, n_spanwise=self.ns, section_s_loc=section_s_loc)

        for key, value in control_points.items():
            guide.generate_by_interp1d(section_s_loc, value, key=key, kind='linear')
        
        return guide



    def prepare(self):
        
        #* Update section profile
        for i in range(self.n_section):
            
            theta = np.linspace(0, 2*np.pi, self.nn, endpoint=True)
            
            self.sections[i].xx = np.cos(theta)
            self.sections[i].yy = np.sin(theta)
        
        #* Define guide curve
        guide = self.get_default_guide_curve()

        profiles = self.get_profiles()
        
        self.lofting = Lofting(profiles, global_guide_curve=guide, 
                                is_guide_curve_at_LE=self.is_guide_curve_at_LE)


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])


    #* Delta wing
    if True:
        
        wing = Surface(n_sec=4, name='wing', nn=201, ns=21)

        wing.read_setting(fname=os.path.join(path, 'aircraft.txt'), tail=0.005)
        wing.prepare()
        wing.geo()

        wing.output_tecplot(fname=os.path.join(path, 'wing-1.dat'), one_piece=True, split=False)
        wing.output_plot3d(fname=os.path.join(path, 'wing-1.grd'), split=False)
        
        wing.flip(plane='XY')
        wing.output_tecplot(fname=os.path.join(path, 'wing-2.dat'), one_piece=True, split=False)
        wing.output_plot3d(fname=os.path.join(path, 'wing-2.grd'), split=False)

    #* V-tail
    if True:
        
        vTail = Surface(n_sec=2, name='vTail', nn=201, ns=21)

        vTail.read_setting(fname=os.path.join(path, 'aircraft.txt'), tail=0.005)
        vTail.prepare()
        vTail.geo()

        vTail.output_tecplot(fname=os.path.join(path, 'vTail-1.dat'), one_piece=True, split=False)
        vTail.output_plot3d(fname=os.path.join(path, 'vTail-1.grd'), split=False)
        
        vTail.flip(plane='XY')
        vTail.output_tecplot(fname=os.path.join(path, 'vTail-2.dat'), one_piece=True, split=False)
        vTail.output_plot3d(fname=os.path.join(path, 'vTail-2.grd'), split=False)

    #* Fuselage
    if True:

        fuselage = Fuselage(n_sec=4, name='fuselage', nn=201, ns=21)
        
        fuselage.read_setting(fname=os.path.join(path, 'aircraft.txt'))
        fuselage.prepare()
        fuselage.geo()

        fuselage.output_tecplot(fname=os.path.join(path, 'fuselage.dat'), one_piece=True)
        fuselage.output_plot3d(fname=os.path.join(path, 'fuselage.grd'))

