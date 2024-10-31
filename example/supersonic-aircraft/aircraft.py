'''
Lockheed Martin 1021 test case (AIAA Sonic Boom Workshop)

Note: this is only a simplified version of the original geometry.
The planform is not accurate, and geometry is not a perfect match.

https://lbpw.larc.nasa.gov/
'''

import os
import sys
sys.path.append('.')

import copy
import numpy as np
from cst_modeling.operation import GuideCurve
from cst_modeling.surface2 import Surface, BasicSurface
from cst_modeling.io import plot3d_to_igs, output_plot3d_for_parts


class Fuselage(BasicSurface):
    '''
    The major difference between Fuselage and Surface is that
    the 'span-wise' direction changes from the 'z' direction to the 'x' direction.
    
    In consequence, the 3D section profile is in y-z plane.
    '''
    
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
        Read in Surface layout from file

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

    def get_default_guide_curve(self, **kwargs) -> GuideCurve:
        '''
        Initialize the default guide curve object.
        It has a piecewise linear distribution along the span, defined by the section parameters.
        
        Returns
        --------
        guide: GuideCurve
            default guide curve object.
        '''
        custom_control_points = {
            'x':        [sec.xLE    for sec in self.sections],
            'y':        [sec.yLE    for sec in self.sections],
            'z':        [sec.zLE    for sec in self.sections],
            'scale':    [sec.chord  for sec in self.sections],
            'rot_x':    [0.0        for _   in self.sections],
            'rot_y':    [90.0       for _   in self.sections],
            'rot_z':    [0.0        for _   in self.sections],
            'rot_axis': [0.0        for _   in self.sections],
        }
                   
        guide = super().get_default_guide_curve(custom_control_points, smooth_keys=['x', 'y', 'z', 'scale'])
        
        return guide


if __name__ == "__main__":

    path = os.path.dirname(sys.argv[0])


    #* Geometry
    if True:
        
        #* Delta wing
        wing1 = Surface(n_sec=4, name='wing', nn=201, ns=101, smooth_surface=True)
        wing1.read_setting(fname=os.path.join(path, 'aircraft.txt'), tail=0.005)
        wing1.prepare(smooth_keys=['x', 'y', 'z', 'scale', 'rot_z'])
        wing1.geo()
        
        wing2 = copy.deepcopy(wing1)
        wing2.flip(plane='XY')
        

        #* V-tail
        vTail1 = Surface(n_sec=2, name='vTail', nn=201, ns=21)
        vTail1.read_setting(fname=os.path.join(path, 'aircraft.txt'), tail=0.005)
        vTail1.prepare()
        vTail1.geo()

        vTail2 = copy.deepcopy(vTail1)
        vTail2.flip(plane='XY')


        #* Fuselage
        fuselage = Fuselage(n_sec=6, name='fuselage', nn=201, ns=21, smooth_surface=True,
                    smooth_sections=[(0,2), (3,5)])
        fuselage.read_setting(fname=os.path.join(path, 'aircraft.txt'))
        fuselage.prepare()
        fuselage.geo()

        
    #* Output each part
    if False:
        
        wing1.output_plot3d(fname=os.path.join(path, 'wing-1.xyz'), one_piece=True, split=False)
        wing2.output_plot3d(fname=os.path.join(path, 'wing-2.xyz'), one_piece=True, split=False)
        
        vTail1.output_plot3d(fname=os.path.join(path, 'vTail-1.xyz'), one_piece=True, split=False)
        vTail2.output_plot3d(fname=os.path.join(path, 'vTail-2.xyz'), one_piece=True, split=False)
        
        fuselage.output_plot3d(fname=os.path.join(path, 'fuselage.xyz'), one_piece=True)
        
        plot3d_to_igs(fname=os.path.join(path, 'wing-1'))
        plot3d_to_igs(fname=os.path.join(path, 'wing-2'))
        plot3d_to_igs(fname=os.path.join(path, 'vTail-1'))
        plot3d_to_igs(fname=os.path.join(path, 'vTail-2'))
        plot3d_to_igs(fname=os.path.join(path, 'fuselage'))
        
        
    #* Output parts in one file
    if True:
        
        fname = os.path.join(path, 'aircraft')
        
        output_plot3d_for_parts(fname+'.xyz', 
                wing1.surfaces, wing2.surfaces, vTail1.surfaces, vTail2.surfaces, fuselage.surfaces)
        
        plot3d_to_igs(fname)
        
        
        
        