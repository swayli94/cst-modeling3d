'''
#TODO: Modify Surface Class to use Lofting



'''
import os
import copy
import numpy as np

from typing import List, Tuple

import matplotlib.pyplot as plt

from .basic import BasicSection
from .section import Section, OpenSection
from .operation import Lofting, GuideCurve
from .io import output_plot3d


class BasicSurface():
    '''
    Multi-section surface based on `BasicSection` and `Lofting`.
    
    Parameters
    ----------
    n_section : int
        number of control sections.
        
    name : str
        name of the object, by default 'Surf'.
        
    nn : int
        number of points in the unit 2D curve's `xx`, by default 1001.
        
    ns : int
        number of points in the sweep direction between sections, by default 101.
            
    smooth_surface : bool
        whether to smooth the surface, by default False.
        
    smooth_sections : List[Tuple[int, int]]
        surfaces to be smoothed, by default None. None means all surfaces are smoothed.
        The tuple is (start, end) index of the sections.
        For example, [(0, 1), (2, 4)] means the surfaces between the 0-1 and 2-4 sections are smoothed.
        
    rotate_x_section : bool
        whether to rotate the section about the x-axis, by default False.
        
    rotation_sections : List[Tuple[int, int]]
        sections to be rotated about the x-axis, by default None. None means all sections are rotated.
        The tuple is (start, end) index of the sections.
        For example, [(0, 1), (2, 4)] means the sections between the 0-1 and 2-4 sections are rotated.
    
    is_guide_curve_at_LE : bool
        whether the guide curve runs through the leading edge, by default True.
        Details can be found in the `Lofting` class.
    
    
    Notes
    -------
    - +x:     flow direction (m)
    - +y:     upside (m)
    - +z:     span-wise (m)

        
    Attributes
    ------------        
    sections : list of BasicSection
        section objects.
        
    surfaces : List[List[np.ndarray]]  [n_section-1][3][ns, nn]
        surface coordinates, i.e., a list of [surf_x, surf_y, surf_z].
        The surface coordinates `surf_*` are 2D arrays with shape (ns, nn).
        
    section_s_loc : List[float]
        span-wise parametric locations of the sections.
        
    half_size : float
        half of the size of the entire surface for plotting.
        
    center : ndarray
        surface center for plotting.
    '''
    def __init__(self, n_sec=1, name='Surf', nn=1001, ns=101, 
                    smooth_surface=False, smooth_sections : List[Tuple[int, int]] = None,
                    rotate_x_section=False, rotation_sections : List[Tuple[int, int]] = None,
                    is_guide_curve_at_LE=True):
        
        self.name = name
        self.nn = nn
        self.ns = ns
        self.smooth_surface = smooth_surface
        self.smooth_sections = smooth_sections
        self.rotate_x_section = rotate_x_section
        self.rotation_sections = rotation_sections
        self.is_guide_curve_at_LE = is_guide_curve_at_LE
        
        #* Attributes
        
        self.sections = [ BasicSection() for _ in range(max(1, n_sec)) ]
        self.surfaces : List[List[np.ndarray]] = []
        
        self.lofting : Lofting = None
        
        #* Parameters for plot
        self.half_size = 0.5
        self.center = np.array([0.5, 0.0, 0.5])
        
    @property
    def is_2d(self) -> bool:
        '''
        Whether this is a 3D surface for a 2D curve (unit span). 
        '''
        return self.n_section == 1
        
    @property
    def n_section(self) -> int:
        '''
        Number of sections
        '''
        return len(self.sections)
        
    @property
    def spanwise_locations(self) -> np.ndarray:
        '''
        Span-wise locations of the sections, i.e., `zLE` for each section.
        '''
        return np.array([self.sections[i].zLE for i in range(self.n_section)])
        
    def read_setting(self, fname: str) -> None:
        '''
        Read the surface layout parameters from file
        
        Parameters
        -----------
        fname : str
            parameter file name.
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface setting')
        
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
                        self.sections[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.sections[i].specified_thickness = float(line[5])

                        if self.is_2d:
                            self.sections[i].zLE = 0.0

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                i_line += 1
        
        self.update_layout_center()

    def update_layout_center(self) -> None:
        '''
        Update the layout center for plotting
        '''
        x_range = [self.sections[0].xLE, self.sections[0].xLE + self.sections[0].chord]
        y_range = [self.sections[0].yLE, self.sections[0].yLE]
        z_range = [self.sections[0].zLE, self.sections[0].zLE + self.is_2d]
        
        if not self.is_2d:
            
            for i in range(self.n_section):
                x_range[0] = min(x_range[0], self.sections[i].xLE)
                x_range[1] = max(x_range[1], self.sections[i].xLE + self.sections[i].chord)
                y_range[0] = min(y_range[0], self.sections[i].yLE)
                y_range[1] = max(y_range[1], self.sections[i].yLE)
                z_range[0] = min(z_range[0], self.sections[i].zLE)
                z_range[1] = max(z_range[1], self.sections[i].zLE)
        
        self.half_size = np.max([x_range[1]-x_range[0], y_range[1]-y_range[0], z_range[1]-z_range[0]])*0.5
        
        self.center[0] = 0.5*(x_range[1]+x_range[0])
        self.center[1] = 0.5*(y_range[1]+y_range[0])
        self.center[2] = 0.5*(z_range[1]+z_range[0])
        
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
            'rot_y':    [0.0        for _   in self.sections],
            'rot_z':    [sec.twist  for sec in self.sections],
            'rot_axis': [0.0        for _   in self.sections],
        }
                   

        #* Generate the default (piecewise linear) guide curve
        
        guide = GuideCurve(self.n_section, n_spanwise=self.ns, section_s_loc=section_s_loc)

        for key, value in control_points.items():
            guide.generate_by_interp1d(section_s_loc, value, key=key, kind='linear')
        
        
        #* Update the guide curve for smoothing
        
        if self.smooth_surface:
                            
            if self.smooth_sections is None:
                
                for key in ['x', 'y', 'scale', 'rot_z']:
                
                    guide.generate_by_spline(section_s_loc, control_points[key], key=key, slope_s0=None, slope_s1=None, periodic=False)
                
            else:
                
                for start, end in self.smooth_sections:
                    
                    control_s = section_s_loc[start:end+1]

                    for key in ['x', 'y', 'scale', 'rot_z']:
                        
                        control_v = control_points[key][start:end+1]
                    
                        if start == 0:
                            slope_s0 = None
                        else:
                            slope_s0 = (control_points[key][start]-control_points[key][start-1])/(section_s_loc[start]-section_s_loc[start-1])
                            
                        if end == self.n_section-1:
                            slope_s1 = None
                        else:
                            slope_s1 = (control_points[key][end+1]-control_points[key][end])/(section_s_loc[end+1]-section_s_loc[end])
                        
                        guide.update_by_spline(control_s, control_v, key=key, slope_s0=slope_s0, slope_s1=slope_s1, periodic=False)
                    
        
        #* Update the guide curve for rotating section about the x-axis
        
        if self.rotate_x_section:
            
            guide.update_rotation_angle_with_tangent(key='rot_x', sections=self.rotation_sections)
        
        return guide

    def get_profiles(self) -> List[List[np.ndarray]]:
        '''
        Get the 2D profiles for all sections.
        
        Returns
        --------
        profiles : List[List[np.ndarray]]
            2D profiles for each section, i.e., a list of [xx, yu, yl].
        '''
        profiles = []
        
        for i in range(self.n_section):
            profiles.append(self.sections[i].get_profile())
            
        return profiles

    def update_section(self) -> None:
        '''
        Update all sections' profile curves.
        '''
        raise NotImplementedError('This method should be implemented in the subclass.')

    def prepare(self, guide: GuideCurve = None, update_section_profile = True) -> None:
        '''
        Prepare the profiles, guide curve and lofting object for surface generation.
        
        Parameters
        -----------
        guide : GuideCurve
            user-defined GuideCurve object, by default None.
        
        update_section_profile : bool
            whether to update the section profile curves, by default True.
            
        Notes
        -------
        In the `BasicSurface` class, the `update_section` method is not implemented.
        Therefore, the profile curves are not generated.
        '''
        if update_section_profile:
            self.update_section()
        
        if guide is None:
            guide = self.get_default_guide_curve()
            
        profiles = self.get_profiles()
        
        self.lofting = Lofting(profiles, global_guide_curve=guide, 
                                is_guide_curve_at_LE=self.is_guide_curve_at_LE)

    def geo(self) -> None:
        '''
        Interpolate the 3D surface.
        
        Notes
        ------
        The `prepare()` method should be called first to prepare the profiles and lofting object.
        '''
        if self.lofting is None:
            #raise Exception('Please call `prepare()` first, to prepare the profiles and lofting object.')
            self.prepare()

        kind = None
        
        if self.smooth_surface:
            
            if self.n_section == 2:
                kind = 'linear'
                
            elif self.n_section == 3:
                kind = 'quadratic'
                
            elif self.n_section > 3:
                kind = 'cubic'
                
            else:
                raise Exception('The number of sections should be greater than 1.')

        self.surfaces = self.lofting.sweep(interp_profile_kind=kind)

    def flip(self, axis='None', plane='None') -> None:
        '''
        Flip surfaces and layout center. This should be the last action before output.
        
        Parameters
        -----------
        axis : str
            turn 90 degrees about axis: +X, -X, +Y, -Y, +Z, -Z.
            
        plane : str
            get symmetry about plane: 'XY', 'YZ', 'ZX'.

        Notes
        ---------
        The `axis` and `plane` can be a single phrase, 
        or a string contains multiple actions to take in order, e.g., '+X  +Y'.
        '''
        for axis_ in axis.split():
            if '+X' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][2]
                    self.surfaces[i_sec][2] = copy.deepcopy(self.surfaces[i_sec][1])
                    self.surfaces[i_sec][1] = copy.deepcopy(temp)

                temp = self.center[2]*1.0
                self.center[2] = self.center[1]*1.0
                self.center[1] = -temp

            if '-X' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][1]
                    self.surfaces[i_sec][1] = copy.deepcopy(self.surfaces[i_sec][2])
                    self.surfaces[i_sec][2] = copy.deepcopy(temp)

                temp = self.center[1]*1.0
                self.center[1] = self.center[2]
                self.center[2] = -temp

            if '+Y' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][0]
                    self.surfaces[i_sec][0] = copy.deepcopy(self.surfaces[i_sec][2])
                    self.surfaces[i_sec][2] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[2]
                self.center[2] = -temp

            if '-Y' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][2]
                    self.surfaces[i_sec][2] = copy.deepcopy(self.surfaces[i_sec][0])
                    self.surfaces[i_sec][0] = copy.deepcopy(temp)

                temp = self.center[2]
                self.center[2] = self.center[0]
                self.center[0] = -temp

            if '+Z' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][1]
                    self.surfaces[i_sec][1] = copy.deepcopy(self.surfaces[i_sec][0])
                    self.surfaces[i_sec][0] = copy.deepcopy(temp)

                temp = self.center[1]
                self.center[1] = self.center[0]
                self.center[0] = -temp

            if '-Z' in axis_:
                for i_sec in range(len(self.surfaces)):
                    temp = -self.surfaces[i_sec][0]
                    self.surfaces[i_sec][0] = copy.deepcopy(self.surfaces[i_sec][1])
                    self.surfaces[i_sec][1] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[1]
                self.center[1] = -temp

        if 'XY' in plane:
            for i_sec in range(len(self.surfaces)):
                self.surfaces[i_sec][2] = -self.surfaces[i_sec][2]
            self.center[2] = - self.center[2]

        if 'YZ' in plane:
            for i_sec in range(len(self.surfaces)):
                self.surfaces[i_sec][0] = -self.surfaces[i_sec][0]
            self.center[0] = - self.center[0]

        if 'ZX' in plane:
            for i_sec in range(len(self.surfaces)):
                self.surfaces[i_sec][1] = -self.surfaces[i_sec][1]
            self.center[1] = - self.center[1]

    def translate(self, dX=0.0, dY=0.0, dZ=0.0) -> None:
        '''
        Translate surface coordinates
        '''
        for surf in self.surfaces:
            surf[0] += dX
            surf[1] += dY
            surf[2] += dZ

        self.center[0] += dX
        self.center[1] += dY
        self.center[2] += dZ

    def split(self, index_splitting_point: List[int]) -> None:
        '''
        Split each surface into several pieces in the nn direction.
        
        Parameters
        -----------
        index_splitting_point : list of int
            index of the split points on each section curves
        '''
        nn = self.surfaces[0][0].shape[1]  # i.e., self.nn
        index_splitting_point.sort()
        index_splitting_point = [1] + index_splitting_point + [nn]
        
        surfaces = copy.deepcopy(self.surfaces)
        self.surfaces = []

        for i in range(len(index_splitting_point)-1):

            for surf in surfaces:

                surf_x = surf[0]
                surf_y = surf[1]
                surf_z = surf[2]
                
                self.surfaces.append([
                    surf_x[:,index_splitting_point[i]-1:index_splitting_point[i+1]],
                    surf_y[:,index_splitting_point[i]-1:index_splitting_point[i+1]],
                    surf_z[:,index_splitting_point[i]-1:index_splitting_point[i+1]]])

    #* Output
    def output_tecplot(self, fname=None, one_piece=False) -> None:
        '''
        Output the surface to `*.dat` in Tecplot format.
        
        Parameters
        ------------
        fname : str
            name of the output file.
            
        one_piece : bool
            if True, combine the span-wise sections into one piece.
        '''
        # surf_x[ns,nt], ns => spanwise

        if fname is None:
            fname = self.name + '.dat'

        n_sec   = 1 if self.is_2d else self.n_section-1
        n_piece = len(self.surfaces)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            ns = self.ns

            if not one_piece:

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfaces[i_sec][0]
                    surf_y = self.surfaces[i_sec][1]
                    surf_z = self.surfaces[i_sec][2]
                    
                    nt = surf_x.shape[1]

                    f.write('zone T="sec %d" i= %d j= %d \n'%(i_sec, nt, ns))

                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))
                            
            else:
                
                n_point = n_sec*(self.ns-1) + 1
                nt = self.surfaces[0][0].shape[1]
                
                f.write('zone T="sec" i= %d j= %d \n'%(nt, n_point))

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfaces[i_sec][0]
                    surf_y = self.surfaces[i_sec][1]
                    surf_z = self.surfaces[i_sec][2]
                    
                    nt = surf_x.shape[1]

                    if i_sec>=n_piece-1:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))

    def output_plot3d(self, fname=None, scale=1.0) -> None:
        '''
        Output the surface to `*.grd` in plot3d format.

        Parameters
        ------------
        fname : str
            name of the output file.
            
        scale : float
            scaling factor for the coordinates
        '''
        Xs = [self.surfaces[i][0] for i in range(len(self.surfaces))]
        Ys = [self.surfaces[i][1] for i in range(len(self.surfaces))]
        Zs = [self.surfaces[i][2] for i in range(len(self.surfaces))]
        
        if fname is None:
            fname = self.name + '.grd'
        
        output_plot3d(Xs, Ys, Zs, fname=fname, scale=scale)

    def output_guide_curve(self, fname=None) -> None:
        '''
        Output the guide curve to `*.dat` in Tecplot format.

        Parameters
        ------------
        fname : str
            name of the output file.
        '''
        if fname is None:
            fname = self.name + '-guide.dat'

        self.lofting.guide_curve.output(fname)

    def plot(self, fig_id=1, type='wireframe', show=True):
        '''
        Plot surface (the figure is not closed).

        Parameters
        ------------
        fig_id : int
            ID of the figure
        type : str
            'wireframe', or 'surface'
        show : bool
            whether plot on screen
            
        Return
        ---------
        ax
            figure (subplot) handle
        '''
        self.update_layout_center()
        
        fig = plt.figure(fig_id)
        ax = fig.add_subplot(projection='3d')

        for surf in self.surfaces:
            if type in 'wireframe':
                ax.plot_wireframe(surf[0], surf[1], surf[2])
            else:
                ax.plot_surface(surf[0], surf[1], surf[2])

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(self.center[0]-self.half_size, self.center[0]+self.half_size)
        ax.set_ylim3d(self.center[1]-self.half_size, self.center[1]+self.half_size)
        ax.set_zlim3d(self.center[2]-self.half_size, self.center[2]+self.half_size)

        if show:
            plt.show()
            
        return ax

    #* Obsolete functions
    def add_sec(self, *args, **kwargs) -> None:
        '''
        Obsolete function.
        '''
        print('>>> [Error] The `add_sec` method is not implemented.')
        raise NotImplementedError
    
    def smooth(self, *args, **kwargs) -> None:
        '''
        Obsolete function.
        '''
        print('>>> [Warning] The `smooth` method is obsolete.')
        print('    Please use `smooth_surface` and `smooth_sections` in the constructor.')
        print('    Or provide a user-defined guide curve `guide` in the `prepare()` method.')
        
    def bend(self, *args, **kwargs) -> None:
        '''
        Obsolete function.
        '''
        print('>>> [Warning] The `bend` method is obsolete.')
        print('    Please use `rotate_x_section` and `rotation_sections` in the constructor.')
        print('    Or provide a user-defined guide curve `guide` in the `prepare()` method.')


class OpenSurface(BasicSurface):
    '''
    Open surface defined by multiple OpenSection objects.
    '''
    def __init__(self, n_sec=1, name='Surf', nn=1001, ns=101, 
                    smooth_surface=False, smooth_sections : List[Tuple[int, int]] = None,
                    rotate_x_section=False, rotation_sections : List[Tuple[int, int]] = None,
                    is_guide_curve_at_LE=True):
        
        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, 
                            smooth_surface=smooth_surface, smooth_sections=smooth_sections,
                            rotate_x_section=rotate_x_section, rotation_sections=rotation_sections,
                            is_guide_curve_at_LE=is_guide_curve_at_LE)

        self.sections = [ OpenSection() for _ in range(max(1, n_sec)) ]

    def read_setting(self, fname) -> None:
        '''
        Read in Surface layout and CST parameters from file.

        Parameters
        ----------
        fname : str
            settings file name

        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3}

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
                        self.sections[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.sections[i].specified_thickness = float(line[5])

                        if self.is_2d:
                            self.sections[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_section):
                        i_line += 2
                        line = lines[i_line].split()
                        self.sections[i].cst = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    i_line += 2
                    line = lines[i_line].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        i_line += self.n_section*3
                        found_key = 0
                        continue

                    for i in range(self.n_section):

                        i_line += 2
                        line1 = lines[i_line].split()
                        cst_r = np.zeros(n_cst_refine)

                        i1 = 0

                        for j in range(n_cst_refine):
                            if j>=i_cst_start-1 and i1<len(line1):
                                cst_r[j] = float(line1[i1])
                                i1 += 1

                        self.sections[i].refine = cst_r

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                i_line += 1
        
        print('Read surface [%s] settings'%(self.name))

        self.update_layout_center()

    def update_section(self) -> None:
        '''
        Update all sections' profile curves.
        '''
        for i in range(self.n_section):
            self.sections[i].section(nn=self.nn, projection=False)


class Surface(BasicSurface):
    '''
    Surface defined by multiple Section objects, i.e., foils
    '''
    def __init__(self, n_sec=1, name='Surf', nn=1001, ns=101, 
                    smooth_surface=False, smooth_sections : List[Tuple[int, int]] = None,
                    rotate_x_section=False, rotation_sections : List[Tuple[int, int]] = None,
                    is_guide_curve_at_LE=True):
        
        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, 
                            smooth_surface=smooth_surface, smooth_sections=smooth_sections,
                            rotate_x_section=rotate_x_section, rotation_sections=rotation_sections,
                            is_guide_curve_at_LE=is_guide_curve_at_LE)
        
        self.sections = [ Section() for _ in range(max(1, n_sec)) ]

    def read_setting(self, fname, tail=0.0) -> None:
        '''
        Read in Surface layout and CST parameters from file

        Parameters
        ----------
        fname : str
            settings file name.
        tail : float or list
            tail thickness (m) of each section.
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3, 'CST_flip:': 4}

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
                        self.sections[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.sections[i].specified_thickness = float(line[5])

                        if isinstance(tail, float):
                            self.sections[i].tail  = tail/self.sections[i].chord
                        elif len(tail)==self.n_section:
                            self.sections[i].tail  = tail[i]/self.sections[i].chord
                        else:
                            raise Exception('tail must be a float or a list with length = section number')
                        
                        if self.sections[i].specified_thickness <= 0.0:
                            self.sections[i].specified_thickness = None

                        if self.is_2d:
                            self.sections[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_section):
                        i_line += 2
                        line = lines[i_line].split()
                        self.sections[i].cst_u = np.array([float(aa) for aa in line])

                        i_line += 1
                        line = lines[i_line].split()
                        self.sections[i].cst_l = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    i_line += 2
                    line = lines[i_line].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        i_line += self.n_section*3
                        found_key = 0
                        continue

                    for i in range(self.n_section):

                        i_line += 2
                        line1 = lines[i_line].split()

                        i_line += 1
                        line2 = lines[i_line].split()

                        cst_ur = np.zeros(n_cst_refine)
                        cst_lr = np.zeros(n_cst_refine)

                        i1 = 0
                        i2 = 0
                        for j in range(n_cst_refine):
                            if j>=i_cst_start-1 and i1<len(line1):
                                cst_ur[j] = float(line1[i1])
                                i1 += 1
                            if j>=i_cst_start-1 and i2<len(line2):
                                cst_lr[j] = float(line2[i2])
                                i2 += 1

                        self.sections[i].refine_u = cst_ur
                        self.sections[i].refine_l = cst_lr

                    found_key = 0

                elif found_surf and found_key == 4:
                    i_line += 2
                    line = lines[i_line].split()
                    n_cst_refine = int(line[0])

                    if n_cst_refine <= 0:
                        i_line += self.n_section*3
                        found_key = 0
                        continue

                    for i in range(self.n_section):

                        i_line += 2
                        line1 = lines[i_line].split()

                        i_line += 1
                        line2 = lines[i_line].split()

                        cst_ur = np.zeros(n_cst_refine)
                        cst_lr = np.zeros(n_cst_refine)

                        i1 = 0
                        i2 = 0
                        for j in range(n_cst_refine):
                            if i1<len(line1):
                                cst_ur[j] = float(line1[i1])
                                i1 += 1
                            if i2<len(line2):
                                cst_lr[j] = float(line2[i2])
                                i2 += 1

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                i_line += 1
        
        print('Read surface [%s] settings'%(self.name))

        # Locate layout center for plot
        self.update_layout_center()

    def update_section(self) -> None:
        '''
        Update all sections' profile curves.
        '''
        for i in range(self.n_section):
            self.sections[i].section(nn=self.nn, projection=False)

    def output_tecplot(self, fname=None, one_piece=False, split=False) -> None:
        '''
        Output the surface to `*.dat` in Tecplot format.

        Parameters
        ------------
        fname : str
            name of the output file.
        one_piece : bool
            if True, combine the span-wise sections into one piece.
        split : bool
            if True, split to upper and lower surfaces.
        '''
        if not split:
            super().output_tecplot(fname=fname, one_piece=one_piece)
            return

        if fname is None:
            fname = self.name + '.dat'

        n_sec   = 1 if self.is_2d else self.n_section-1
        n_piece = len(self.surfaces)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            if not one_piece:

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfaces[i_sec][0]
                    surf_y = self.surfaces[i_sec][1]
                    surf_z = self.surfaces[i_sec][2]

                    # surf_x[ns,nt], ns => spanwise
                    ns = self.ns
                    nt = int((surf_x.shape[1]+1)/2)

                    f.write('zone T="sec-u %d" i= %d j= %d \n'%(i_sec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                    f.write('zone T="sec-l %d" i= %d j= %d \n'%(i_sec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))

            else:
                
                n_point = n_sec*(self.ns-1) + 1
                ns = self.ns
                nt = int((self.surfaces[0][0].shape[1]+1)/2)
                
                f.write('zone T="sec-u" i= %d j= %d \n'%(nt, n_point))

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfaces[i_sec][0]
                    surf_y = self.surfaces[i_sec][1]
                    surf_z = self.surfaces[i_sec][2]

                    if i_sec>=n_piece-1:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                f.write('zone T="sec-l" i= %d j= %d \n'%(nt, n_point))

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfaces[i_sec][0]
                    surf_y = self.surfaces[i_sec][1]
                    surf_z = self.surfaces[i_sec][2]
                    
                    if i_sec>=n_piece-1:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))

    def output_plot3d(self, fname=None, split=False) -> None:
        '''
        Output the surface to `*.grd` in plot3d format.

        Parameters
        ------------
        fname : str
            name of the output file.
        split : bool
            if True, split to upper and lower surfaces.
        '''
        if not split:
            super().output_plot3d(fname=fname)
            return

        if fname is None:
            fname = self.name + '.grd'

        n_piece = len(self.surfaces)

        # surf_x[ns,nt], ns => spanwise
        ns = self.ns

        with open(fname, 'w') as f:

            f.write('%d \n '%(n_piece*2))   # Number of surfaces
            for i_sec in range(n_piece):
                nt = int((self.surfaces[i_sec][0].shape[1]+1)/2)
                f.write('%d %d 1\n '%(nt, ns))
                f.write('%d %d 1\n '%(nt, ns))

            for i_sec in range(n_piece):

                X = self.surfaces[i_sec][0]
                Y = self.surfaces[i_sec][1]
                Z = self.surfaces[i_sec][2]
                nt = int((X.shape[1]+1)/2)

                #* Upper surface
                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(X[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Y[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Z[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                #* Lower surface
                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(X[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Y[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Z[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

