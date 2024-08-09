'''
Classes and functions for the operation of surfaces.
'''

import numpy as np

from typing import Dict, Tuple, List

from scipy.interpolate import CubicSpline, interp1d
from scipy.spatial.transform import Rotation

from .math import transform_curve, angle_between_vectors
from .basic import BasicSection


class GuideCurve():
    '''
    Guide curve for surface lofting.
    
    Parameters
    ----------
    n_section : int
        The number of sections to loft.
        
    n_spanwise : int
        The number of spanwise points for the surface between the two sections.
        
    section_s_loc : List[float]
        The parametric coordinates of the sections along the guide curve, range in [0, 1].
        
    Attributes
    ----------
    n_total : int
        The total number of points for the guide curve.
        
    global_guide_curve : Dict[str, np.ndarray] [n_total]
        The global guide curve to sweep the sections.
    '''
    def __init__(self, n_section=2, n_spanwise=101, section_s_loc=[0.0, 1.0]) -> None:
        
        self.n_section = n_section
        self.n_spanwise = n_spanwise
        self.section_s_loc = section_s_loc
        
        self.n_total = (n_section-1) * (n_spanwise-1) + 1
        
        self.init_guide_curve()

    def init_guide_curve(self) -> None:
        '''
        Initialize the 3D guide curve, a straight line segment by default.
        '''
        self.global_guide_curve = {
            's':        np.zeros(self.n_total),
            'x':        np.zeros(self.n_total),
            'y':        np.zeros(self.n_total),
            'z':        np.zeros(self.n_total),
            'scale':    np.ones(self.n_total),
            'rot_x':    np.zeros(self.n_total),
            'rot_y':    np.zeros(self.n_total),
            'rot_z':    np.zeros(self.n_total),
        }
        
        for i_sec in range(self.n_section-1):
            
            s0 = self.section_s_loc[i_sec]
            s1 = self.section_s_loc[i_sec+1]
            
            for i in range(self.n_spanwise):
                
                ii = i_sec * (self.n_spanwise-1) + i
                
                self.global_guide_curve['s'][ii] = s0 + (s1 - s0) * i / (self.n_spanwise-1.0)
        
        self.global_guide_curve['z'] = self.global_guide_curve['s']

    def get_local_parametric_coordinate(self, i_sec0: int) -> Tuple[np.ndarray, int, int]:
        '''
        Get the local parametric coordinates of the guide curve.
        
        Parameters
        ----------
        i_sec0 : int
            The index of the section with smaller index of the two sections that construct a surface.
        
        Returns
        -------
        local_s : np.ndarray [n_spanwise]
            The local parametric coordinates of the guide curve, range in [0, 1].
            
        index0, index1 : int
            The indices of the two sections in the global guide curve.
        '''
        s0 = self.section_s_loc[i_sec0]
        s1 = self.section_s_loc[i_sec0+1]
        
        local_s = np.zeros(self.n_spanwise)

        for i in range(self.n_spanwise):
            
            ii = i_sec0 * (self.n_spanwise-1) + i
            
            s = self.global_guide_curve['s'][ii]
            
            local_s[i] = (s - s0) / (s1 - s0)

        index0 = i_sec0 * (self.n_spanwise-1)
        index1 = (i_sec0+1) * (self.n_spanwise-1)
        
        return local_s, index0, index1

    def get_local_guide_curve(self, i_sec0: int) -> Dict[str, np.ndarray]:
        '''
        Get the local guide curve.
        
        Parameters
        ----------
        i_sec0 : int
            The index of the section with smaller index of the two sections that construct a surface.
        '''
        local_s, index0, index1 = self.get_local_parametric_coordinate(i_sec0)
        
        guide_curve = {
            's':        local_s,
            'x':        self.global_guide_curve['x'][index0:index1+1],
            'y':        self.global_guide_curve['y'][index0:index1+1],
            'z':        self.global_guide_curve['z'][index0:index1+1],
            'scale':    self.global_guide_curve['scale'][index0:index1+1],
            'rot_x':    self.global_guide_curve['rot_x'][index0:index1+1],
            'rot_y':    self.global_guide_curve['rot_y'][index0:index1+1],
            'rot_z':    self.global_guide_curve['rot_z'][index0:index1+1],
        }
                
        return guide_curve
        
    def update_with_value(self, **kwargs) -> None:
        '''
        Update the global guide curve.
        
        Parameters
        ----------
        kwargs : Dict[str, np.ndarray]
            The keyword arguments to update the guide curve.
            Including 'x', 'y', 'z', 'scale', 'rot_x', 'rot_y', 'rot_z'.
        '''
        for key, value in kwargs.items():
            if key in self.global_guide_curve:
                self.global_guide_curve[key] = value

    def generate_by_spline(self, global_control_s: np.ndarray, global_values: np.ndarray, 
                                slope_s0=None, slope_s1=None, key='x', periodic=False) -> None:
        '''
        Update the global guide curve by a spline interpolation.
        
        Parameters
        ----------
        global_control_s : np.ndarray [n_point]
            The global parametric coordinates of the control points, range in [0, 1].
            
        global_values : np.ndarray [n_point]
            The global values of the control points.
            
        slope_s0, slope_s1 : float
            The slope of the curve at the start and end points.
            
        key : str
            The key of the global guide curve to update.
        '''        
        if periodic:
            
            bc_type = 'periodic'
            
        else:
            
            bcx0 = (2, 0.0)
            bcx1 = (2, 0.0)
            
            if slope_s0 is not None:
                bcx0 = (1, slope_s0)
                
            if slope_s1 is not None:
                bcx1 = (1, slope_s1)
            
            bc_type = (bcx0, bcx1)
        
        func = CubicSpline(global_control_s, global_values, bc_type=bc_type)
        
        self.global_guide_curve[key] = func(self.global_guide_curve['s'])
        
    def generate_by_interp1d(self, global_control_s: np.ndarray, global_values: np.ndarray, key='x', kind='linear') -> None:
        '''
        Update the global guide curve by a linear interpolation.
        
        Parameters
        ----------
        global_control_s : np.ndarray [n_point]
            The global parametric coordinates of the control points, range in [0, 1].
            
        global_values : np.ndarray [n_point]
            The global values of the control points.
            
        key : str
            The key of the global guide curve to update.
            
        kind : str
            The kind of interpolation. See scipy.interpolate.interp1d.
        '''
        func = interp1d(global_control_s, global_values, kind=kind, fill_value='extrapolate')
        
        self.global_guide_curve[key] = func(self.global_guide_curve['s'])
        
    def update_rotation_angle_with_tangent(self, key='all') -> None:
        '''
        Update the rotation angle with the tangent of the guide curve.
        
        Parameters
        ----------
        key : str
            The key of the global guide curve to update.
            It can be 'all', 'rot_x', 'rot_y', 'rot_z'.
        '''
        for i_global in range(self.n_total-1):
            
            x0 = self.global_guide_curve['x'][i_global]
            y0 = self.global_guide_curve['y'][i_global]
            z0 = self.global_guide_curve['z'][i_global]
            
            x1 = self.global_guide_curve['x'][i_global+1]
            y1 = self.global_guide_curve['y'][i_global+1]
            z1 = self.global_guide_curve['z'][i_global+1]
            
            if key == 'all':
                
                vec0 = np.array([0, 0, 1.0])
                vec1 = np.array([x1-x0, y1-y0, z1-z0])
                
                rot, _ = Rotation.align_vectors(vec0, vec1)
                
                angles = rot.as_euler('zxy', degrees=True)
                
                self.global_guide_curve['rot_z'][i_global] = angles[0]
                self.global_guide_curve['rot_x'][i_global] = angles[1]
                self.global_guide_curve['rot_y'][i_global] = angles[2]
                
            elif key == 'rot_z':
                
                vec0 = np.array([1.0, 0.0, 0.0])
                vec1 = np.array([x1-x0, y1-y0, 0.0])
                angle = angle_between_vectors(vec0, vec1)
                sign = np.sign(vec1[1])

                self.global_guide_curve['rot_z'][i_global] = angle * sign
                
            elif key == 'rot_x':
                
                vec0 = np.array([0.0, 0.0, 1.0])
                vec1 = np.array([0.0, y1-y0, z1-z0])
                angle = angle_between_vectors(vec0, vec1)
                sign = - np.sign(vec1[1])
                
                self.global_guide_curve['rot_x'][i_global] = angle * sign
                
            elif key == 'rot_y':
                
                vec0 = np.array([0.0, 0.0, 1.0])
                vec1 = np.array([x1-x0, 0.0, z1-z0])
                angle = angle_between_vectors(vec0, vec1)
                sign = np.sign(vec1[0])
                
                self.global_guide_curve['rot_y'][i_global] = angle * sign
                
            else:
                
                raise ValueError('Invalid key.', key)                
                
        #* The last point
        if key == 'all' or key == 'rot_z':
            self.global_guide_curve['rot_z'][-1] = self.global_guide_curve['rot_z'][-2]
            
        if key == 'all' or key == 'rot_x':
            self.global_guide_curve['rot_x'][-1] = self.global_guide_curve['rot_x'][-2]
            
        if key == 'all' or key == 'rot_y':
            self.global_guide_curve['rot_y'][-1] = self.global_guide_curve['rot_y'][-2]
        
    def output(self, fname='guide-curve.dat'):
        '''
        Output the guide curve to a file.
        '''
        with open(fname, 'w') as f:
            
            f.write('Variables= "s", "x", "y", "z", "scale", "rot_x", "rot_y", "rot_z" \n')
            f.write('zone T="guide-curve"  i= %d \n' % self.n_total)
            
            for i in range(self.n_total):
                
                f.write('%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n' % (
                    self.global_guide_curve['s'][i], 
                    self.global_guide_curve['x'][i], 
                    self.global_guide_curve['y'][i], 
                    self.global_guide_curve['z'][i], 
                    self.global_guide_curve['scale'][i], 
                    self.global_guide_curve['rot_x'][i], 
                    self.global_guide_curve['rot_y'][i], 
                    self.global_guide_curve['rot_z'][i]))

        
class Lofting_2Profile():
    '''
    Create a surface by sweeping and blending two 2D profiles (unit curves),
    using a guide curve that runs through the leading/trailing edge point of each profile. 
    
    Parameters
    ----------
    sec0, sec1 : BasicSection
        The section objects to provide the 2D profiles.
        
    n_spanwise : int
        The number of spanwise points for the surface between the two sections.
        
    is_guide_curve_at_LE : bool
        If True, the guide curve runs through the leading edge point (0,0) of each 2D profile.
        Otherwise, it runs through the point (1,0) of each 2D profile, which is the trailing edge of a unit curve.
        
        This changes the scaling and rotation center of the sections, 
        which consequently affects the actual x,y-coordinates of the leading edge.
        When the rotation center is at the leading edge, 
        the coordinates of the leading edge are directly defined by the xLE, yLE, zLE attributes of the sections.
        However, when the rotation center is at the trailing edge, 
        the coordinates of the leading edge are defined by the xLE, yLE, zLE attributes of the sections, and the scaling factor and rotation angle.
        
    Attributes
    ----------
    parametric_coord : np.ndarray [n_spanwise]
        The parametric coordinates of the guide curve, range in [0, 1].
    
    guide_curve : Dict[str, np.ndarray] [n_spanwise]
        The guide curve to sweep the sections.
        
        - 's' : The parametric coordinates of the guide curve, range in [0, 1].
        - 'x' : The x-coordinates of the guide curve.
        - 'y' : The y-coordinates of the guide curve.
        - 'z' : The z-coordinates of the guide curve.
        - 'scale' : The scaling factor of the section at `s`.
        - 'rot_x' : The rotation angle about the x-axis of the section at `s`.
        - 'rot_y' : The rotation angle about the y-axis of the section at `s`.
        - 'rot_z' : The rotation angle about the z-axis of the section at `s`.
        
    profile_0, profile_1 : List[np.ndarray] [2][n_point]
        The 2D profiles for interpolation, [profile_x, profile_y].
        
    n_point : int
        The number of points in the 2D profiles.
    '''

    def __init__(self, sec0: BasicSection, sec1: BasicSection, n_spanwise=101, is_guide_curve_at_LE=True) -> None:
        
        self.sec0 = sec0
        self.sec1 = sec1
        
        self.n_spanwise = n_spanwise
        
        self.is_guide_curve_at_LE = is_guide_curve_at_LE
        
        self.check_sections()
        
        self.init_guide_curve()
        
        self.init_profile()

    def check_sections(self) -> None:
        '''
        Check the two sections.
        '''
        if not self.sec0.has_profile or not self.sec1.has_profile:
            raise ValueError('The two 2D profiles must be provided.')
        
        if self.sec0.n_point_profile != self.sec1.n_point_profile:
            raise ValueError('The two 2D profiles must have the same number of points.')
        
    def init_guide_curve(self) -> None:
        '''
        Initialize the 3D guide curve, a straight line segment by default.
        '''
        self.guide_curve = {
            's':        np.linspace(0.0,                1.0,                self.n_spanwise, endpoint=True),
            'x':        np.linspace(self.sec0.xLE,      self.sec1.xLE,      self.n_spanwise, endpoint=True),
            'y':        np.linspace(self.sec0.yLE,      self.sec1.yLE,      self.n_spanwise, endpoint=True),
            'z':        np.linspace(self.sec0.zLE,      self.sec1.zLE,      self.n_spanwise, endpoint=True),
            'scale':    np.linspace(self.sec0.scale,    self.sec1.scale,    self.n_spanwise, endpoint=True),
            'rot_x':    np.linspace(self.sec0.rot_x,    self.sec1.rot_x,    self.n_spanwise, endpoint=True),
            'rot_y':    np.linspace(self.sec0.rot_y,    self.sec1.rot_y,    self.n_spanwise, endpoint=True),
            'rot_z':    np.linspace(self.sec0.rot_z,    self.sec1.rot_z,    self.n_spanwise, endpoint=True),
        }

        if not self.is_guide_curve_at_LE:
            
            self.guide_curve['x'] = self.guide_curve['x'] + self.guide_curve['scale']

    def init_profile(self) -> None:
        '''
        Initialize the 2D profiles for interpolation.
        '''
        self.profile_0 : List[np.ndarray] = [None, None]
        self.profile_1 : List[np.ndarray] = [None, None]
        
        #* Combine the upper and lower surface curves into one curve.
        if not self.sec0.is_open_curve:
        
            self.profile_0[0] = np.concatenate((np.flip(self.sec0.xx), self.sec0.xx[1:]), axis=0)
            self.profile_0[1] = np.concatenate((np.flip(self.sec0.yl), self.sec0.yu[1:]), axis=0)
            self.profile_1[0] = np.concatenate((np.flip(self.sec1.xx), self.sec1.xx[1:]), axis=0)
            self.profile_1[1] = np.concatenate((np.flip(self.sec1.yl), self.sec1.yu[1:]), axis=0)
            
        else:
            
            self.profile_0[0] = self.sec0.xx
            self.profile_0[1] = self.sec0.yy
            self.profile_1[0] = self.sec1.xx
            self.profile_1[1] = self.sec1.yy
            
        self.n_point = self.profile_0[0].shape[0]
        
    def update_guide_curve(self, **kwargs) -> None:
        '''
        Update the guide curve.
        
        Parameters
        ----------
        kwargs : Dict[str, np.ndarray]
            The keyword arguments to update the guide curve.
            Including 'x', 'y', 'z', 'scale', 'rot_x', 'rot_y', 'rot_z'.
        '''
        for key, value in kwargs.items():
            if key in self.guide_curve:
                self.guide_curve[key] = value

    def sweep(self, spanwise_profiles=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Sweep the two 2D profiles along the guide curve.
        
        Parameters
        ----------
        spanwise_profiles : List[List[np.ndarray]] [n_spanwise][2][n_point]
            The profiles for transformation along the spanwise direction.
            If None, the profiles are interpolated linearly from the two end sections.
        
        Returns
        -------
        surf_x, surf_y, surf_z : Tuple[np.ndarray, np.ndarray, np.ndarray]
            The 3D coordinates of the surface.
        '''
        #* Initialize the 3D surface.
        surf_x = np.zeros((self.n_spanwise, self.n_point))
        surf_y = np.zeros((self.n_spanwise, self.n_point))
        surf_z = np.zeros((self.n_spanwise, self.n_point))
        
        #* Interpolate the two profiles and transform the 2D curve to a 3D curve.
        for i in range(self.n_spanwise):
            
            ratio = self.guide_curve['s'][i]
            
            if spanwise_profiles is None:
            
                xx = (1.0 - ratio) * self.profile_0[0] + ratio * self.profile_1[0]
                yy = (1.0 - ratio) * self.profile_0[1] + ratio * self.profile_1[1]
                
            else:
                
                xx = spanwise_profiles[i][0]
                yy = spanwise_profiles[i][1]
            
            gx = self.guide_curve['x'][i]
            gy = self.guide_curve['y'][i]
            gz = self.guide_curve['z'][i]           

            if self.is_guide_curve_at_LE:
                
                dx = gx
                x0 = gx
                xr = gx

            else:
                
                dx = gx - self.guide_curve['scale'][i]
                x0 = dx
                xr = gx

            surf_x[i], surf_y[i], surf_z[i] = transform_curve(xx, yy, 
                    dx=dx, dy=gy, dz=gz,
                    scale=self.guide_curve['scale'][i], x0=x0, y0=gy,
                    rot_x=self.guide_curve['rot_x'][i], rot_y=self.guide_curve['rot_y'][i], rot_z=self.guide_curve['rot_z'][i],
                    xr=xr, yr=gy, zr=gz)
                    
        return surf_x, surf_y, surf_z


class Lofting():
    '''
    Create a surface by sweeping and blending multiple 2D profiles (unit curves),
    using a guide curve that runs through the leading/trailing edge point of each profile. 
    
    Parameters
    ----------
    sections : List[BasicSection]
        The section objects to provide the 2D profiles.
        
    global_guide_curve : GuideCurve
        The global guide curve to sweep the sections.
        
    is_guide_curve_at_LE : bool
        If True, the guide curve runs through the leading edge point (0,0) of each 2D profile.
        Otherwise, it runs through the point (1,0) of each 2D profile, which is the trailing edge of a unit curve.
        
        This changes the scaling and rotation center of the sections, 
        which consequently affects the actual x,y-coordinates of the leading edge.
        When the rotation center is at the leading edge, 
        the coordinates of the leading edge are directly defined by the xLE, yLE, zLE attributes of the sections.
        However, when the rotation center is at the trailing edge, 
        the coordinates of the leading edge are defined by the xLE, yLE, zLE attributes of the sections, and the scaling factor and rotation angle.
        
    Attributes
    ----------
    n_section : int
        The number of sections to loft.
    
    n_spanwise : int
        The number of spanwise points for the surface between the two sections.
        
    n_total : int
        The total number of points for the guide curve.
    
    surfs : List[List[np.ndarray]] [n_section-1][3][n_spanwise, n_point]
        The 3D coordinates of the surfaces.
    '''
    def __init__(self, sections: List[BasicSection], global_guide_curve: GuideCurve, 
                    is_guide_curve_at_LE=True) -> None:
        
        self.sections = sections
        self.guide_curve = global_guide_curve
        
        self.n_section = global_guide_curve.n_section
        self.n_spanwise = global_guide_curve.n_spanwise
        self.n_total = global_guide_curve.n_total
        
        self.is_guide_curve_at_LE = is_guide_curve_at_LE
        
        self.surfs : List[List[np.ndarray]] = []
        
        if len(sections) != self.n_section:
            raise ValueError('The number of sections must be consistent with the guide curve.')
    
    def init_all_profiles(self) -> None:
        '''
        Initialize the 2D profiles for interpolation.
        
        This is only needed when calling `create_spanwise_profiles`.
        '''
        self.profiles : List[List[np.ndarray]] = [[None, None] for _ in range(self.n_section)]
        
        self.n_point = None
        
        for i_sec in range(self.n_section):
            
            sec = self.sections[i_sec]
        
            #* Combine the upper and lower surface curves into one curve.
            if not self.sections[i_sec].is_open_curve:
            
                self.profiles[i_sec][0] = np.concatenate((np.flip(sec.xx), sec.xx[1:]), axis=0)
                self.profiles[i_sec][1] = np.concatenate((np.flip(sec.yl), sec.yu[1:]), axis=0)

            else:
                
                self.profiles[i_sec][0] = sec.xx
                self.profiles[i_sec][1] = sec.yy

            if self.n_point is None:
                self.n_point = self.profiles[i_sec][0].shape[0]
                
            elif self.n_point != self.profiles[i_sec][0].shape[0]:
                raise ValueError('The number of points in the 2D profiles must be consistent.')
    
    def create_spanwise_profiles(self, kind='linear') -> List[List[np.ndarray]]:
        '''
        Create the spanwise profiles for the surface.
        
        Returns
        -------
        spanwise_profiles : List[List[np.ndarray, np.ndarray]] [n_total][2][n_point]
            The profiles for transformation along the spanwise direction.
        '''
        ss = self.guide_curve.global_guide_curve['s']
        
        self.init_all_profiles()
        
        spanwise_profiles : List[List[np.ndarray]] = [[np.zeros(self.n_point), np.zeros(self.n_point)] for _ in range(self.n_total)]
    
        for i_point in range(self.n_point):
            
            control_point_s = self.guide_curve.section_s_loc
            
            control_point_x = [self.profiles[i_sec][0][i_point] for i_sec in range(self.n_section)]
            control_point_y = [self.profiles[i_sec][1][i_point] for i_sec in range(self.n_section)]
            
            func_x = interp1d(control_point_s, control_point_x, kind=kind, fill_value='extrapolate')
            func_y = interp1d(control_point_s, control_point_y, kind=kind, fill_value='extrapolate')
            
            xx = func_x(ss)
            yy = func_y(ss)
            
            for i_span in range(self.n_total):
                
                spanwise_profiles[i_span][0][i_point] = xx[i_span]
                spanwise_profiles[i_span][1][i_point] = yy[i_span]
                
        return spanwise_profiles
        
    def sweep(self, interp_profile_kind=None) -> List[List[np.ndarray]]:
        '''
        Sweep the profiles along the guide curve.
        
        Parameters
        ----------
        interp_profile_kind : str
            The kind of interpolation for the spanwise profiles.
            If None, the profiles are interpolated linearly from the two end sections.
            Otherwise, the profiles are interpolated by the specified kind. See scipy.interpolate.interp1d.
        
        Returns
        -------
        surfs : List[List[np.ndarray]] [n_section-1][3][n_spanwise, n_point]
            The 3D coordinates of the surfaces.
        '''
        if interp_profile_kind is not None:
            spanwise_profiles = self.create_spanwise_profiles(kind=interp_profile_kind)
        
        for i_sec in range(self.n_section-1):
            
            sec0 = self.sections[i_sec]
            sec1 = self.sections[i_sec+1]
            
            loft = Lofting_2Profile(sec0, sec1, n_spanwise=self.n_spanwise, is_guide_curve_at_LE=self.is_guide_curve_at_LE)
            
            guide_curve = self.guide_curve.get_local_guide_curve(i_sec)
            
            loft.update_guide_curve(**guide_curve)
            
            if interp_profile_kind is not None:
                
                index0 = i_sec * (self.n_spanwise-1)
                index1 = (i_sec+1) * (self.n_spanwise-1)
                
                spanwise_profiles_2 = spanwise_profiles[index0:index1+1]
                
            else:
                
                spanwise_profiles_2 = None
                    
            surf_x, surf_y, surf_z = loft.sweep(spanwise_profiles_2)
            
            self.surfs.append([surf_x, surf_y, surf_z])
            
        return self.surfs



