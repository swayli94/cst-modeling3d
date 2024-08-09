'''
Classes and functions for the operation of surfaces.
'''

import numpy as np

from typing import List, Tuple

from .math import transform_curve
from .basic import BasicSection


class Lofting_2Profile():
    '''
    Create a surface by sweeping and blending two 2D profiles (unit curves),
    using a guide curve that runs through the leading edge point of each profile. 
    
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
    '''

    def __init__(self, sec0: BasicSection, sec1: BasicSection, n_spanwise=101, is_guide_curve_at_LE=True) -> None:
        
        self.sec0 = sec0
        self.sec1 = sec1
        
        self.n_spanwise = n_spanwise
        
        self.is_guide_curve_at_LE = is_guide_curve_at_LE
        
        self.check_sections()
        
        self.init_guide_curve()

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

    def sweep(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        '''
        Sweep the two 2D profiles along the guide curve.
        
        Returns
        -------
        surf_x, surf_y, surf_z : Tuple[np.ndarray, np.ndarray, np.ndarray]
            The 3D coordinates of the surface.
        '''
        #* Combine the upper and lower surface curves into one curve.
        if not self.sec0.is_open_curve:
        
            sec_x0 = np.concatenate((np.flip(self.sec0.xx), self.sec0.xx[1:]), axis=0)
            sec_y0 = np.concatenate((np.flip(self.sec0.yl), self.sec0.yu[1:]), axis=0)
            sec_x1 = np.concatenate((np.flip(self.sec1.xx), self.sec1.xx[1:]), axis=0)
            sec_y1 = np.concatenate((np.flip(self.sec1.yl), self.sec1.yu[1:]), axis=0)
            
        else:
            
            sec_x0 = self.sec0.xx
            sec_y0 = self.sec0.yy
            sec_x1 = self.sec1.xx
            sec_y1 = self.sec1.yy
        
        
        #* Initialize the 3D surface.
        n_point = len(sec_x0)
        surf_x = np.zeros((self.n_spanwise, n_point))
        surf_y = np.zeros((self.n_spanwise, n_point))
        surf_z = np.zeros((self.n_spanwise, n_point))
        
        
        #* Interpolate the two profiles and transform the 2D curve to a 3D curve.
        for i in range(self.n_spanwise):
            
            ratio = self.guide_curve['s'][i]
            
            xx = (1.0 - ratio) * sec_x0 + ratio * sec_x1
            yy = (1.0 - ratio) * sec_y0 + ratio * sec_y1
            
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



