'''
Basic classes for sections and surfaces, and fundamental functions
'''
import copy
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import CubicSpline

from .math import rotate, transform_foil, stretch_fixed_point, to_cylinder, scaled_sigmoid
from .section import BasicSection, interpolate_basic_sections  # re-exported for backward compatibility
from .section import Section, OpenSection
from .io import read_surface_setting_blocks, apply_surface_layout




class BasicSurface():
    '''
    Multi-section surface based on `BasicSection`.
    
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
        
    projection : bool
        whether keeps the projection length the same when rotating the section, by default True.

    Examples
    --------
    >>> surf = BasicSurface(n_section=1, name='Surf', nn=1001, ns=101, projection=True)
    
    
    Notes
    -------
    - +x:     flow direction (m)
    - +y:     upside (m)
    - +z:     span-wise (m)

        
    Attributes
    ------------
    is_2d : bool
        whether this is a 3D surface for a 2D curve (unit span). 
        It is `True` when `n_section` is 0 or 1.
        
    sections : list of BasicSection
        section objects.
        
    surfaces : list of list of ndarray
        surface coordinates. List `surfaces` contains `n_section`-1 sub-lists.
        Each sub-list is the coordinates of the 3D surface, which contains 3 `ndarray`.
        The 3 arrays are the X, Y, Z coordinates of the surface.   
        For example, surfaces = [[x0, y0, z0], [x1, y1, z1]] when n_section=3.
        
    half_size : float
        half span for plotting.
        
    center : ndarray
        surface center for plotting.

    '''

    def __init__(self, n_section=1, name='Surf', nn=1001, ns=101, projection=True):
        '''
        Construct a BasicSurface object.
        
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
        projection : bool
            whether keeps the projection length the same when rotating the section, by default True.
            
        Examples
        --------
        >>> surf = BasicSurface(n_section=1, name='Surf', nn=1001, ns=101, projection=True)
        '''
        
        n_ = max(1, n_section)
        self.is_2d   = n_ == 1    # type: bool
        self.name  = name       # type: str
        self.nn    = nn         # type: int
        self.ns    = ns         # type: int
        self.sections  = [ BasicSection() for _ in range(n_) ]
        self.surfaces = []         # type: list[list]
        self.projection = projection  # type: bool

        # Parameters for plot
        self.half_size = 0.5    # type: float
        self.center = np.array([0.5, 0.5, 0.5])

    @property
    def n_section(self) -> int:
        '''
        Number of sections
        '''
        return len(self.sections)

    @property
    def spanwise_locations(self) -> List[float]:
        '''
        List of section zLE
        '''
        return [round(sec.zLE,5) for sec in self.sections]

    def read_setting(self, fname: str) -> None:
        '''
        Read in Surface layout parameters from file
        
        Parameters
        -----------
        fname : str
            parameter file name.

        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface setting')
        
        blocks = read_surface_setting_blocks(fname, self.name, self.n_section, {'Layout:': 1}, 1)
        
        if 'Layout:' in blocks:
            apply_surface_layout(self.sections, blocks['Layout:'], self.is_2d)
        
        self.update_layout_center()

    def update_layout_center(self) -> None:
        '''
        Locate layout center for plot
        '''
        x_range = [self.sections[0].xLE, self.sections[0].xLE+self.sections[0].chord]
        y_range = [self.sections[0].yLE, self.sections[0].yLE]
        z_range = [self.sections[0].zLE, self.sections[0].zLE+self.is_2d]
        
        if not self.is_2d:
            
            for i in range(self.n_section):
                x_range[0] = min(x_range[0], self.sections[i].xLE)
                x_range[1] = max(x_range[1], self.sections[i].xLE+self.sections[i].chord)
                y_range[0] = min(y_range[0], self.sections[i].yLE)
                y_range[1] = max(y_range[1], self.sections[i].yLE)
                z_range[0] = min(z_range[0], self.sections[i].zLE)
                z_range[1] = max(z_range[1], self.sections[i].zLE)
        
        span = np.array([x_range[1]-x_range[0], y_range[1]-y_range[0], z_range[1]-z_range[0]])
        self.half_size = span.max()/2.0
        self.center[0] = 0.5*(x_range[1]+x_range[0])
        self.center[1] = 0.5*(y_range[1]+y_range[0])
        self.center[2] = 0.5*(z_range[1]+z_range[0])

    def linear_interpolate_z(self, z: float, key='x') -> dict:
        '''
        Linear interpolation of key by a given z.

        Parameters
        ----------
        z : float
            location of the value
        key : str
            the value to be interpolated:
            'x' or 'X', 'y' or 'Y',
            'c' or 'C' or 'chord',
            't' or 'thick' or 'thickness', 'twist'.
        '''
        #* Find the adjacent control sections
        i_sec = self.n_section
        for i in range(self.n_section-1):
            if (z-self.sections[i].zLE)*(z-self.sections[i+1].zLE)<0 or z==self.sections[i].zLE or z==self.sections[i+1].zLE:
                i_sec = i

        if i_sec >= self.n_section:
            raise Exception('z is not within the surface: ', z, self.sections[0].zLE, self.sections[-1].zLE)

        #* Linear interpolation
        tt = (z-self.sections[i_sec].zLE)/(self.sections[i_sec+1].zLE-self.sections[i_sec].zLE)
        key_value = None

        if key == 'x' or key == 'X':
            key_value = (1-tt)*self.sections[i_sec].xLE + tt*self.sections[i_sec+1].xLE
        elif key == 'y' or key == 'Y':
            key_value = (1-tt)*self.sections[i_sec].yLE + tt*self.sections[i_sec+1].yLE
        elif key == 'c' or key == 'C' or  key == 'chord':
            key_value = (1-tt)*self.sections[i_sec].chord + tt*self.sections[i_sec+1].chord
        elif key == 't' or key == 'thick' or key == 'thickness':
            key_value = (1-tt)*self.sections[i_sec].thick + tt*self.sections[i_sec+1].thick
        elif key == 'twist':
            key_value = (1-tt)*self.sections[i_sec].twist + tt*self.sections[i_sec+1].twist
        else:
            raise Exception('Unknown key:', key)

        return key_value


    def update_sections(self, flip_x=False) -> None:
        '''
        Update surface sections, including the construction of 2D unit curves (optional)
        and transforming to 3D curves.

        Parameters
        ----------
        flip_x : bool
            whether flip `xx` in the reverse order, by default False.

        '''
        for i in range(self.n_section):
            self.sections[i].geo(nn=self.nn, flip_x=flip_x, projection=self.projection)

    def geo(self, flip_x=False, update_sec=True) -> None:
        '''
        Generate surface geometry.
        
        First, update sections by calling `sec.geo()` (optional): \n
            1) update the 2D curve `sec.xx, sec.yy, sec.yu, sec.yl`; \n
            2) transform the 2D curve to the 3D curve `sec.x, sec.y, sec.z`; \n
        Then, interpolate the 3D surface `[surf_x, surf_y, sur_z]` from 3D curves.
        
        
        Parameters
        -----------
        flip_x : bool
            whether flip `xx` in the reverse order, by default False.
        update_sec : bool
            whether update sections, by default True.
        
        '''
        #* Update sections
        if update_sec:
            self.update_sections(flip_x=flip_x)

        #* Interpolate the 3D surface from 3D curves.
        self.surfaces = []
        
        if self.is_2d:
            
            sec_ = copy.deepcopy(self.sections[0])
            sec_.zLE = 1.0
            sec_.z = np.ones_like(sec_.x)
            surf = self.sections_to_surface(self.sections[0], sec_, ns=self.ns)
            self.surfaces.append(surf)

        else:
            
            for i in range(self.n_section-1):
                surf = self.sections_to_surface(self.sections[i], self.sections[i+1], ns=self.ns)
                self.surfaces.append(surf)

    def geo_axisymmetric(self, phi, flip_x=False, update_sec=True) -> None:
        '''
        Generate axisymmetric surface geometry.
        
        Parameters
        -----------
        phi : list or ndarray
            position angle of control sections.
        flip_x : bool
            whether flip `xx` in the reverse order, by default False.
        update_sec : bool
            whether update sections, by default True.

        '''
        #* Update sections
        if update_sec:
            self.update_sections(flip_x=flip_x)

        #* Interpolate the 3D surface from 3D curves.
        self.surfaces = []

        if self.is_2d:
            raise Exception('Axisymmetric geometry can not be 2D surface')

        else:
            for i in range(self.n_section-1):
                surf = self.sections_to_axisymmetric_surface(self.sections[i], self.sections[i+1], phi[i], phi[i+1], ns=self.ns)
                self.surfaces.append(surf)


    @staticmethod
    def sections_to_surface(sec0: BasicSection, sec1: BasicSection, ns=101) -> List[np.ndarray]:
        '''
        Interpolate surface from section 3D curves.
        
        Parameters
        ----------
        sec0, sec1 : BasicSection
            sections on both ends of the surface.
        ns : int
            number of points in the interpolation direction.
            
        Returns
        ---------
        surf : list of ndarray
            coordinates of the surface, `[surf_x, surf_y, surf_z]`,
            `surf_x`'s shape is `[ns, nn]`. 
        
        Examples
        ---------
        >>> surf = sections_to_surface(sec0, sec1, ns)

        '''
        if sec0.x.shape[0]<=1 or sec1.x.shape[0]<=1:
            raise Exception('The 3D curve (sec.x, sec.y, sec.z) is not available')
        
        nn = sec0.x.shape[0]
        surf_x = np.zeros((ns,nn))
        surf_y = np.zeros((ns,nn))
        surf_z = np.zeros((ns,nn))
        
        for i in range(ns):
            tt = 1.0*i/(ns-1.0)
            surf_x[i,:] = (1-tt)*sec0.x + tt*sec1.x
            surf_y[i,:] = (1-tt)*sec0.y + tt*sec1.y
            surf_z[i,:] = (1-tt)*sec0.z + tt*sec1.z

        surf = [surf_x, surf_y, surf_z]

        return surf

    @staticmethod
    def sections_to_axisymmetric_surface(sec0: BasicSection, sec1: BasicSection, 
                                  phi0: float, phi1: float, ns=101) -> List[np.ndarray]:
        '''
        Interpolate axisymmetric surface section between curves.
        
        Parameters
        ----------
        sec0, sec1 : BasicSection
            sections on both ends of the surface.
        phi0, phi1 : float
            angle (degree) about X-axis (X-Y plane is 0 degree).
        ns : int
            number of points in the interpolation direction.
            
        Returns
        ---------
        surf : list of ndarray
            coordinates of the surface, `[surf_x, surf_y, surf_z]`,
            `surf_x`'s shape is `[ns, nn]`. 
            
        Examples
        ---------
        >>> surf = sections_to_axisymmetric_surface(sec0, sec1, phi0, phi1, ns)
        '''
        if sec0.x.shape[0]<=1 or sec1.x.shape[0]<=1:
            raise Exception('The 3D curve (sec.x, sec.y, sec.z) is not available')
        
        nn = sec0.x.shape[0]
        surf_x = np.zeros((ns,nn))
        surf_y = np.zeros((ns,nn))
        surf_z = np.zeros((ns,nn))
        xx = np.zeros(nn)
        yy = np.zeros(nn)
        zz = np.zeros(nn)

        R = np.sqrt(sec0.yLE**2+sec0.zLE**2)
        
        for i in range(ns):

            tt    = 1.0*i/(ns-1.0)
            t0    = 1-tt

            xLE   = t0*sec0.xLE + tt*sec1.xLE
            yLE_  = t0*sec0.yLE + tt*sec1.yLE
            zLE_  = t0*sec0.zLE + tt*sec1.zLE

            angle = t0*phi0 + tt*phi1
            yLE   = R*np.cos(angle/180.0*np.pi)
            zLE   = R*np.sin(angle/180.0*np.pi)

            xx = t0*sec0.x + tt*sec1.x
            yy = t0*sec0.y + tt*sec1.y + (yLE-yLE_)
            zz = t0*sec0.z + tt*sec1.z + (zLE-zLE_)

            surf_x[i,:], surf_y[i,:], surf_z[i,:] = rotate(xx, yy, zz, angle=angle, origin=[xLE, yLE, zLE], axis='X')

        return [surf_x, surf_y, surf_z]


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

    def translate(self, dx=0.0, dy=0.0, dz=0.0) -> None:
        '''
        Translate surface coordinates
        '''
        for surf in self.surfaces:
            surf[0] += dx
            surf[1] += dy
            surf[2] += dz

        self.center[0] += dx
        self.center[1] += dy
        self.center[2] += dz

    def scale(self, scale=1.0, X0=0.0, Y0=0.0, Z0=0.0) -> None:
        '''
        Scale surface coordinates about point `[X0, Y0, Z0]`.
        '''
        for surf in self.surfaces:
            surf[0] = (surf[0]-X0)*scale + X0
            surf[1] = (surf[1]-Y0)*scale + Y0
            surf[2] = (surf[2]-Z0)*scale + Z0

        self.center[0] = (self.center[0]-X0)*scale + X0
        self.center[1] = (self.center[1]-Y0)*scale + Y0
        self.center[2] = (self.center[2]-Z0)*scale + Z0


    def split(self, ips: list) -> None:
        '''
        Split each surface `surfaces` into several pieces.
        Length of `surfaces` is multiplied by len(ips)+1.
        
        Parameters
        -----------
        ips : list of int
            split point indices on each section curves

        Notes
        ----------
        `surf` is a list of ndarray, i.e., `[surf_x, surf_y, surf_z]`,
        the shape of each element is `[ns, nn]`.
        '''
        nt = self.surfaces[0][0].shape[1]  # i.e., self.nn
        ips.sort()
        ips = [1] + ips + [nt]
        
        surfaces = copy.deepcopy(self.surfaces)
        self.surfaces = []

        for i in range(len(ips)-1):

            for surf in surfaces:

                surf_x = surf[0]
                surf_y = surf[1]
                surf_z = surf[2]
                
                self.surfaces.append([
                    surf_x[:,ips[i]-1:ips[i+1]],
                    surf_y[:,ips[i]-1:ips[i+1]],
                    surf_z[:,ips[i]-1:ips[i+1]]])

    def smooth(self, i_sec0: int, i_sec1: int, smooth0=False, smooth1=False, 
               dyn0=None, ratio_end=10) -> None:
        '''
        Smooth the span-wise curve between i_sec0 and i_sec1.
        
        Parameters
        -----------
        i_sec0, i_sec1: int
            The index of surface `surfaces` to be smoothed is `i_sec0, ..., i_sec1-1`.
            They usually are the starting and ending section index of the smooth region.
        smooth0, smooth1: bool
            whether have smooth transition to the neighboring surfaces.
        dyn0: {None, float}
            If float, sets the slope of y-z curve at the end of section 0, i.e., (dy/dz)|n.
            If None, the slope at section 0 is not specified.
        ratio_end: {float, list of float}
            the ratio controls how fast changing to the original geometry 
            at both ends of the curve (`ip=0, n_point-1`). \n
            If input a list, `ratio_end=[a0, a1, b]`. \n
            If input a float `a`, do not change to the original geometry when `a <= 0`,
            `ratio_end=[a, a, 1]` when `a` > 0. \n
            `a0`, `a1`: larger the faster. \n
            `b`: controls the width of the original geometry at both ends. \n
            `b<=1`: no width, `b>1`: larger the wider. 

        '''
        #* Do not have neighboring surfaces
        if i_sec0 == 0:
            smooth0 = False
        if i_sec1 == self.n_section-1 or i_sec1 == len(self.surfaces)-1:
            smooth1 = False

        n_point = self.surfaces[i_sec0][0].shape[1]

        #* Smoothly change to the original geometry at both ends of the curve (ip=0, n_point-1)
        _x = np.linspace(0, 1, n_point, endpoint=True)
        
        if not isinstance(ratio_end, list):
            if ratio_end <= 0:
                ratio = np.zeros(n_point)
            else:
                ratio = self.smooth_ratio_function(_x, a0=ratio_end, a1=ratio_end, b=1)
        else:
            ratio = self.smooth_ratio_function(_x, a0=ratio_end[0], a1=ratio_end[1], b=ratio_end[2])


        #* For each point in the section curve (n_point)
        for ip in range(n_point):

            #* Collect the spanwise control points
            xx = []
            yy = []
            zz = []
            for i_surf in range(i_sec0, i_sec1):
                xx.append(self.surfaces[i_surf][0][0,ip])
                yy.append(self.surfaces[i_surf][1][0,ip])
                zz.append(self.surfaces[i_surf][2][0,ip])
            xx.append(self.surfaces[i_sec1-1][0][-1,ip])
            yy.append(self.surfaces[i_sec1-1][1][-1,ip])
            zz.append(self.surfaces[i_sec1-1][2][-1,ip])

            #* Construct spanwise spline curve
            bcx0 = (2,0.0)
            bcx1 = (2,0.0)
            bcy0 = (2,0.0)
            bcy1 = (2,0.0)
            
            if smooth0:
                ii = i_sec0-1
                dz = self.surfaces[ii][2][-1,ip] - self.surfaces[ii][2][-2,ip]
                dxz0 = (self.surfaces[ii][0][-1,ip] - self.surfaces[ii][0][-2,ip])/dz
                dyz0 = (self.surfaces[ii][1][-1,ip] - self.surfaces[ii][1][-2,ip])/dz
                bcx0 = (1,dxz0)
                bcy0 = (1,dyz0)
                
            if smooth1:
                ii = i_sec1+1
                dz = self.surfaces[ii][2][1,ip] - self.surfaces[ii][2][0,ip]
                dxz1 = (self.surfaces[ii][0][1,ip] - self.surfaces[ii][0][0,ip])/dz
                dyz1 = (self.surfaces[ii][1][1,ip] - self.surfaces[ii][1][0,ip])/dz
                bcx1 = (1,dxz1)
                bcy1 = (1,dyz1)

            curve_x = CubicSpline(zz, xx, bc_type=(bcx0, bcx1))
            
            if isinstance(dyn0, float) or isinstance(dyn0, int):

                if abs(dyn0)<=1e-6:
                    
                    if ip < n_point-1:
                        _x1 = self.surfaces[i_sec0][0][0,ip+1] - self.surfaces[i_sec0][0][0,ip]
                        _y1 = self.surfaces[i_sec0][1][0,ip+1] - self.surfaces[i_sec0][1][0,ip]
                        _z2 = self.surfaces[i_sec0][2][1,ip]   - self.surfaces[i_sec0][2][0,ip]
                        _x2 = curve_x(self.surfaces[i_sec0][2][1,ip]) - self.surfaces[i_sec0][0][0,ip]
                        _yz = _y1/_z2 * np.clip(_x2/_x1, -1, 1)
                        bcy0 = (1,_yz)
                    else:
                        #* Explicitly use the slope of the previous point (ip-1)
                        _x1 = self.surfaces[i_sec0][0][0,ip]   - self.surfaces[i_sec0][0][0,ip-1]
                        _y1 = self.surfaces[i_sec0][1][0,ip]   - self.surfaces[i_sec0][1][0,ip-1]
                        _z2 = self.surfaces[i_sec0][2][1,ip-1] - self.surfaces[i_sec0][2][0,ip-1]
                        _x2 = curve_x_prev(self.surfaces[i_sec0][2][1,ip-1]) - self.surfaces[i_sec0][0][0,ip-1]
                        _yz = _y1/_z2 * np.clip(_x2/_x1, -1, 1)
                        bcy0 = (1,_yz)
                else:
                    
                    bcy0 = (1,dyn0)

            curve_x_prev = curve_x

            curve_y = CubicSpline(zz, yy, bc_type=(bcy0, bcy1))
            
            #* Use the spanwise spline to update the spanwise geometry
            for i_surf in range(i_sec0, i_sec1):
                ns = self.surfaces[i_surf][0].shape[0]
                for j in range(ns):
                    zi = self.surfaces[i_surf][2][j,ip]
                    self.surfaces[i_surf][0][j,ip] = (1-ratio[ip])*curve_x(zi) + ratio[ip]*self.surfaces[i_surf][0][j,ip]
                    self.surfaces[i_surf][1][j,ip] = (1-ratio[ip])*curve_y(zi) + ratio[ip]*self.surfaces[i_surf][1][j,ip]
   
    def smooth_axisymmetric(self, i_sec0: int, i_sec1: int, phi, linear_TEx=True, 
                            r_TE=None, r_TE_inner=None, func_trans=None):
        '''
        Smooth the axisymmetric curve between i_sec0 and i_sec1
        
        Parameters
        ----------
        i_sec0, i_sec1 : int
            the starting and ending section index of the smooth region
        phi : list or ndarray
            position angle (degree) of control sections: `i_sec0` ~ `i_sec1`
        linear_TEx : bool
            if True, the x coordinates of trailing edge curve are piece-wise linear distribution. 
            Otherwise, they can be nonlinear distribution due to the leading edge curve.
        r_TE : {None, float}
            if None, the trailing edge curve in YZ plane is generated by the layout parameters.
            If provided a float, then the trailing edge curve in YZ plane is set to a circle. 
            Its origin is (0,0), radius is `r_TE`.
        r_TE_inner : {None, float}
            if `r_TE_inner` is provided, it means the control sections are close sections, i.e., 
            both upper and lower surfaces of the control section exist.
            Then, `r_TE_inner` is the inner circle radius.
        func_trans : {None, function}
            if None, `ratio` = `tx`. \n
            If a function `ratio = func_trans(tx)` is provided:
            `ratio` is a float (0~1), representing how much the YZ-plane curve is similar to a circle. \n
            When `ratio` is 1, the curve is the specified circle of which the radius is `r_TE`.
            `tx` is a float (0~1), representing the relative x-axis location of the YZ-plane curve.

        '''
        periodic = False
        if np.abs(phi[0]+phi[-1]-360.0)<1E-3:
            periodic = True

        #* First, smooth the X-axis position of each section
        xx = []
        for i in range(i_sec0, i_sec1+1):
            xx.append(self.sections[i].xLE)

        if periodic:
            curve_x = CubicSpline(phi, xx, bc_type='periodic')
        else:
            curve_x = CubicSpline(phi, xx)

        for i_surf in range(i_sec0, i_sec1):

            sec0 = self.sections[i_surf]
            sec1 = self.sections[i_surf+1]

            for j in range(self.ns):

                tt    = 1.0*j/(self.ns-1.0)
                xLE_  = (1-tt)*sec0.xLE   + tt*sec1.xLE
                chord = (1-tt)*sec0.chord + tt*sec1.chord

                angle = (1-tt)*phi[i_surf] + tt*phi[i_surf+1]
                xLE   = curve_x(angle)  # type: float

                if linear_TEx:
                    self.surfaces[i_surf][0][j,:] = (self.surfaces[i_surf][0][j,:]-xLE_)/chord*(chord-xLE+xLE_) + xLE
                else:
                    self.surfaces[i_surf][0][j,:] += xLE - xLE_

        #* Second, smooth the radius distribution in the circumferential direction
        #  For each point in the section curve (nn)
        nn = self.sections[0].x.shape[0]
        for ip in range(nn):

            # Collect the circumferential control points
            # Must use surfaces data instead of sections data, since only the surfaces data is rotated
            rr = []
            for i_surf in range(i_sec0, i_sec1):
                y_ = self.surfaces[i_surf][1][0,ip]
                z_ = self.surfaces[i_surf][2][0,ip]
                r_ = np.sqrt(y_**2+z_**2)
                rr.append(r_)
            y_ = self.surfaces[i_surf][1][-1,ip]
            z_ = self.surfaces[i_surf][2][-1,ip]
            r_ = np.sqrt(y_**2+z_**2)
            rr.append(r_)

            if periodic:
                curve_r = CubicSpline(phi, rr, bc_type='periodic')
            else:
                curve_r = CubicSpline(phi, rr)
        
            # Use the circumferential spline to update the circumferential geometry
            for i_surf in range(i_sec0, i_sec1):
                for j in range(self.ns):

                    tt    = 1.0*j/(self.ns-1.0)
                    angle = (1-tt)*phi[i_surf] + tt*phi[i_surf+1]
                    R     = curve_r(angle)  # type: float

                    if isinstance(r_TE, float):
                        chord = (1-tt)*self.sections[i_surf].chord + tt*self.sections[i_surf+1].chord
                        xLE_  = (1-tt)*self.sections[i_surf].xLE   + tt*self.sections[i_surf+1].xLE
                        xLE   = curve_x(angle)  # type: float
                        tx = (self.surfaces[i_surf][0][j,ip]-xLE)/(chord-xLE+xLE_)
                        
                        if func_trans is not None:
                            tx = func_trans(tx)

                        if isinstance(r_TE_inner, float):
                            if ip>nn/2.0:
                                R = (1-tx)*R + tx*r_TE
                            else:
                                R = (1-tx)*R + tx*r_TE_inner
                        else:
                            R = (1-tx)*R + tx*r_TE

                    self.surfaces[i_surf][1][j,ip] = R*np.cos(angle/180.0*np.pi)
                    self.surfaces[i_surf][2][j,ip] = R*np.sin(angle/180.0*np.pi)
    
    @staticmethod
    def smooth_ratio_function(x: np.ndarray, a0=4, a1=4, b=1) -> np.ndarray:
        '''
        Smooth ratio function, `x` in [0,1], `ratio` in [0,1].
        
        A larger `a` gives a steeper ramp, 
        a larger `b` gives a longer plateau at both ends.
        '''
        r0 = BasicSurface.smooth_ramp_function(-10*x, a0)
        r1 = BasicSurface.smooth_ramp_function(-10*(1-x), a1)
        ratio = scaled_sigmoid(r0+r1, b)
        return ratio
    
    @staticmethod
    def smooth_ramp_function(x: np.ndarray, a=1) -> np.ndarray:
        '''
        Smooth ramp function, `x`<=0, `y` in [0,1].
        '''
        y1 = 1.0/(1.0+np.exp(-a*x-2))
        y2 = 1.0 + a/4*x
        rr = x < -2/a
        return rr*y1 + (1-rr)*y2


    def bend(self, i_sec0: int, i_sec1: int, leader=None, 
             kx=None, ky=None, kc=None, rot_x=False) -> None:
        '''
        Bend surfaces by a guide curve, i.e., leader.
        
        Parameters
        ------------
        i_sec0, i_sec1 : int
            the index of the start section and the end section.
            
        leader : List[List[float]] or None
            coordinates of the leader curve control points (and chord length). \n
            The leader is a spline curve defined by a list of control points.
            The leading edge point at both ends are automatically included in leader. \n
            `leader = [[x,y,z(,c)], [x,y,z(,c)], ...]`.
            
        kx : {None, float}
            X-axis slope (dx/dz) at both ends [kx0, kx1].
            
        ky : {None, float}
            Y-axis slope (dy/dz) at both ends [ky0, ky1].
            
        kc : {None, float}
            chord slope (dc/dz) at both ends [kc0, kc1].
            
        rot_x : bool
            if True, rotate sections in x-axis to make the section vertical to the leader.

        Notes
        ------
        Regenerate the surface between section i_sec0 and i_sec1. \n
        X is the flow direction (chord direction).

        Examples
        ---------
        >>> bend(i_sec0: int, i_sec1: int, leader=None, 
        >>>         kx=None, ky=None, kc=None, rot_x=False)

        '''
        if self.is_2d:
            print('No bending for 2D cases')
            return

        def sort_z(loc):
            return loc[2]

        #* Control points of the leader curve
        leader_points = []
        spline_chord = False
        if not kc is None:
            spline_chord = True
        elif not leader is None:
            if len(leader[0])==4:
                spline_chord = True

        if spline_chord:
            for i in range(i_sec0, i_sec1+1):
                leader_points.append([self.sections[i].xLE, self.sections[i].yLE, self.sections[i].zLE, self.sections[i].chord])
        else:
            for i in range(i_sec0, i_sec1+1):
                leader_points.append([self.sections[i].xLE, self.sections[i].yLE, self.sections[i].zLE])

        #* Manually provided leader points
        if not leader is None:
            if (spline_chord and len(leader[0])==4) or (not spline_chord and len(leader[0])==3):
                # Need c and provide c // Don't need c and have no c
                for point in leader:
                    leader_points.append(point)
            elif spline_chord and len(leader[0])==3:
                # Need c but have no c
                for point in leader:
                    chord  = self.linear_interpolate_z(point[2], key='chord')
                    point_ = point + [chord]
                    leader_points.append(point_)

            else:
                print('spline_chord', spline_chord)
                print('len(leader[0])', len(leader[0]))
                print('kc', kc)
                raise Exception('Should not happen')

        leader_points.sort(key=sort_z)

        n_point = len(leader_points)

        #* Generating leader curve
        u = np.zeros(n_point)   # independent variable list
        v = np.zeros(n_point)   # dependent variable list
        w = np.zeros(n_point)   # dependent variable list
        c = np.zeros(n_point)   # chord list
        for i in range(n_point):
            u[i] = leader_points[i][2]  # z
            v[i] = leader_points[i][0]  # x
            w[i] = leader_points[i][1]  # y
            if spline_chord:
                c[i] = leader_points[i][3]  # chord
        
        if kx is None:
            leader_x = CubicSpline(u, v)
        else:
            leader_x = CubicSpline(u, v, bc_type=((1,kx[0]), (1,kx[1])))

        if ky is None:
            leader_y = CubicSpline(u, w)
        else:
            leader_y = CubicSpline(u, w, bc_type=((1,ky[0]), (1,ky[1])))

        if spline_chord and kc is None:
            leader_c = CubicSpline(u, c)
        elif not kc is None:
            leader_c = CubicSpline(u, c, bc_type=((1,kc[0]), (1,kc[1])))

        #* Bend surfaces
        i0 = i_sec0
        i1 = i_sec1
        
        for i_surf in range(i0, i1):

            sec0 = self.sections[i_surf]
            sec1 = self.sections[i_surf+1]

            ns = self.surfaces[i_surf][0].shape[0]
            for j in range(ns):

                # Transition of inner sections
                if i_sec0!=0 and j==0:
                    if i_surf==i0:
                        continue

                if i_sec1!=self.n_section-1 and j==ns-1:
                    if i_surf==i1-1:
                        continue

                # Start bending
                xx  = self.surfaces[i_surf][0][j,:]
                yy  = self.surfaces[i_surf][1][j,:]
                zz  = self.surfaces[i_surf][2][j,:]
                nn  = xx.shape[0]

                zLE = zz[0]
                xLE = leader_x(zLE)
                yLE = leader_y(zLE)

                # Original leading edge coordinates
                tt  = 1.0*j/(ns-1.0)
                x0  = (1-tt)*sec0.xLE + tt*sec1.xLE
                y0  = (1-tt)*sec0.yLE + tt*sec1.yLE
                c0  = (1-tt)*sec0.chord + tt*sec1.chord

                #*  Rotation of x-axis (dy/dz)
                if rot_x:
                    angle = -np.arctan(leader_y(zLE, 1))/np.pi*180.0
                    #xx, yy, zz = rotate(xx, yy, zz, angle=angle, origin=[xLE, yLE, zLE])
                    xx, yy, zz = rotate(xx, yy, zz, angle=angle, origin=[x0, y0, zLE])

                #*  Translation
                if spline_chord:
                    xx, _, yy, _ = transform_foil(xx, xx, yy, yy, dx=xLE-x0, dy=yLE-y0, 
                                        x0=xLE, y0=yLE, scale=leader_c(zLE)/c0)

                else:
                    
                    i_half = int(np.floor(nn/2.0))

                    if abs(xx[i_half]-x0)>1e-6 or abs(yy[i_half]-y0)>1e-6:
                        #*  The location of curve end is fixed
                        #   Single piece of open curve to be bent
                        xx, yy = stretch_fixed_point(xx, yy, dx=xLE-x0, dy=yLE-y0, 
                                        xm=x0, ym=y0, xf=xx[-1], yf=yy[-1])

                    else:
                        #*  The locations of the trailing edge of upper and lower surface are fixed
                        #   An airfoil (containing both upper and lower surfaces) to be bent
                        #   Original leading edge:  x0, xu[0], xl[-1]
                        #   New leading edge:       xLE
                        #   Original trailing edge: xu[-1], xl[0]
                        xu = xx[i_half:]
                        xl = xx[:i_half+1]
                        yu = yy[i_half:]
                        yl = yy[:i_half+1]

                        xu, yu = stretch_fixed_point(xu, yu, dx=xLE-x0, dy=yLE-y0, 
                                        xm=x0, ym=y0, xf=xu[-1], yf=yu[-1])

                        xl, yl = stretch_fixed_point(xl, yl, dx=xLE-x0, dy=yLE-y0, 
                                        xm=x0, ym=y0, xf=xl[0], yf=yl[0])

                        xx = np.concatenate((xl,xu[1:]), axis=0)
                        yy = np.concatenate((yl,yu[1:]), axis=0)

                self.surfaces[i_surf][0][j,:] = xx.copy()
                self.surfaces[i_surf][1][j,:] = yy.copy()
                self.surfaces[i_surf][2][j,:] = zz.copy()


    def surf_to_cylinder(self, flip=True, origin=None) -> None:
        '''
        Bend the surface (surfaces) to cylinder (turbomachinery).
        The original surface is constructed by 2D sections.
        
        Parameters
        -----------
        flip : bool
            whether flip `xx` in the reverse order, by default False.
        origin : {None, list of ndarray}
            the cylinder origin axis, by default None.
            If None, the cylinder origin axis is Z-axis for all sections.
            Otherwise, `origin` = [O0, O1, ...], list length is the number of sections.
            Each element is the cylinder origin of that section, i.e., [xOi, yOi] (i=0,1,...),
            it can be a ndarray or list.
        '''

        if origin is None:
            for surf in self.surfaces:
                ns = surf[0].shape[0]
                for j in range(ns):
                    x, y, z = to_cylinder(surf[0][j,:], surf[1][j,:], surf[2][j,:], flip=flip)
                    surf[0][j,:] = x.copy()
                    surf[1][j,:] = y.copy()
                    surf[2][j,:] = z.copy()

            for sec in self.sections:
                sec.x, sec.y, sec.z = to_cylinder(sec.x, sec.y, sec.z, flip=flip)

        else:

            for i in range(len(self.surfaces)):

                surf = self.surfaces[i]
                ns = surf[0].shape[0]

                for j in range(ns):

                    #! This linear interpolation of origins
                    #! causes non-smooth surface even when the smooth function is used
                    tt = j/(ns-1.0)
                    x0 = (1-tt)*origin[i][0] + tt*origin[i+1][0]
                    y0 = (1-tt)*origin[i][1] + tt*origin[i+1][1]

                    x, y, z = to_cylinder(surf[0][j,:], surf[1][j,:], surf[2][j,:], flip=flip, origin=[x0,y0])

                    surf[0][j,:] = x.copy()
                    surf[1][j,:] = y.copy()
                    surf[2][j,:] = z.copy()

            for i in range(self.n_section):
                sec = self.sections[i]
                sec.x, sec.y, sec.z = to_cylinder(sec.x, sec.y, sec.z, flip=flip, origin=origin[i])

    def read_cylinder_origins(self, fname: str) -> None:
        '''
        Read in origins of each section from file.
        
        Parameters
        -------------
        fname: str
            settings file name
        
        Examples
        -----------
        >>> origins = read_cylinder_origins(fname)

        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'CylinderOrigin:': 9}

        origins = []

        found_surf = False
        found_key = 0
        with open(fname, 'r') as f:

            lines = f.readlines()
            iL = 0

            while iL<len(lines):

                line = lines[iL].split()

                if len(line) < 1:
                    iL += 1
                    continue
                
                if not found_surf and len(line) > 1:
                    if '[Surf]' in line[0] and self.name == line[1]:
                        found_surf = True

                elif found_surf and '[Surf]' in line[0]:
                    break

                elif found_surf and found_key == 0:
                    if line[0] in key_dict:
                        found_key = key_dict[line[0]]

                elif found_surf and found_key == 9:
                    for i in range(self.n_section):
                        iL += 1
                        line = lines[iL].split()
                        origins.append([float(line[0]), float(line[1])])

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1

        return origins


    def add_section(self, location: list, axis='Z') -> None:
        '''
        Add sections to the surface, the new sections are interpolated from existed ones.
        
        Parameters
        ---------------
        location : list 
            span-wise locations (must within current sections).
        axis : str
            the direction for interpolation, i.e., Y, Z.

        Notes
        -------
        1. Must run before update_sections(), geo(), geo_axisymmetric() and flip() \n
        2. This will automatically update the curves of all sections \n
        3. X is the flow direction (chord direction) \n

        '''
        if self.is_2d:
            print('Can not add sections in 2D case')
            return

        if len(location) == 0:
            print('Must specify locations when adding sections')
            return
        
        if self.sections[0].xx is None:
            self.update_sections()

        #* Find new section's location
        for loc in location:
            
            found = False
            
            for j in range(self.n_section-1):
                
                if axis == 'Y':
                    if (self.sections[j].yLE-loc)*(self.sections[j+1].yLE-loc)<0.0:
                        rr = (loc - self.sections[j].yLE)/(self.sections[j+1].yLE-self.sections[j].yLE)
                        found = True

                if axis == 'Z':
                    if (self.sections[j].zLE-loc)*(self.sections[j+1].zLE-loc)<0.0:
                        rr = (loc - self.sections[j].zLE)/(self.sections[j+1].zLE-self.sections[j].zLE)
                        found = True
                
                if found:
                    sec_add = interpolate_basic_sections(self.sections[j], self.sections[j+1], ratio=abs(rr))
                    self.sections.insert(j+1, sec_add)
                    break
            
            if not found:
                print('Warning: [Surface.add_section] location %.3f in %s axis is not valid'%(loc, axis))


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

        n_section   = 1 if self.is_2d else self.n_section-1
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
                
                n_point = n_section*(self.ns-1) + 1
                
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

    def output_plot3d(self, fname=None) -> None:
        '''
        Output the surface to `*.xyz` in plot3d format.

        Parameters
        ------------
        fname : str
            name of the output file.
        '''
        
        # Note: X[ns][nn], ns => spanwise
        
        if fname is None:
            fname = self.name + '.xyz'

        n_piece = len(self.surfaces)

        with open(fname, 'w') as f:
            
            f.write('%d \n '%(n_piece))     # Number of surfaces
            
            for i_sec in range(n_piece):
                
                X = self.surfaces[i_sec][0]
                ns = X.shape[0]
                nn = X.shape[1]
                f.write('%d %d 1\n '%(nn, ns))

            for i_sec in range(n_piece):
                
                X = self.surfaces[i_sec][0]
                ns = X.shape[0]
                nn = X.shape[1]
                
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(X[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

                Y = self.surfaces[i_sec][1]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Y[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

                Z = self.surfaces[i_sec][2]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Z[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

    def output_section(self, fname=None, is_2d=True) -> None:
        '''
        Output the control sections.

        Parameters
        ------------
        fname : str
            name of the output file.
        is_2d : bool
            if True, output the 2D unit curves.
            Otherwise, output the 3D control sections.

        '''
        if fname is None:
            fname = self.name + '-section.dat'

        f = open(fname, 'w')

        if is_2d:
            f.write('Variables= X  Y \n ')
            nn = self.sections[0].xx.shape[0]
            for i in range(self.n_section):
                f.write('zone T="sec-u %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f \n'%(self.sections[i].xx[j], self.sections[i].yu[j]))
                f.write('zone T="sec-l %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f \n'%(self.sections[i].xx[j], self.sections[i].yl[j]))

        else:
            f.write('Variables= X  Y  Z \n ')
            nn = self.sections[0].x.shape[0]
            for i in range(self.n_section):
                f.write('zone T="sec %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f  %20.10f \n'%(
                        self.sections[i].x[j], self.sections[i].y[j], self.sections[i].z[j]))

        f.close()

    def plot(self, fig_id=1, type='wireframe', show=True):
        '''
        Plot surface (the figure is not closed).

        Parameters
        ------------
        fig_id : int
            zone_id of the figure
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
            if type == 'wireframe':
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


#* ===========================================
#* CST surfaces
#* ===========================================

class OpenSurface(BasicSurface):
    '''
    Open surface defined by multiple OpenSection objects.
    
    Examples
    ----------
    >>> surf = OpenSurface(n_section=0, name='Patch', nn=1001, ns=101, projection=True)
    '''
    def __init__(self, n_section=0, name='Patch', nn=1001, ns=101, projection=True):

        super().__init__(n_section=n_section, name=name, nn=nn, ns=ns, projection=projection)

        n_ = max(1, n_section)
        self.sections = [ OpenSection() for _ in range(n_) ]

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

        blocks = read_surface_setting_blocks(fname, self.name, self.n_section, key_dict, 1)
        
        if 'Layout:' in blocks:
            apply_surface_layout(self.sections, blocks['Layout:'], self.is_2d)
        
        if 'CST_coefs:' in blocks:
            for i in range(self.n_section):
                self.sections[i].cst = np.array([float(aa) for aa in blocks['CST_coefs:'][i][0]])
        
        if 'CST_refine:' in blocks:
            header, sec_refines = blocks['CST_refine:']
            n_cst_refine = int(header[0])
            i_cst_start  = int(header[1])

            if n_cst_refine > 0:
                for i in range(self.n_section):

                    line1 = sec_refines[i][0]
                    cst_r = np.zeros(n_cst_refine)

                    i1 = 0

                    for j in range(n_cst_refine):
                        if j>=i_cst_start-1 and i1<len(line1):
                            cst_r[j] = float(line1[i1])
                            i1 += 1

                    self.sections[i].refine = cst_r
        
        print('Read surface [%s] settings'%(self.name))

        self.update_layout_center()


class Surface(BasicSurface):
    '''
    Surface defined by multiple Section objects, i.e., foils

    Examples
    -----------
    >>> surf = Surface(n_section=0, name='Wing', nn=1001, ns=101, projection=True)
    '''
    def __init__(self, n_section=0, name='Wing', nn=1001, ns=101, projection=True):

        super().__init__(n_section=n_section, name=name, nn=nn, ns=ns, projection=projection)

        self.sections = [ Section() for _ in range(max(1, n_section)) ]

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

        blocks = read_surface_setting_blocks(fname, self.name, self.n_section, key_dict, 2)
        
        if 'Layout:' in blocks:
            for i, line in enumerate(blocks['Layout:']):

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
        
        if 'CST_coefs:' in blocks:
            for i in range(self.n_section):
                self.sections[i].cst_u = np.array([float(aa) for aa in blocks['CST_coefs:'][i][0]])
                self.sections[i].cst_l = np.array([float(aa) for aa in blocks['CST_coefs:'][i][1]])
        
        if 'CST_refine:' in blocks:
            header, sec_refines = blocks['CST_refine:']
            n_cst_refine = int(header[0])
            i_cst_start  = int(header[1])

            if n_cst_refine > 0:
                for i in range(self.n_section):

                    line1 = sec_refines[i][0]
                    line2 = sec_refines[i][1]

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
        
        print('Read surface [%s] settings'%(self.name))

        # Locate layout center for plot
        self.update_layout_center()

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

        n_section   = 1 if self.is_2d else self.n_section-1
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
                
                npoint = n_section*(self.ns-1) + 1
                ns = self.ns
                nt = int((self.surfaces[0][0].shape[1]+1)/2)
                
                f.write('zone T="sec-u" i= %d j= %d \n'%(nt, npoint))

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

                f.write('zone T="sec-l" i= %d j= %d \n'%(nt, npoint))

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
        Output the surface to `*.xyz` in plot3d format.

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
            fname = self.name + '.xyz'

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

