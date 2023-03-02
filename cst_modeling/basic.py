'''
Basic classes for sections and surfaces, and fundamental functions
'''
import copy
import os
import re
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import spatial
from scipy.interpolate import CubicSpline, interp1d
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation


class BasicSection():
    '''
    Coordinates of the 2D unit curve and the 3D curve.
    
    Parameters
    ----------
    thick : {float, None}, optional
        specified maximum relative thickness, by default None.
    chord : float
        chord length (m), by default 1.
    twist : float
        twist angle, by default 0.
    lTwistAroundLE : bool
        whether the twist center is LE, otherwise TE, by default True.
    
    Examples
    --------
    >>> sec = BasicSection(thick=None, chord=1.0, twist=0.0, lTwistAroundLE=True)
    

    Notes
    -------
    - +x:     flow direction (m)
    - +y:     upside (m)
    - +z:     span-wise (m)
    - twist:  +z direction (deg)
    - chord:  chord length (m)
    - thick:  relative maximum thickness
    - tail:   absolute tail thickness (m)

    
    Attributes
    ------------
    xLE, yLE, zLE : float
        coordinates of the leading edge.
    thick : float
        actual maximum relative thickness.
    chord : float
        chord length.
    twist : float
        twist angle.
    thick_set : {float, None}
        specified maximum relative thickness, by default None.
    lTwistAroundLE : bool
        whether the twist center is LE, otherwise TE.
    xx, yy, yu, yl : ndarray
        the 2D unit curve coordinates. 
        When it is an open section, only `yy` is available.
        When it is a closed section, only `yu` and `yl` are available. 
        They are then concatenated to one curve, e.g., an airfoil.
    x, y, z : ndarray
        the 3D curve coordinates. The 3D curve is generated from the 2D curve through
        translation, scale, and rotation.
    '''

    def __init__(self, thick=None, chord=1.0, twist=0.0, lTwistAroundLE=True) -> None:
        '''
        Create a BasicSection object

        Parameters
        ----------
        thick : {float, None}, optional
            specified maximum relative thickness, by default None.
        chord : float
            chord length (m), by default 1.
        twist : float
            twist angle, by default 0.
        lTwistAroundLE : bool
            whether the twist center is LE, otherwise TE, by default True.
        
        Examples
        --------
        >>> sec = BasicSection(thick=None, chord=1.0, twist=0.0, lTwistAroundLE=True)

        '''
        
        self.xLE = 0.0
        self.yLE = 0.0
        self.zLE = 0.0
        self.chord = chord
        self.twist = twist
        self.thick = 0.0
        self.thick_set = thick
        
        self.lTwistAroundLE = lTwistAroundLE

        #* 2D unit curve
        self.xx  = None
        self.yy  = None     # open curve
        self.yu  = None     # upper surface of closed curve
        self.yl  = None     # lower surface of closed curve

        #* 3D section
        self.x = np.zeros(1)
        self.y = np.zeros(1)
        self.z = np.zeros(1)

    def section(self, nn=1001, flip_x=False, projection=True) -> None:
        '''
        Calculate the 3D curve coordinates from the known 2D curve.

        Parameters
        ------------
        nn : int
            number of points in `xx`, `yy`, `yu`, and `yl`. 
            It's here for the consistency with `Section.section` and `BasicSurface.update_sections`.
        flip_x : bool
            whether flip `xx` in the reverse order, by default False.
        projection : bool
            whether keeps the projection length the same when rotating the section, by default True.
        
        Examples
        ------------
        >>> sec.section(nn=1001, flip_x=False, projection=True)

        '''
        if not isinstance(self.xx, np.ndarray):
            raise Exception('The 2D curve (sec.xx, sec.yy, sec.yu, sec.yl) has not been constructed')

        #* Flip xx
        if flip_x:
            self.xx = np.flip(self.xx)
        
        #* Twist center (rotation after translation and scale)
        if self.lTwistAroundLE:
            xr = None
            yr = None
        else:
            xr = self.xLE + self.chord

        #* Transform to 3D for open section
        if isinstance(self.yy, np.ndarray):
            
            if not self.lTwistAroundLE:
                yr = self.yLE + self.yy[-1]*self.chord
            
            self.x, _, self.y, _ = transform(self.xx, self.xx, self.yy, self.yy, 
                scale=self.chord, rot=self.twist, xr=xr, yr=yr,
                dx=self.xLE, dy=self.yLE, projection=projection)

            self.z = np.ones_like(self.x)*self.zLE

        #* Transform to 3D for closed section
        if isinstance(self.yu, np.ndarray):
            
            if not self.lTwistAroundLE:
                yr = self.yLE + 0.5*(self.yu[-1]+self.yl[-1])*self.chord
            
            xu_, xl_, yu_, yl_ = transform(self.xx, self.xx, self.yu, self.yl, 
                scale=self.chord, rot=self.twist, xr=xr, yr=yr,
                dx=self.xLE, dy=self.yLE, projection=projection)

            self.x = np.concatenate((np.flip(xl_),xu_[1:]), axis=0)
            self.y = np.concatenate((np.flip(yl_),yu_[1:]), axis=0)
            self.z = np.ones_like(self.x)*self.zLE

        if self.x.shape[0] <= 1:
            raise Exception('The 3D curve (sec.x, sec.y, sec.z) is not successfully constructed')


class BasicSurface():
    '''
    Multi-section surface based on `BasicSection`.
    
    Parameters
    ----------
    n_sec : int
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
    >>> surf = BasicSurface(n_sec=1, name='Surf', nn=1001, ns=101, projection=True)
    
    
    Notes
    -------
    - +x:     flow direction (m)
    - +y:     upside (m)
    - +z:     span-wise (m)

        
    Attributes
    ------------
    l2d : bool
        whether this is a 3D surface for a 2D curve (unit span). 
        It is `True` when `n_sec` is 0 or 1.
    secs : list of BasicSection
        section objects.
    surfs : list of list of ndarray
        surface coordinates. List `surfs` contains `n_sec`-1 sub-lists.
        Each sub-list is the coordinates of the 3D surface, which contains 3 `ndarray`.
        The 3 arrays are the X, Y, Z coordinates of the surface.   
        For example, surfs = [[x0, y0, z0], [x1, y1, z1]] when n_sec=3.
    half_span : float
        half span for plotting.
    center : ndarray
        surface center for plotting.

    '''

    def __init__(self, n_sec=1, name='Surf', nn=1001, ns=101, projection=True):
        '''
        Construct a BasicSurface object.
        
        Parameters
        ----------
        n_sec : int
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
        >>> surf = BasicSurface(n_sec=1, name='Surf', nn=1001, ns=101, projection=True)
        '''
        
        n_ = max(1, n_sec)
        self.l2d   = n_ == 1    # type: bool
        self.name  = name       # type: str
        self.nn    = nn         # type: int
        self.ns    = ns         # type: int
        self.secs  = [ BasicSection() for _ in range(n_) ]
        self.surfs = []         # type: list[list]
        self.projection = projection  # type: bool

        # Parameters for plot
        self.half_span = 0.5    # type: float
        self.center = np.array([0.5, 0.5, 0.5])

    @property
    def n_sec(self) -> int:
        '''
        Number of sections
        '''
        return len(self.secs)

    @property
    def zLEs(self) -> List[float]:
        '''
        List of section zLE
        '''
        return [round(sec.zLE,5) for sec in self.secs]

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
        
        key_dict = {'Layout:': 1}

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

                elif found_surf and found_key == 1:
                    for i in range(self.n_sec):
                        iL += 1
                        line = lines[iL].split()
                        self.secs[i].xLE   = float(line[0])
                        self.secs[i].yLE   = float(line[1])
                        self.secs[i].zLE   = float(line[2])
                        self.secs[i].chord = float(line[3])
                        self.secs[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.secs[i].thick_set = float(line[5])

                        if self.l2d:
                            self.secs[i].zLE = 0.0

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1
        
        self.layout_center()

    def layout_center(self) -> None:
        '''
        Locate layout center for plot
        '''
        x_range = [self.secs[0].xLE, self.secs[0].xLE]
        y_range = [self.secs[0].yLE, self.secs[0].yLE]
        z_range = [self.secs[0].zLE, self.secs[0].zLE]
        for i in range(self.n_sec):
            x_range[0] = min(x_range[0], self.secs[i].xLE)
            x_range[1] = max(x_range[1], self.secs[i].xLE+self.secs[i].chord)
            y_range[0] = min(y_range[0], self.secs[i].yLE)
            y_range[1] = max(y_range[1], self.secs[i].yLE)
            z_range[0] = min(z_range[0], self.secs[i].zLE)
            z_range[1] = max(z_range[1], self.secs[i].zLE)
        
        span = np.array([x_range[1]-x_range[0], y_range[1]-y_range[0], z_range[1]-z_range[0]])
        self.half_span = span.max()/2.0
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
        i_sec = self.n_sec
        for i in range(self.n_sec-1):
            if (z-self.secs[i].zLE)*(z-self.secs[i+1].zLE)<0 or z==self.secs[i].zLE:
                i_sec = i

        if i_sec >= self.n_sec:
            raise Exception('z is not within the surface: ', z, self.secs[0].zLE, self.secs[-1].zLE)

        #* Linear interpolation
        tt = (z-self.secs[i_sec].zLE)/(self.secs[i_sec+1].zLE-self.secs[i_sec].zLE)
        key_value = None

        if key == 'x' or key == 'X':
            key_value = (1-tt)*self.secs[i_sec].xLE + tt*self.secs[i_sec+1].xLE
        elif key == 'y' or key == 'Y':
            key_value = (1-tt)*self.secs[i_sec].yLE + tt*self.secs[i_sec+1].yLE
        elif key == 'c' or key == 'C' or  key == 'chord':
            key_value = (1-tt)*self.secs[i_sec].chord + tt*self.secs[i_sec+1].chord
        elif key == 't' or key == 'thick' or key == 'thickness':
            key_value = (1-tt)*self.secs[i_sec].thick + tt*self.secs[i_sec+1].thick
        elif key == 'twist':
            key_value = (1-tt)*self.secs[i_sec].twist + tt*self.secs[i_sec+1].twist
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
        for i in range(self.n_sec):
            self.secs[i].section(nn=self.nn, flip_x=flip_x, projection=self.projection)

    def geo(self, flip_x=False, update_sec=True) -> None:
        '''
        Generate surface geometry.
        
        First, update sections by calling `sec.section()` (optional): \n
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
        self.surfs = []
        
        if self.l2d:
            
            sec_ = copy.deepcopy(self.secs[0])
            sec_.zLE = 1.0
            sec_.z = np.ones_like(sec_.x)
            surf = self.section2surf(self.secs[0], sec_, ns=self.ns)
            self.surfs.append(surf)

        else:
            
            for i in range(self.n_sec-1):
                surf = self.section2surf(self.secs[i], self.secs[i+1], ns=self.ns)
                self.surfs.append(surf)

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
        self.surfs = []

        if self.l2d:
            raise Exception('Axisymmetric geometry can not be 2D surface')

        else:
            for i in range(self.n_sec-1):
                surf = self.section_surf_axisymmetric(self.secs[i], self.secs[i+1], phi[i], phi[i+1], ns=self.ns)
                self.surfs.append(surf)


    @staticmethod
    def section2surf(sec0: BasicSection, sec1: BasicSection, ns=101) -> List[np.ndarray]:
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
        >>> surf = section2surf(sec0, sec1, ns)

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
    def section_surf_axisymmetric(sec0: BasicSection, sec1: BasicSection, 
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
        >>> surf = section_surf_axisymmetric(sec0, sec1, phi0, phi1, ns)
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
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][2]
                    self.surfs[i_sec][2] = copy.deepcopy(self.surfs[i_sec][1])
                    self.surfs[i_sec][1] = copy.deepcopy(temp)

                temp = self.center[2]*1.0
                self.center[2] = self.center[1]*1.0
                self.center[1] = -temp

            if '-X' in axis_:
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][1]
                    self.surfs[i_sec][1] = copy.deepcopy(self.surfs[i_sec][2])
                    self.surfs[i_sec][2] = copy.deepcopy(temp)

                temp = self.center[1]*1.0
                self.center[1] = self.center[2]
                self.center[2] = -temp

            if '+Y' in axis_:
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][0]
                    self.surfs[i_sec][0] = copy.deepcopy(self.surfs[i_sec][2])
                    self.surfs[i_sec][2] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[2]
                self.center[2] = -temp

            if '-Y' in axis_:
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][2]
                    self.surfs[i_sec][2] = copy.deepcopy(self.surfs[i_sec][0])
                    self.surfs[i_sec][0] = copy.deepcopy(temp)

                temp = self.center[2]
                self.center[2] = self.center[0]
                self.center[0] = -temp

            if '+Z' in axis_:
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][1]
                    self.surfs[i_sec][1] = copy.deepcopy(self.surfs[i_sec][0])
                    self.surfs[i_sec][0] = copy.deepcopy(temp)

                temp = self.center[1]
                self.center[1] = self.center[0]
                self.center[0] = -temp

            if '-Z' in axis_:
                for i_sec in range(len(self.surfs)):
                    temp = -self.surfs[i_sec][0]
                    self.surfs[i_sec][0] = copy.deepcopy(self.surfs[i_sec][1])
                    self.surfs[i_sec][1] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[1]
                self.center[1] = -temp

        if 'XY' in plane:
            for i_sec in range(len(self.surfs)):
                self.surfs[i_sec][2] = -self.surfs[i_sec][2]
            self.center[2] = - self.center[2]

        if 'YZ' in plane:
            for i_sec in range(len(self.surfs)):
                self.surfs[i_sec][0] = -self.surfs[i_sec][0]
            self.center[0] = - self.center[0]

        if 'ZX' in plane:
            for i_sec in range(len(self.surfs)):
                self.surfs[i_sec][1] = -self.surfs[i_sec][1]
            self.center[1] = - self.center[1]

    def translate(self, dX=0.0, dY=0.0, dZ=0.0) -> None:
        '''
        Translate surface coordinates
        '''
        for surf in self.surfs:
            surf[0] += dX
            surf[1] += dY
            surf[2] += dZ

        self.center[0] += dX
        self.center[1] += dY
        self.center[2] += dZ

    def scale(self, scale=1.0, X0=0.0, Y0=0.0, Z0=0.0) -> None:
        '''
        Scale surface coordinates about point `[X0, Y0, Z0]`.
        '''
        for surf in self.surfs:
            surf[0] = (surf[0]-X0)*scale + X0
            surf[1] = (surf[1]-Y0)*scale + Y0
            surf[2] = (surf[2]-Z0)*scale + Z0

        self.center[0] = (self.center[0]-X0)*scale + X0
        self.center[1] = (self.center[1]-Y0)*scale + Y0
        self.center[2] = (self.center[2]-Z0)*scale + Z0


    def split(self, ips: list) -> None:
        '''
        Split each surface `surfs` into several pieces.
        Length of `surfs` is multiplied by len(ips)+1.
        
        Parameters
        -----------
        ips : list of int
            split point indices on each section curves

        Notes
        ----------
        `surf` is a list of ndarray, i.e., `[surf_x, surf_y, surf_z]`,
        the shape of each element is `[ns, nn]`.
        '''
        nt = self.surfs[0][0].shape[1]  # i.e., self.nn
        ips.sort()
        ips = [1] + ips + [nt]
        
        surfs = copy.deepcopy(self.surfs)
        self.surfs = []

        for i in range(len(ips)-1):

            for surf in surfs:

                surf_x = surf[0]
                surf_y = surf[1]
                surf_z = surf[2]
                
                self.surfs.append([
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
            The index of surface `surfs` to be smoothed is `i_sec0, ..., i_sec1-1`.
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
        if i_sec1 == self.n_sec-1 or i_sec1 == len(self.surfs)-1:
            smooth1 = False

        n_point = self.surfs[i_sec0][0].shape[1]

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
                xx.append(self.surfs[i_surf][0][0,ip])
                yy.append(self.surfs[i_surf][1][0,ip])
                zz.append(self.surfs[i_surf][2][0,ip])
            xx.append(self.surfs[i_sec1-1][0][-1,ip])
            yy.append(self.surfs[i_sec1-1][1][-1,ip])
            zz.append(self.surfs[i_sec1-1][2][-1,ip])

            #* Construct spanwise spline curve
            bcx0 = (2,0.0)
            bcx1 = (2,0.0)
            bcy0 = (2,0.0)
            bcy1 = (2,0.0)
            
            if smooth0:
                ii = i_sec0-1
                dz = self.surfs[ii][2][-1,ip] - self.surfs[ii][2][-2,ip]
                dxz0 = (self.surfs[ii][0][-1,ip] - self.surfs[ii][0][-2,ip])/dz
                dyz0 = (self.surfs[ii][1][-1,ip] - self.surfs[ii][1][-2,ip])/dz
                bcx0 = (1,dxz0)
                bcy0 = (1,dyz0)
                
            if smooth1:
                ii = i_sec1+1
                dz = self.surfs[ii][2][1,ip] - self.surfs[ii][2][0,ip]
                dxz1 = (self.surfs[ii][0][1,ip] - self.surfs[ii][0][0,ip])/dz
                dyz1 = (self.surfs[ii][1][1,ip] - self.surfs[ii][1][0,ip])/dz
                bcx1 = (1,dxz1)
                bcy1 = (1,dyz1)

            curve_x = CubicSpline(zz, xx, bc_type=(bcx0, bcx1))
            
            if isinstance(dyn0, float) or isinstance(dyn0, int):

                if abs(dyn0)<=1e-6:
                    
                    if ip < n_point-1:
                        _x1 = self.surfs[i_sec0][0][0,ip+1] - self.surfs[i_sec0][0][0,ip]
                        _y1 = self.surfs[i_sec0][1][0,ip+1] - self.surfs[i_sec0][1][0,ip]
                        _z2 = self.surfs[i_sec0][2][1,ip]   - self.surfs[i_sec0][2][0,ip]
                        _x2 = curve_x(self.surfs[i_sec0][2][1,ip]) - self.surfs[i_sec0][0][0,ip]
                        _yz = _y1/_z2 * np.clip(_x2/_x1, -1, 1)
                        bcy0 = (1,_yz)
                    else:
                        bcy0 = (1,_yz)
                else:
                    
                    bcy0 = (1,dyn0)

            curve_y = CubicSpline(zz, yy, bc_type=(bcy0, bcy1))
            
            #* Use the spanwise spline to update the spanwise geometry
            for i_surf in range(i_sec0, i_sec1):
                ns = self.surfs[i_surf][0].shape[0]
                for j in range(ns):
                    zi = self.surfs[i_surf][2][j,ip]
                    self.surfs[i_surf][0][j,ip] = (1-ratio[ip])*curve_x(zi) + ratio[ip]*self.surfs[i_surf][0][j,ip]
                    self.surfs[i_surf][1][j,ip] = (1-ratio[ip])*curve_y(zi) + ratio[ip]*self.surfs[i_surf][1][j,ip]
   
    def smooth_axisymmetric(self, i_sec0: int, i_sec1: int, phi, linear_TEx=True, 
                            RTE=None, RTE_=None, func_trans=None):
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
        RTE : {None, float}
            if None, the trailing edge curve in YZ plane is generated by the layout parameters.
            If provided a float, then the trailing edge curve in YZ plane is set to a circle. 
            Its origin is (0,0), radius is `RTE`.
        RTE_ : {None, float}
            if `RTE_` is provided, it means the control sections are close sections, i.e., 
            both upper and lower surfaces of the control section exist.
            Then, `RTE_` is the inner circle radius.
        func_trans : {None, function}
            if None, `ratio` = `tx`. \n
            If a function `ratio = func_trans(tx)` is provided:
            `ratio` is a float (0~1), representing how much the YZ-plane curve is similar to a circle. \n
            When `ratio` is 1, the curve is the specified circle of which the radius is `RTE`.
            `tx` is a float (0~1), representing the relative x-axis location of the YZ-plane curve.

        '''
        periodic = False
        if np.abs(phi[0]+phi[-1]-360.0)<1E-3:
            periodic = True

        #* First, smooth the X-axis position of each section
        xx = []
        for i in range(i_sec0, i_sec1+1):
            xx.append(self.secs[i].xLE)

        if periodic:
            curve_x = CubicSpline(phi, xx, bc_type='periodic')
        else:
            curve_x = CubicSpline(phi, xx)

        for i_surf in range(i_sec0, i_sec1):

            sec0 = self.secs[i_surf]
            sec1 = self.secs[i_surf+1]

            for j in range(self.ns):

                tt    = 1.0*j/(self.ns-1.0)
                xLE_  = (1-tt)*sec0.xLE   + tt*sec1.xLE
                chord = (1-tt)*sec0.chord + tt*sec1.chord

                angle = (1-tt)*phi[i_surf] + tt*phi[i_surf+1]
                xLE   = curve_x(angle)  # type: float

                if linear_TEx:
                    self.surfs[i_surf][0][j,:] = (self.surfs[i_surf][0][j,:]-xLE_)/chord*(chord-xLE+xLE_) + xLE
                else:
                    self.surfs[i_surf][0][j,:] += xLE - xLE_

        #* Second, smooth the radius distribution in the circumferential direction
        #  For each point in the section curve (nn)
        nn = self.secs[0].x.shape[0]
        for ip in range(nn):

            # Collect the circumferential control points
            # Must use surfs data instead of secs data, since only the surfs data is rotated
            rr = []
            for i_surf in range(i_sec0, i_sec1):
                y_ = self.surfs[i_surf][1][0,ip]
                z_ = self.surfs[i_surf][2][0,ip]
                r_ = np.sqrt(y_**2+z_**2)
                rr.append(r_)
            y_ = self.surfs[i_surf][1][-1,ip]
            z_ = self.surfs[i_surf][2][-1,ip]
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

                    if isinstance(RTE, float):
                        chord = (1-tt)*self.secs[i_surf].chord + tt*self.secs[i_surf+1].chord
                        xLE_  = (1-tt)*self.secs[i_surf].xLE   + tt*self.secs[i_surf+1].xLE
                        xLE   = curve_x(angle)  # type: float
                        tx = (self.surfs[i_surf][0][j,ip]-xLE)/(chord-xLE+xLE_)
                        
                        if func_trans is not None:
                            tx = func_trans(tx)

                        if isinstance(RTE_, float):
                            if ip>nn/2.0:
                                R = (1-tx)*R + tx*RTE
                            else:
                                R = (1-tx)*R + tx*RTE_
                        else:
                            R = (1-tx)*R + tx*RTE

                    self.surfs[i_surf][1][j,ip] = R*np.cos(angle/180.0*np.pi)
                    self.surfs[i_surf][2][j,ip] = R*np.sin(angle/180.0*np.pi)
    
    @staticmethod
    def smooth_ratio_function(x: np.ndarray, a0=4, a1=4, b=1) -> np.ndarray:
        '''
        Smooth ratio function, `x` in [0,1], `ratio` in [0,1].
        
        A larger `a` gives a steeper ramp, 
        a larger `b` gives a longer plateau at both ends.
        '''
        r0 = BasicSurface.smooth_ramp_function(-10*x, a0)
        r1 = BasicSurface.smooth_ramp_function(-10*(1-x), a1)
        ratio = BasicSurface.scaled_sigmoid(r0+r1, b)
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
    
    @staticmethod
    def scaled_sigmoid(x, b=1):
        y = 1.0/(1.0+np.exp(-b*x))
        return (y-np.min(y))/(np.max(y)-np.min(y))


    def bend(self, i_sec0: int, i_sec1: int, leader=None, 
             kx=None, ky=None, kc=None, rot_x=False) -> None:
        '''
        Bend surfaces by a guide curve, i.e., leader.
        
        Parameters
        ------------
        i_sec0, i_sec1 : int
            the index of the start section and the end section.
        leader : {None, list of list}
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
        if self.l2d:
            print('No bending for 2D cases')
            return

        def sortZ(loc):
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
                leader_points.append([self.secs[i].xLE, self.secs[i].yLE, self.secs[i].zLE, self.secs[i].chord])
        else:
            for i in range(i_sec0, i_sec1+1):
                leader_points.append([self.secs[i].xLE, self.secs[i].yLE, self.secs[i].zLE])

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
                    point_ = point.append(chord)
                    leader_points.append(point_)

            else:
                print('spline_chord', spline_chord)
                print('len(leader[0])', len(leader[0]))
                print('kc', kc)
                raise Exception('Should not happen')

        leader_points.sort(key=sortZ)

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

            sec0 = self.secs[i_surf]
            sec1 = self.secs[i_surf+1]

            ns = self.surfs[i_surf][0].shape[0]
            for j in range(ns):

                # Transition of inner sections
                if i_sec0!=0 and j==0:
                    if i_surf==i0:
                        continue

                if i_sec1!=self.n_sec-1 and j==ns-1:
                    if i_surf==i1-1:
                        continue

                # Start bending
                xx  = self.surfs[i_surf][0][j,:]
                yy  = self.surfs[i_surf][1][j,:]
                zz  = self.surfs[i_surf][2][j,:]
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
                    xx, _, yy, _ = transform(xx, xx, yy, yy, dx=xLE-x0, dy=yLE-y0, 
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

                self.surfs[i_surf][0][j,:] = xx.copy()
                self.surfs[i_surf][1][j,:] = yy.copy()
                self.surfs[i_surf][2][j,:] = zz.copy()


    def surf_to_cylinder(self, flip=True, origin=None) -> None:
        '''
        Bend the surface (surfs) to cylinder (turbomachinery).
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
            for surf in self.surfs:
                ns = surf[0].shape[0]
                for j in range(ns):
                    x, y, z = toCylinder(surf[0][j,:], surf[1][j,:], surf[2][j,:], flip=flip)
                    surf[0][j,:] = x.copy()
                    surf[1][j,:] = y.copy()
                    surf[2][j,:] = z.copy()

            for sec in self.secs:
                sec.x, sec.y, sec.z = toCylinder(sec.x, sec.y, sec.z, flip=flip)

        else:

            for i in range(len(self.surfs)):

                surf = self.surfs[i]
                ns = surf[0].shape[0]

                for j in range(ns):

                    #! This linear interpolation of origins
                    #! causes non-smooth surface even when the smooth function is used
                    tt = j/(ns-1.0)
                    x0 = (1-tt)*origin[i][0] + tt*origin[i+1][0]
                    y0 = (1-tt)*origin[i][1] + tt*origin[i+1][1]

                    x, y, z = toCylinder(surf[0][j,:], surf[1][j,:], surf[2][j,:], flip=flip, origin=[x0,y0])

                    surf[0][j,:] = x.copy()
                    surf[1][j,:] = y.copy()
                    surf[2][j,:] = z.copy()

            for i in range(self.n_sec):
                sec = self.secs[i]
                sec.x, sec.y, sec.z = toCylinder(sec.x, sec.y, sec.z, flip=flip, origin=origin[i])

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
                    for i in range(self.n_sec):
                        iL += 1
                        line = lines[iL].split()
                        origins.append([float(line[0]), float(line[1])])

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1

        return origins


    def add_sec(self, location: list, axis='Z') -> None:
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
        if self.l2d:
            print('Can not add sections in 2D case')
            return

        if len(location) == 0:
            print('Must specify locations when adding sections')
            return
        

        #* Find new section's location
        for loc in location:
            
            found = False
            
            for j in range(self.n_sec-1):
                
                if axis in 'Y':
                    if (self.secs[j].yLE-loc)*(self.secs[j+1].yLE-loc)<0.0:
                        rr = (loc - self.secs[j].yLE)/(self.secs[j+1].yLE-self.secs[j].yLE)
                        found = True

                if axis in 'Z':
                    if (self.secs[j].zLE-loc)*(self.secs[j+1].zLE-loc)<0.0:
                        rr = (loc - self.secs[j].zLE)/(self.secs[j+1].zLE-self.secs[j].zLE)
                        found = True
                
                if found:
                    sec_add = interp_basic_sec(self.secs[j], self.secs[j+1], ratio=abs(rr))
                    self.secs.insert(j+1, sec_add)
                    break
            
            if not found:
                print('Warning: [Surface.add_sec] location %.3f in %s axis is not valid'%(loc, axis))


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

        n_sec   = 1 if self.l2d else self.n_sec-1
        n_piece = len(self.surfs)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            ns = self.ns

            if not one_piece:

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfs[i_sec][0]
                    surf_y = self.surfs[i_sec][1]
                    surf_z = self.surfs[i_sec][2]
                    
                    nt = surf_x.shape[1]

                    f.write('zone T="sec %d" i= %d j= %d \n'%(i_sec, nt, ns))

                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))
                            
            else:
                
                n_point = n_sec*(self.ns-1) + 1
                
                f.write('zone T="sec" i= %d j= %d \n'%(nt, n_point))

                for i_sec in range(n_piece):
                    
                    surf_x = self.surfs[i_sec][0]
                    surf_y = self.surfs[i_sec][1]
                    surf_z = self.surfs[i_sec][2]
                    
                    nt = surf_x.shape[1]

                    if i_sec>=n_piece-2:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))

    def output_plot3d(self, fname=None) -> None:
        '''
        Output the surface to `*.grd` in plot3d format.

        Parameters
        ------------
        fname : str
            name of the output file.
        '''
        
        # Note: X[ns][nn], ns => spanwise
        
        if fname is None:
            fname = self.name + '.grd'

        n_piece = len(self.surfs)

        with open(fname, 'w') as f:
            
            f.write('%d \n '%(n_piece))     # Number of surfaces
            
            for i_sec in range(n_piece):
                
                X = self.surfs[i_sec][0]
                ns = X.shape[0]
                nn = X.shape[1]
                f.write('%d %d 1\n '%(nn, ns))

            for i_sec in range(n_piece):
                
                X = self.surfs[i_sec][0]
                ns = X.shape[0]
                nn = X.shape[1]
                
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(X[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

                Y = self.surfs[i_sec][1]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Y[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

                Z = self.surfs[i_sec][2]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Z[i,j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nn-1):
                            f.write(' \n ')

    def output_section(self, fname=None, TwoD=True) -> None:
        '''
        Output the control sections.

        Parameters
        ------------
        fname : str
            name of the output file.
        TwoD : bool
            if True, output the 2D unit curves.
            Otherwise, output the 3D control sections.

        '''
        if fname is None:
            fname = self.name + '-section.dat'

        f = open(fname, 'w')

        if TwoD:
            f.write('Variables= X  Y \n ')
            nn = self.secs[0].xx.shape[0]
            for i in range(self.n_sec):
                f.write('zone T="sec-u %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f \n'%(self.secs[i].xx[j], self.secs[i].yu[j]))
                f.write('zone T="sec-l %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f \n'%(self.secs[i].xx[j], self.secs[i].yl[j]))

        else:
            f.write('Variables= X  Y  Z \n ')
            nn = self.secs[0].x.shape[0]
            for i in range(self.n_sec):
                f.write('zone T="sec %d" i= %d \n'%(i, nn))
                for j in range(nn):
                    f.write('  %20.10f  %20.10f  %20.10f \n'%(
                        self.secs[i].x[j], self.secs[i].y[j], self.secs[i].z[j]))

        f.close()

    def plot(self, fig_id=1, type='wireframe') -> None:
        '''
        Plot surface

        Parameters
        ------------
        fig_id : int
            ID of the figure
        type : str
            'wireframe', or 'surface'
        '''
        fig = plt.figure(fig_id)
        ax = fig.add_subplot(projection='3d')

        for surf in self.surfs:
            if type in 'wireframe':
                ax.plot_wireframe(surf[0], surf[1], surf[2])
            else:
                ax.plot_surface(surf[0], surf[1], surf[2])

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(self.center[0]-self.half_span, self.center[0]+self.half_span)
        ax.set_ylim3d(self.center[1]-self.half_span, self.center[1]+self.half_span)
        ax.set_zlim3d(self.center[2]-self.half_span, self.center[2]+self.half_span)
        plt.show()


#* ===========================================
#* Math
#* ===========================================

def curve_curvature(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    '''
    Calculate curvature of points in the curve
    
    Parameters
    ----------
    x, y: ndarray
        coordinates of the curve
    
    Returns
    --------
    curvature: ndarray
        curvature distribution

    Examples
    -----------
    >>> curvature = curve_curvature(x, y)

    '''
    nn = x.shape[0]
    if nn<3:
        raise Exception('curvature needs at least 3 points')
    
    curvature = np.zeros(nn)
    for i in range(1, nn-1):
        X1 = np.array([x[i-1], y[i-1]])
        X2 = np.array([x[i  ], y[i  ]])
        X3 = np.array([x[i+1], y[i+1]])

        a = np.linalg.norm(X1-X2)
        b = np.linalg.norm(X2-X3)
        c = np.linalg.norm(X3-X1)
        p = 0.5*(a+b+c)
        t = p*(p-a)*(p-b)*(p-c)
        R = a*b*c
        if R <= 1.0E-12:
            curv_ = 0.0
        else:
            curv_ = 4.0*np.sqrt(t)/R

        a1 = X2[0] - X1[0]
        a2 = X2[1] - X1[1]
        b1 = X3[0] - X1[0]
        b2 = X3[1] - X1[1]
        if a1*b2 < a2*b1:
            curv_ = -curv_

        curvature[i] = curv_

    curvature[0] = curvature[1]
    curvature[-1] = curvature[-2]

    return curvature
    
def dis_matrix(xs1: np.ndarray, xs2: np.ndarray) -> np.ndarray:
    '''
    Calculate the distance between vectors in xs1 and xs2.

    Parameters
    ------------
    xs1 : ndarray [n1, nx]
        vectors of all samples.
    xs2 : ndarray [n2, nx]
        vectors of all samples.
    
    Returns
    ---------
    RR : ndarray [n1, n2]
        `dis=sqrt(sum((x1-x2)**2)/nx)`

    Examples
    -----------
    >>> RR = dis_matrix(xs1, xs2)

    Notes
    -----------
    Suggest each components of vectors in x1 and x2 is 0~1.

    '''
    nx = xs1.shape[1]
    RR = cdist(xs1, xs2, metric='euclidean')
    RR = RR/np.sqrt(1.0*nx)
    return RR


#* ===========================================
#* Transformation
#* ===========================================

def transform(xu: np.ndarray, xl: np.ndarray, yu: np.ndarray, yl: np.ndarray, 
              scale=1.0, rot=None, x0=None, y0=None, xr=None, yr=None, dx=0.0, dy=0.0, 
              projection=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
    Apply chord length, twist angle(deg) and leading edge position to a 2D curve.

    Parameters
    -------------
    xu, xl, yu, yl : ndarray
        current 2D curve or unit 2D airfoil.
    scale : bool
        scale factor, e.g., chord length.
    rot : {None, float}
        rotate angle (deg), +z direction for x-y plane, e.g., twist angle.
    x0, y0 : float
        coordinates of the scale center.
    xr, yr : float
        coordinates of the rotation center (rotate after translation and scale).
    dx, dy : float
        translation vector, e.g., leading edge location.
    projection : bool
        whether keeps the projection length the same when rotating the section, by default True.

    Returns
    ---------
    xu_new, xl_new, yu_new, yl_new : ndarray
        coordinates of the new 2D curve.

    Examples
    ---------
    >>> xu_new, xl_new, yu_new, yl_new = transform()

    '''
    #* Translation
    xu_new = dx + xu
    xl_new = dx + xl
    yu_new = dy + yu
    yl_new = dy + yl

    #* Scale center
    if x0 is None:
        x0 = xu_new[0]
    if y0 is None:
        y0 = 0.5*(yu_new[0]+yl_new[0])
    
    #* Scale (keeps the same projection length)
    rr = 1.0
    if projection and not rot is None:
        angle = rot/180.0*np.pi  # rad
        rr = np.cos(angle)

    xu_new = x0 + (xu_new-x0)*scale/rr
    xl_new = x0 + (xl_new-x0)*scale/rr
    yu_new = y0 + (yu_new-y0)*scale/rr
    yl_new = y0 + (yl_new-y0)*scale/rr

    #* Rotation center
    if xr is None:
        xr = x0
    if yr is None:
        yr = y0

    #* Rotation
    if not rot is None:
        xu_new, yu_new, _ = rotate(xu_new, yu_new, None, angle=rot, origin=[xr, yr, 0.0], axis='Z')
        xl_new, yl_new, _ = rotate(xl_new, yl_new, None, angle=rot, origin=[xr, yr, 0.0], axis='Z')

    return xu_new, xl_new, yu_new, yl_new

def rotate(x: np.ndarray, y: np.ndarray, z: np.ndarray,
           angle=0.0, origin=[0.0, 0.0, 0.0], axis='X') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Rotate the 3D curve according to origin
    
    Parameters
    ----------
    x, y, z : ndarray
        coordinates of the curve
    angle : float
        rotation angle (deg)
    origin : list of float
        rotation origin
    axis : {'X', 'Y', 'Z'}
        rotation axis (angle is defined by the right-hand rule along this axis)

    Returns
    --------
    x_, y_, z_ : ndarray
        coordinates of the rotated curve
        
    Examples
    --------
    >>> x_, y_, z_ = rotate(x, y, z, angle=0.0, origin=[0.0, 0.0, 0.0], axis='X')
    
    '''
    cc = np.cos( angle/180.0*np.pi )
    ss = np.sin( angle/180.0*np.pi )
    x_ = copy.deepcopy(x)
    y_ = copy.deepcopy(y)
    z_ = copy.deepcopy(z)

    if axis in 'X':
        y_ = origin[1] + (y-origin[1])*cc - (z-origin[2])*ss
        z_ = origin[2] + (y-origin[1])*ss + (z-origin[2])*cc

    if axis in 'Y':
        z_ = origin[2] + (z-origin[2])*cc - (x-origin[0])*ss
        x_ = origin[0] + (z-origin[2])*ss + (x-origin[0])*cc

    if axis in 'Z':
        x_ = origin[0] + (x-origin[0])*cc - (y-origin[1])*ss
        y_ = origin[1] + (x-origin[0])*ss + (y-origin[1])*cc

    return x_, y_, z_

def stretch_fixed_point(x: np.ndarray, y: np.ndarray, dx=0.0, dy=0.0, 
                        xm=None, ym=None, xf=None, yf=None) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Linearly stretch a 2D curve when a certain point (on the curve) is fixed.

    Parameters
    ------------------
    x, y : ndarray
        coordinates of the 2D curve
    dx, dy : float
        movement of the stretched point
    xm, ym : {None, float}
        coordinates of the stretched point.
        If None, the stretched point is the first element of the curve.
    xf, yf : {None, float}
        coordinates of the fixed point.
        If None, the fixed point is the last element of the curve.

    Returns
    -------------
    x_, y_ : ndarray
        coordinates of the stretched curve
    
    Examples
    ------------
    >>> x_, y_ = stretch_fixed_point(x, y, dx, dy, xm, ym, xf, yf)
    
    '''
    x_ = x.copy()
    y_ = y.copy()

    if xf is None or yf is None:
        xf = x[-1]
        yf = y[-1]

    if xm is None or ym is None:
        xm = x[0]
        ym = y[0]

    lm = np.linalg.norm([xm-xf, ym-yf])

    for i in range(x.shape[0]):
        rr  = np.linalg.norm([x[i]-xf, y[i]-yf]) / lm
        x_[i] = x_[i] + rr*dx
        y_[i] = y_[i] + rr*dy

    return x_, y_

def fromCylinder(x: np.ndarray, y: np.ndarray, z: np.ndarray, 
                 flip=True, origin=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Bend the cylinder curve to a 2D plane curve.

    Parameters
    ----------
    x, y ,z : ndarray
        coordinates of the curve on a cylinder. `x` and `y` must not be 0 at the same time.
    flip : bool
        if True, flip the X of the extracted plane curve.
    origin: {None, list of float}
        if provided a list [x0, y0], the cylinder origin is [x0, y0].

    Returns
    ---------
    X, Y, Z :  ndarray
        coordinates of the curve bent to the 2D X-Y plane

    Notes
    ----------
    The cylinder's default origin is `(0,0,0)`, axis is z-axis.
    
    The origin of cylinder and plane curves is the same (0,0,0).
    
    Cylinder: x, y, z ~~ r, theta, z \n
    Plane:    X, Y, Z \n

    theta = arctan(y/x) \n
    r = sqrt(x^2+y^2) \n
    z = z \n

    X = r*theta \n
    Y = z \n
    Z = r \n
    '''
    coef = -1.0 if flip else 1.0

    if origin is not None:
        x = x - origin[0]
        y = y - origin[1]

    rr = np.sqrt(x*x+y*y)
    tt = np.arctan2(y, x) * coef

    X = rr*tt
    Y = z.copy()
    Z = rr

    return X, Y, Z

def toCylinder(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, 
               flip=True, origin=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Bend the plane sections to curves on a cylinder.

    Parameters
    ----------
    X, Y, Z : ndarray
        coordinates of the curve on a plane. `Z` must not be 0.
    flip : bool
        if True, flip the X of the extracted plane curve.
    origin: {None, list of float}
        if provided a list [x0, y0], the cylinder origin is [x0, y0].

    Returns
    ---------
    x, y ,z :  ndarray
        coordinates of the curve bent to a cylinder.

    Notes
    ----------
    The cylinder's default origin is `(0,0,0)`, axis is z-axis.
    
    The origin of cylinder and plane curves is the same (0,0,0).
    
    Cylinder: x, y, z ~~ r, theta, z \n
    Plane:    X, Y, Z \n

    theta = arctan(y/x) \n
    r = sqrt(x^2+y^2) \n
    z = z \n

    X = r*theta \n
    Y = z \n
    Z = r \n
    '''
    coef = -1.0 if flip else 1.0

    nn = X.shape[0]
    x = np.zeros(nn)
    y = np.zeros(nn)
    z = Y.copy()

    for i in range(nn):
        r = Z[i]
        theta = X[i]/r * coef
        x[i] = r*np.cos(theta)
        y[i] = r*np.sin(theta)

    if origin is not None:
        x = x + origin[0]
        y = y + origin[1]

    return x, y, z

def rotate_vector(x, y, z, angle=0, origin=[0, 0, 0], axis_vector=[0,0,1]) -> np.ndarray:
    '''
    Rotate 3D points (vectors) by axis-angle representation.

    Parameters
    ----------
    x, y, z : float or ndarray [:]
        coordinates of the points.
    angle : float
        rotation angle (deg) about the axis (right-hand rule).
    origin : ndarray [3]
        origin of the rotation axis.
    axis_vector : ndarray [3]
        indicating the direction of an axis of rotation.
        The input `axis_vector` will be normalized to a unit vector `e`.
        The rotation vector, or Euler vector, is `angle*e`.

    Returns
    --------
    points : ndarray [3] or [:,3]
        coordinates of the rotated points
        
    Examples
    --------
    >>> points = rotate_vector(x, y, z, angle=0, origin=[0, 0, 0], axis_vector=[0,0,1])
    
    References
    ----------
    
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html#scipy.spatial.transform.Rotation
    
    https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
    
    https://en.wikipedia.org/wiki/Rotation_matrix
    
    '''
    origin = np.array(origin)
    vector = np.transpose(np.array([x, y, z]))  # [3] or [:,3]
    vector = vector - origin
    
    rotation_vector = np.array(axis_vector)/np.linalg.norm(axis_vector)

    rot = Rotation.from_rotvec(angle*rotation_vector, degrees=True)
    
    # In terms of rotation matricies, this application is the same as rot.as_matrix().dot(vector).
    points = rot.apply(vector) + origin
    
    return points

#* ===========================================
#* Interpolation
#* ===========================================

def interp_basic_sec(sec0: BasicSection, sec1: BasicSection, ratio: float) -> BasicSection:
    '''
    Interpolate a basic section by ratio.
    
    Parameters
    ------------
    sec0, sec1 : BasicSection
        sections at both ends.
    ratio : float
        interpolation ratio.
    
    Examples
    --------------
    >>> sec = interp_basic_sec(sec0, sec1, ratio)
    '''
    
    sec = copy.deepcopy(sec0)

    sec.xLE   = (1-ratio)*sec0.xLE   + ratio*sec1.xLE
    sec.yLE   = (1-ratio)*sec0.yLE   + ratio*sec1.yLE
    sec.zLE   = (1-ratio)*sec0.zLE   + ratio*sec1.zLE
    sec.chord = (1-ratio)*sec0.chord + ratio*sec1.chord
    sec.twist = (1-ratio)*sec0.twist + ratio*sec1.twist
    sec.thick = (1-ratio)*sec0.thick + ratio*sec1.thick

    sec.xx = (1-ratio)*sec0.xx + ratio*sec1.xx
    
    if isinstance(sec.yy, np.ndarray):
        sec.yy = (1-ratio)*sec0.yy + ratio*sec1.yy
    else:
        sec.yu = (1-ratio)*sec0.yu + ratio*sec1.yu
        sec.yl = (1-ratio)*sec0.yl + ratio*sec1.yl

    sec.x  = (1-ratio)*sec0.x + ratio*sec1.x
    sec.y  = (1-ratio)*sec0.y + ratio*sec1.y
    sec.z  = (1-ratio)*sec0.z + ratio*sec1.z

    return sec

def interp_from_curve(x0, x: np.ndarray, y: np.ndarray):
    '''
    Interpolate points from curve represented points [x, y].
    
    Parameters
    ----------
    x0 : float or ndarray
        ndarray/value of x locations to be interpolated.
    x, y : ndarray
        coordinates of the curve.

    Returns
    ----------
    y0 : float or ndarray
        interpolated coordinates

    Examples
    ---------
    >>> y0 = interp_from_curve(x0, x, y)
    '''
    f  = interp1d(x, y, kind='cubic')
    y0 = f(x0)

    return y0

def interpolate_IDW(x0: np.ndarray, xs: np.ndarray, ys: np.ndarray, eps=1e-10) -> np.ndarray:
    '''
    Inverse distance weighted interpolation.
    
    Parameters
    ----------
    x0 : ndarray
        coordinates to be interpolated, shape [n0,3].
    xs : ndarray
        coordinates of the reference points, shape [n, 3].
    ys : ndarray
        values at the reference points, shape [n, ny].
    eps : float
        critical distance between `x0` to `xs`.
        
    Returns
    ----------
    y0 : ndarray
        interpolated values, shape [n0,ny].
    
    Examples
    ---------
    >>> y0 = interpolate_IDW(x0, xs, ys, eps=1e-10)
    
    '''
    n0 = x0.shape[0]
    n  = xs.shape[0]
    ny = ys.shape[1]
    y0 = np.zeros([n0,ny])
    ds = dis_matrix(x0, xs) # [n0, n]
    
    for i0 in range(n0):
        
        if np.min(ds[i0,:]) <= eps:
            j = np.argmin(ds[i0,:])
            y0[i0,:] = ys[j,:]
            continue
        
        ws = ds[i0,:]**-1   # [n]
        w_all = np.sum(ws)  # float
        y0[i0,:] = np.dot(np.transpose(ys), ws)/w_all
        
    return y0

    
#* ===========================================
#* Intersection
#* ===========================================
    
def curve_intersect(x1, y1, x2, y2):
    '''
    Find the intersect index between two curves.
    
    Parameters
    ----------
    x1, y1 : list or ndarray
        coordinates of curve 1. 
    x2, y2 : list or ndarray
        coordinates of curve 2. 
        
    Returns
    ---------
    i1, i2 : int
        index of the closest points in curve 1 & 2.
    points : tuple of ndarray
        tuple of two closest points in curve 1 & 2.
        
    Examples
    -----------
    >>> i1, i2, points = curve_intersect(x1, y1, x2, y2)

    '''

    arr1 = np.vstack((np.array(x1),np.array(y1))).T
    arr2 = np.vstack((np.array(x2),np.array(y2))).T

    tree = spatial.KDTree(arr2)
    distance, arr2_index = tree.query(arr1)
    i1 = distance.argmin()  # type: int
    i2 = arr2_index[i1]     # type: int
    points = (arr1[i1], arr2[i2])

    return i1, i2, points

def intersect_point(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray):
    '''
    Calculate intersection point of two segments p1p2 & p3p4.
    
    Parameters
    --------------
    p1, p2, p3, p4 : ndarray
        ndarray [2] or [:,2]
    
    Returns
    ----------
    pi : ndarray
        ndarray [2] or [:,2]
    '''
    if len(p1.shape)==1:
        a1 = p2[1]-p1[1]
        b1 = p1[0]-p2[0]
        c1 = p1[0]*p2[1]-p2[0]*p1[1]
        a2 = p4[1]-p3[1]
        b2 = p3[0]-p4[0]
        c2 = p3[0]*p4[1]-p4[0]*p3[1]
        dd = a1*b2-a2*b1

        if dd==0:
            print('Parallel segments')
            return None
        else:
            x0 = (c1*b2-c2*b1)/dd
            y0 = (c2*a1-c1*a2)/dd
            return np.array([x0,y0])
    
    else:
        
        a1 = p2[:,1]-p1[:,1]
        b1 = p1[:,0]-p2[:,0]
        c1 = p1[:,0]*p2[:,1]-p2[:,0]*p1[:,1]
        a2 = p4[:,1]-p3[:,1]
        b2 = p3[:,0]-p4[:,0]
        c2 = p3[:,0]*p4[:,1]-p4[:,0]*p3[:,1]
        dd = a1*b2-a2*b1

        if np.any(dd==0):
            print('Parallel segments')
            return None
        else:
            x0 = (c1*b2-c2*b1)/dd
            y0 = (c2*a1-c1*a2)/dd
            pi = np.concatenate((x0, y0), axis=1)
            return pi

def intersect_vec_plane(V0: np.ndarray, V1: np.ndarray, 
                        P0: np.ndarray, P1: np.ndarray,
                        P3: np.ndarray) -> Tuple[np.ndarray, float, float, float]:
    '''
    Calculate the intersection point of a vector and a plane.
    
    Parameters
    -----------
    V0, V1: ndarray [3]
        coordinates of vector: V01.
    P0, P1, P3 : ndarray [3]
        coordinates of three points of plane P0123.
    
    Returns
    -------
    xi : ndarray [3]
        intersection point
    t1, t3 : float
        ratio of xi in P01, P03 direction.
    rv : float
        ratio of xi in V01 direction.
    
    Examples
    ---------
    >>> xi, t1, t3, rv = intersect_vec_plane(V0, V1, P0, P1, P3)

    '''
    nR  = V1 - V0
    l0  = np.linalg.norm(nR) + 1E-20
    nR  = nR / l0
    
    A = np.zeros((3,3))
    A[:,0] = P1-P0
    A[:,1] = P3-P0
    A[:,2] = - nR
    B = V0-P0
    
    Sol = np.linalg.solve(A, B)
    
    t1 = Sol[0]
    t3 = Sol[1]
    rv = Sol[2]/l0
    xi = V0 + nR*Sol[2]

    return xi, t1, t3, rv

def intersect_surface_plane(surface: np.ndarray, P0: np.ndarray, P1: np.ndarray, P3: np.ndarray, 
                            within_bounds=False, original_order=False):
    '''
    Calculate the intersection curve of a surface and a plane.
    
    Parameters
    ----------
    surface : ndarray [ni,nj,3]
        coordinates of surface.
    P0, P1, P3 : ndarray [3]
        coordinates of three points of plane P0123
    within_bounds : bool
        if True, only keep the curve within the bounds of P0123.
    original_order : bool
        if False, rearrange points to form a smooth curve.
    
    Returns
    ---------
    curve : list of ndarray [3]
        intersection curve.
    ij_curve : list of [i,j]
        the index of nearest point in surface to each point of curve.
    xi_curve, yt_curve: ndarray [:]
        relative coordinates in the plane P0123, range in [0,1].
    
    Examples
    ------------
    >>> curve, ij_curve, xi_curve, yt_curve = intersect_surface_plane(surface, P0, P1, P3)

    '''

    ni = surface.shape[0]
    nj = surface.shape[1]
    norm = np.cross(P1-P0, P3-P0)
    norm = norm/np.linalg.norm(norm)
    
    curve = []
    ij_curve = []
    xi_curve = []
    yt_curve = []

    #* To locate points in both sides of the plane
    norm_dis = np.dot(surface-P0, norm) # [ni,nj]
    
    for j in range(nj):
        for i in range(ni):
            
            if i<ni-1:
                if norm_dis[i,j]*norm_dis[i+1,j]<0 or norm_dis[i,j]==0:
                    
                    xi, t1, t3, rv = intersect_vec_plane(surface[i,j,:], surface[i+1,j,:], P0, P1, P3)
                    
                    if rv<=0.0 or rv>=1.0:
                        raise Exception('norm product should guarantee rv in (0,1)')
                    elif within_bounds and (t1<0.0 or t1>1.0 or t3<0.0 or t3>1.0):
                        continue
                    else:
                        ij_curve.append([i,j])
                        curve.append(xi.copy())
                        xi_curve.append(t1)
                        yt_curve.append(t3)
                        continue

            if j<nj-1:
                if norm_dis[i,j]*norm_dis[i,j+1]<0 or norm_dis[i,j]==0:
                    
                    xi, t1, t3, rv = intersect_vec_plane(surface[i,j,:], surface[i,j+1,:], P0, P1, P3)
                    
                    if rv<=0.0 or rv>=1.0:
                        raise Exception('norm product should guarantee rv in (0,1)')
                    elif within_bounds and (t1<0.0 or t1>1.0 or t3<0.0 or t3>1.0):
                        continue
                    else:
                        ij_curve.append([i,j])
                        curve.append(xi.copy())
                        xi_curve.append(t1)
                        yt_curve.append(t3)
                        continue

    #* Rearrange points in correct order
    xi_curve = np.array(xi_curve)
    yt_curve = np.array(yt_curve)
    
    if len(curve)>2 and not original_order:
        _, old_index = rearrange_points(xi_curve, yt_curve)
        
        curve    = [curve[ii]    for ii in old_index]
        ij_curve = [ij_curve[ii] for ii in old_index]
        xi_curve = np.array([xi_curve[ii] for ii in old_index])
        yt_curve = np.array([yt_curve[ii] for ii in old_index])

    return curve, ij_curve, xi_curve, yt_curve

def rearrange_points(xi: np.ndarray, yt: np.ndarray, avg_dir=None, 
                     cri_break=0.02, cri_dup=1e-6) -> Tuple[np.ndarray, List[int]]:
    '''
    Rearrange a list of points in a 2D curve.
    
    Parameters
    ----------
    xi, yt : ndarray [n]
        2D coordinates of the points.
    avg_dir : None or ndarray [2]
        if ndarray, specified average direction, the start point is fixed for the curve.
    cri_break : float
        critical ratio to decide whether the point is the end point of the curve.
    cri_dup : float
        critical distance to drop duplicated points.
    
    Returns
    ----------
    new_curve : ndarray [n,2]
        new curve.
    old_index : list of int
        the index of point in the original curve.
    
    Examples
    ----------
    >>> new_curve, old_index = rearrange_points(xi, yt, avg_dir=None, cri_break=0.1)
    
    Notes
    -------
    There are a few assumptions: 
    
    1. it is an open curve with no intersections
    
    2. most of the points are in the correct (local) order,
    this gives us a average direction of the curve, which can 
    help us find the starting/end point of the curve
    
    3. the next point is its closest point or the closest point in the average direction
    
    4. drop duplicated points 
    
    '''
    indexes = np.arange(0.0, len(xi), 1.0)
    points  = np.array([xi, yt, indexes]).transpose()
    points  = points.copy().tolist()    # [n,3]
    n_point = len(points)
    
    cri_break = max(cri_break, 2.0/n_point)
    
    #* Calculate average direction
    if not isinstance(avg_dir, np.ndarray):
        
        avg_dir = np.zeros(2)
        
        for i in range(len(xi)-1):
            dxi = xi[i+1]-xi[i]
            dyt = yt[i+1]-yt[i]
            ll  = np.sqrt(dxi**2+dyt**2)
            if ll > cri_dup:
                avg_dir += np.array([dxi,dyt])

        la = np.linalg.norm(avg_dir)
        lx = abs(xi[-1]-xi[0])
        ly = abs(yt[-1]-yt[0])

        if la > 0.2*(lx+ly):
            avg_dir = avg_dir/la
        elif lx > ly:
            avg_dir = np.array([1., 0.])
        else:
            avg_dir = np.array([0., 1.])

        ii = np.argmax(np.abs(avg_dir))
        if avg_dir[ii]<0:
            avg_dir = -avg_dir

        fix_start = False
        
    else:
        
        fix_start = True
            
    #* Find the potential start point
    dd = np.dot(np.array(points)[:,:2], avg_dir)
    ii = np.argmin(dd)
    new_curve = [points[ii]]
    points.pop(ii)
    
    #* Get the length scale of the curve
    jj = np.argmax(dd)
    ls = dd[jj]-dd[ii]

    #* Append curve in the average direction
    while len(points)>0:
        
        data = np.array(points)[:,:2]   # [:,2]
        
        # calculate the distance to the last point
        last_point = np.array(new_curve[-1])[None,:2]   # [1,2]
        d2l = np.linalg.norm(data-last_point, axis=1)   # [:]
        i_l = np.argmin(d2l)
        min_dis2last = d2l[i_l]
        
        if min_dis2last<cri_dup:
            points.pop(i_l)
            continue
        
        # calculate the distance to the start point
        start_point = np.array(new_curve[0])[None,:2]   # [1,2]
        d2s = np.linalg.norm(data-start_point, axis=1)  # [:]
        i_s = np.argmin(d2s)
        min_dis2start = d2s[i_s]
        
        if d2s[i_s]<cri_dup:
            points.pop(i_s)
            continue
        
        direction_l = np.dot(data[i_l,:]-last_point,  avg_dir)[0]
        direction_s = np.dot(data[i_s,:]-start_point, avg_dir)[0]
        
        if (min_dis2last<=min_dis2start or fix_start) and (direction_l>0 or min_dis2last<=cri_break*ls):
            # Append to the last point in the average direction
            new_curve.append(points[i_l])
            points.pop(i_l)
            continue
            
        if min_dis2start<=min_dis2last and (direction_s<0 or min_dis2start<=cri_break*ls) and not fix_start:
            # Add before the start point in the opposite of the average direction
            new_curve = [points[i_s]] + new_curve
            points.pop(i_s)
            continue
        
        cri_break = cri_break * 1.1
        
    new_curve = np.array(new_curve)
    old_index = new_curve[:,2].astype(int)

    return new_curve[:,:2], old_index.tolist()

def join_curves(curves: list, cri_dup=1e-6) -> np.ndarray:
    '''
    Join several curves into one piece.
    
    Parameters
    -----------
    curves : list of ndarray [:,3 or 3+nv]
        coordinates and data of the curves.
    cri_dup : float
        critical distance to drop duplicated points.
        
    Returns
    ----------
    new_curve : ndarray [:,3 or 3+nv]
        new curve
    
    Examples
    ----------
    >>> new_curve = join_curves(curves: list, cri_dup=1e-6)
    
    '''
    new_curve = curves[0].copy()    # [:,3]
    
    curves = copy.deepcopy(curves)
    curves.pop(0)
    
    while len(curves)>0:
        
        d00 = []
        d01 = []
        d10 = []
        d11 = []
        
        for cur in curves:
            
            d00.append(np.linalg.norm(new_curve[ 0,:3]-cur[ 0,:3]))
            d01.append(np.linalg.norm(new_curve[ 0,:3]-cur[-1,:3]))
            d10.append(np.linalg.norm(new_curve[-1,:3]-cur[ 0,:3]))
            d11.append(np.linalg.norm(new_curve[-1,:3]-cur[-1,:3]))
            
        min_ds = [np.min(d00), np.min(d01), np.min(d10), np.min(d11)]
        ii_min = np.argmin(min_ds)

        if ii_min == 0:
            
            jj_min    = np.argmin(d00)
            add_curve = curves[jj_min].copy()
            new_curve = np.flip(new_curve, axis=0)
            
        elif ii_min == 1:
            
            jj_min    = np.argmin(d01)
            add_curve = curves[jj_min].copy()
            new_curve = np.flip(new_curve, axis=0)
            add_curve = np.flip(add_curve, axis=0)
            
        elif ii_min == 2:
            
            jj_min    = np.argmin(d10)
            add_curve = curves[jj_min].copy()
            
        elif ii_min == 3:
            
            jj_min    = np.argmin(d11)
            add_curve = curves[jj_min].copy()
            add_curve = np.flip(add_curve, axis=0)
            
        else:
            raise Exception()
        
        if np.min(min_ds)<cri_dup:
            new_curve = np.concatenate((new_curve, add_curve[1:,:]),axis=0)
        else:
            new_curve = np.concatenate((new_curve, add_curve),axis=0)

        curves.pop(jj_min)
        
    return new_curve

def reconstruct_curve_by_length(curve: np.ndarray, n:int) -> np.ndarray:
    '''
    Reconstruct the curve with equidistant points.
    
    Parameters
    ----------
    curve : ndarray [:,3]
        coordinates of the curve
    n : int
        number of points
    
    Returns
    -------------
    new_curve : ndarray [n,3]
        coordinates of the new curve
    '''
    
    #* Parametric curve: x(t), y(t), z(t), t in [0,1]
    n0 = curve.shape[0]
    l0 = 0.0
    tt = np.zeros(n0)
    for i in range(n0-1):
        l0 += np.linalg.norm(curve[i+1,:]-curve[i,:])
        tt[i+1] = l0
    tt = tt/l0
    
    #* Reconstruction
    fx = interp1d(tt, curve[:,0], kind='cubic')
    fy = interp1d(tt, curve[:,1], kind='cubic')
    fz = interp1d(tt, curve[:,2], kind='cubic')
    
    new_curve = np.zeros((n,3))
    
    for i in range(n):
        t = i/(n-1.0)
        new_curve[i,0] = fx(t)
        new_curve[i,1] = fy(t)
        new_curve[i,2] = fz(t)
        
    return new_curve
    
def extract_slice(locations: list, Pref: np.ndarray, dir_norm: np.ndarray, dir_ref=np.array([1.,0.,0.]),
                    fname='surface-aircraft.dat', zone_id=[], index_xyz=[0,1,2], arrange_method='join'):
    '''
    Extract data sliced by planes.
    
    Parameters
    --------------
    locations : list of float
        list of distances to the reference point in the given direction.
    Pref : ndarray [3]
        coordinates of the reference point.
    dir_norm : ndarray [3]
        direction vector normal to the slice plane (will be normalized).
    dir_ref : ndarray [3]
        direction vector that roughly sets the xi-axis in the slice plane.
    fname : str
        file name.
    zone_id : list of int
        index of zones in the tecplot format file, start from 0.
    index_xyz : list of int 
        index of variables in file for X, Y and Z.
    arrange_method : str
        if 'join', keeps the original order of points (suitable for surface with a few blocks).
        If 'rearrange', rearrange points by minimal distance.
    
    Returns
    ------------
    sections : list of ndarray [:,3+nv]
        coordinates and data on the slice.
    name_var : list or str
        name of variables
    
    Examples
    ----------
    >>> sections, name_var = extract_slice(locations, Pref, dir_norm, dir_ref=np.array([1.,0.,0.]),
                    fname='surface-aircraft.dat', zone_id=[], index_xyz=[0,1,2], arrange_method='join')
    '''
    #* Read surface data
    data_, name_var, _ = read_tecplot(fname)
    index_var = [i for i in range(len(name_var))]
    for i in index_xyz:
        index_var.remove(i)
    
    if len(zone_id)==0:
        data = data_
    else:
        data = [data_[i] for i in zone_id]
    
    #* Intersect sections
    dn = dir_norm/np.linalg.norm(dir_norm)
    dr = dir_ref - np.dot(dir_ref, dn)*dn
    dr = dr/np.linalg.norm(dir_norm)
    dt = np.cross(dn, dr)
    dt = dt/np.linalg.norm(dt)
    
    sections = []

    for loc in locations:
        
        P0 = Pref + loc*dn
        P1 = P0 + dr
        P3 = P0 + dt
        
        curves = []
        xi_curves = []
        yt_curves = []
        
        for data_ in data:
            
            surface = np.concatenate((data_[:,:,:,index_xyz[0]:index_xyz[0]+1],
                        data_[:,:,:,index_xyz[1]:index_xyz[1]+1],
                        data_[:,:,:,index_xyz[2]:index_xyz[2]+1]), axis=3)
            surface = surface.squeeze()
            curve, ij_curve, xi_curve, yt_curve = intersect_surface_plane(surface,
                        P0, P1, P3, within_bounds=False, original_order=(arrange_method=='join'))
            
            surface_var = []
            for iv in index_var:
                surface_var.append(data_[:,:,:,iv])
            surface_var = np.transpose(np.array(surface_var), [1,2,3,0]).squeeze()  # [:,:,nv]
            
            if len(curve) == 0:
                continue

            new_curve = []
            for i in range(len(curve)):
                
                ii, jj = ij_curve[i]
                ii = min(surface.shape[0]-2, ii)
                jj = min(surface.shape[1]-2, jj)
                
                xs = [surface[ii,jj,:], surface[ii+1,jj,:], surface[ii,jj+1,:], surface[ii+1,jj+1,:]]
                ys = [surface_var[ii,jj,:], surface_var[ii+1,jj,:], surface_var[ii,jj+1,:], surface_var[ii+1,jj+1,:]]
                
                xyz = curve[i][None,:]
                var = interpolate_IDW(xyz, np.array(xs), np.array(ys))
                tmp = np.concatenate((xyz, var), axis=1).squeeze()

                new_curve.append(tmp)
            
            if arrange_method == 'join':
                curves.append(np.array(new_curve))
            else:
                curves += new_curve
                xi_curves += xi_curve.tolist()
                yt_curves += yt_curve.tolist()

        if arrange_method == 'join':
            curve = join_curves(curves)
        else:
            _, old_index = rearrange_points(np.array(xi_curves), np.array(yt_curves), avg_dir=np.array([1.,0.]))
            curve = np.array([curves[ii] for ii in old_index])

        sections.append(curve.copy())

    return sections, name_var

    
#* ===========================================
#* I/O and format transfer
#* ===========================================

def output_curve(x: np.ndarray, y: np.ndarray, fname='curve.dat', ID=0) -> None:
    '''
    Output airfoil data to tecplot ASCII format file.

    Parameters
    -----------
    x, y : ndarray
        coordinates of the curve.
    ID : int
        if `ID`=0, create new file and write header.
        If `ID`>0, append to existed file.
    '''
    nn = x.shape[0]

    if ID == 0:
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  \n ')

    with open(fname, 'a') as f:
        f.write('zone T="%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            f.write('   %20.9f  %20.9f \n'%(x[i], y[i]))
        f.write('\n')

def output_foil(x: np.ndarray, yu: np.ndarray, yl: np.ndarray, fname='airfoil.dat', ID=0, info=False) -> None:
    '''
    Output airfoil data to tecplot ASCII format file

    Parameters
    -----------
    x, yu, yl : ndarray
        coordinates of the baseline airfoil.
    ID : int
        if `ID`=0, create new file and write header.
        If `ID`>0, append to existed file.
    info: bool
        if True, include curvature, thickness and camber
    '''
    nn = x.shape[0]
    curv_u = np.zeros(nn)
    curv_l = np.zeros(nn)
    camber = np.zeros(nn)
    thickness = np.zeros(nn)
    
    if ID == 0:
        # Write header
        with open(fname, 'w') as f:
            if info: 
                line = 'Variables= X  Y  Curvature Thickness Camber \n '
            else:
                line = 'Variables= X  Y  \n '
            f.write(line)

    if info:
        
        curv_u = curve_curvature(x, yu)
        curv_l = curve_curvature(x, yl)

        thickness = yu-yl
        camber = 0.5*(yu+yl)
        
    with open(fname, 'a') as f:
        f.write('zone T="Upp-%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            line = '   %20.9f  %20.9f'%(x[i], yu[i])
            if info:
                line = line + '  %20.9f  %20.9f  %20.9f'%(curv_u[i], thickness[i], camber[i])
            f.write(line+'\n')
            
        f.write('zone T="Low-%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            line = '   %20.9f  %20.9f'%(x[i], yl[i])
            if info:
                line = line + '  %20.9f  %20.9f  %20.9f'%(curv_l[i], thickness[i], camber[i])
            f.write(line+'\n')

def read_curves(fname='curve.dat'):
    '''
    Read curves from a tecplot format file.
    
    Parameters
    ------------
    fname : str
        file name.
    
    Returns
    -----------
    xs, ys : list of list of float
        coordinates of multiple curves

    Examples
    -----------
    >>> xs, ys = read_curves(fname='curve.dat')

    '''

    xs = []
    ys = []
    with open(fname, 'r') as f:
        lines = f.readlines()

        for line in lines:

            line = line.split()
            if len(line)<=1:
                continue

            if line[0] in 'zone':
                xs.append([])
                ys.append([])
                continue

            if len(line)!=2:
                continue            
            
            xs[-1].append(float(line[0]))
            ys[-1].append(float(line[1]))

    return xs, ys

def read_tecplot(fname='tecplot.dat'):
    '''
    Read a tecplot format data file.
    
    Parameters
    ------------
    fname : str
        file name.
    
    Returns
    -----------
    data : list of ndarray
        data of all zones, shape [ni,nj,nk,nv]. 
    name_var : list of str
        name of variables.
    titles : list of str
        title of zones
    
    Examples
    -------------
    >>> data, name_var, titles = read_tecplot(fname='tecplot.dat')
    

    '''
    name_var = []
    data = []
    titles = []
    n_var = 0
    
    with open(fname, 'r') as f:
        
        lines = f.readlines()
        nLine = len(lines)
        iLine = 0
    
        while iLine < nLine:
            
            line = lines[iLine].split()
            if len(line) == 0:
                iLine += 1
                continue
            
            if line[0] in 'Variables=' or line[0] in 'VARIABLES=' :
                
                line = re.split(r'[=",\s]', lines[iLine])
                while '' in line:
                    line.remove('')

                name_var = line[1:]
                n_var = len(name_var)
                iLine += 1
                continue
        
            if line[0] in 'zone' or line[0] in 'ZONE' or line[0] in 'Zone':
                line = re.split(r'[=\s]', lines[iLine])
                while '' in line:
                    line.remove('')
                
                if 'i' in line:
                    ni = int(line[line.index('i')+1])
                elif 'I' in line:
                    ni = int(line[line.index('I')+1])
                else:
                    ni = 1
                    
                if 'j' in line:
                    nj = int(line[line.index('j')+1])
                elif 'J' in line:
                    nj = int(line[line.index('J')+1])
                else:
                    nj = 1
                    
                if 'k' in line:
                    nk = int(line[line.index('k')+1])
                elif 'K' in line:
                    nk = int(line[line.index('K')+1])
                else:
                    nk = 1
                    
                if 'T' in line:
                    # 非贪婪模式：寻找最短的可能匹配 https://www.cnblogs.com/baxianhua/p/8571967.html
                    str_pat = re.compile(r'\"(.*?)\"')
                    name = str_pat.findall(lines[iLine])
                    titles.append(name[0])
                else:
                    titles.append('')
                    
                data_ = np.zeros((ni,nj,nk,n_var))
                iLine += 1
                
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            
                            line = ['#']
                            while line[0] == '#':
                                line = lines[iLine].split()
                                iLine += 1
                                
                            for v in range(n_var):
                                data_[i,j,k,v] = float(line[v])
                                
                data.append(data_.copy())
                continue

    return data, name_var, titles

def read_block_plot3d(lines: list, iLine0: int, ni: int, nj: int, nk: int) -> Tuple[np.ndarray, int]:
    '''
    Read block data from lines.
    
    Parameters
    -----------
    lines : list of str
        f.readlines() of the entire plot3d formate file
    iLine0 : int
        the first line of this block is lines[iLine0]
    ni, nj, nk: int
        size of this block
        
    Returns
    ---------
    xyz : ndarray
        coordinates, shape `[ni,nj,nk,3]`.
    iLine0_new : int
        index of line after read.

    Examples
    ----------
    >>> xyz, iLine0_new = read_block_plot3d(lines, iLine0, ni, nj, nk)

    '''
    xyz = np.zeros([ni,nj,nk,3])
    ll  = iLine0
    ii  = 0
    line = []

    for m in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):

                    if ii >= len(line)-1:
                        # Need to read the next line
                        line = lines[ll].split()
                        ii = 0
                        ll += 1
                    else:
                        # Read next value
                        ii += 1

                    xyz[i,j,k,m] = float(line[ii])

    iLine0_new = ll

    return xyz, iLine0_new

def output_plot3d(X: list, Y: list, Z: list, fname: str, scale=1.0) -> None:
    '''
    Output surface to fname in plot3d format.
    
    Parameters
    -------------
    X, Y, Z: list of ndarray [ns,nn]
        coordinates
    fname: str
        the name of the file (`*.grd`)

    '''
    # ns: number of spanwise points
    # nn: number of curve points
    
    n0 = len(X)

    with open(fname, 'w') as f:
        f.write('%d \n '%(n0))     # Number of surfaces
        for i_sec in range(n0):
            ns = X[i_sec].shape[0]
            nn = X[i_sec].shape[1]
            f.write('%d %d 1\n '%(nn, ns))

        for i_sec in range(n0):
            ii = 0
            ns = X[i_sec].shape[0]
            nn = X[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(X[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

            ii = 0
            ns = Y[i_sec].shape[0]
            nn = Y[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(Y[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

            ii = 0
            ns = Z[i_sec].shape[0]
            nn = Z[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(Z[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

def plot3d_to_igs(fname='igs'):
    '''
    Converts Plot3d surface grid file [fname.grd] to IGES file [fname.igs].
    
    Original Fortran version by Prof. Zhang Yufei: zhangyufei@tsinghua.edu.cn.
    '''

    #* Read plot3d format file
    if not os.path.exists(fname+'.grd'):
        raise Exception(fname+' does not exist for format transfermation')
    
    with open(fname+'.grd', 'r') as f:
        lines = f.readlines()
        line  = lines[0].split()
        num_block = int(line[0])
        nIJK  = np.zeros([num_block, 5], dtype=int)

        for i in range(num_block):
            line  = lines[i+1].split()
            nIJK[i,0] = int(line[0])
            nIJK[i,1] = int(line[1])
            nIJK[i,2] = int(line[2])
            nIJK[i,3] = idataline(nIJK[i,0], nIJK[i,1])

            if nIJK[i,2]!=1:
                raise Exception('Wrong input file: dimension K is not 1')

            if nIJK[i,0]<4 or nIJK[i,1]<4:
                raise Exception('Wrong input file: dimension I or J less than 4')

        nIJK[0,4] = 1
        for i in range(1, num_block):
            nIJK[i,4] = nIJK[i-1,3] + nIJK[i-1,4]

        kLine = num_block+1

    #* Output IGES format file
    f = open(fname+'.igs', 'w')

    #* Start section and global section
    f.write('This is igs file generated by ZHANG Yufei. All rights reserved.         S      1\n')
    f.write('1H,,1H;,3Higs,7Higs.igs,44HDASSAULT SYSTEMES CATIA V5 R20 - www.3ds.com,G      1\n')
    f.write('27HCATIA Version 5 Release 20 ,32,75,6,75,15,3Higs,1.0,2,2HMM,1000,1.0, G      2\n')
    f.write('15H20180311.223810,0.001,10000.0,5Hyancy,15HDESKTOP-BEPNROH,11,0,15H2018G      3\n')
    f.write('0311.223810,;                                                           G      4\n')

    #* Index section
    iType = 128
    for ib in range(num_block):
        iLineStart = nIJK[ib, 4]
        iLineEnd   = nIJK[ib, 3]

        f.write(' %7d %7d %7d %7d %7d %7d %7d %7d'%(iType, iLineStart, 0, 0, 0, 0, 0, 0))
        f.write(' %1d %1d %1d %1dD %6d\n'%(0, 0, 0, 0, ib*2+1))
        f.write(' %7d %7d %7d %7d %7d'%(iType, 0, 0, iLineEnd, 0))
        if ib<9:
            f.write('                BSp Surf%1d      0D %6d\n'%(ib+1, ib*2+2))
        else:
            f.write('                BSp Surf%2d     0D %6d\n'%(ib+1, ib*2+2))

    #* Data section
    iLine = 0
    for ib in range(num_block):
        ni = nIJK[ib, 0]
        nj = nIJK[ib, 1]
        nk = nIJK[ib, 2]

        # Starting
        iLine += 1
        f.write(' %4d, %4d, %4d, %4d, %4d,'%(iType, ni-1, nj-1, 3, 3))
        f.write(' %4d, %4d, %4d, %4d, %4d, %11dP %6d\n'%(0, 0, 1, 0, 0, ib*2+1, iLine))

        # Node vector
        xKnot = knotx(ni)
        for ix in range(ni+4):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(xKnot[ix], ib*2+1, iLine))
        ximin = xKnot[0]
        ximax = xKnot[-1]

        xKnot = knotx(nj)
        for ix in range(nj+4):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(xKnot[ix], ib*2+1, iLine))
        xjmin = xKnot[0]
        xjmax = xKnot[-1]

        # Node weight
        for j in range(nj):
            for i in range(ni):
                iLine += 1
                f.write('%19.10e, %51dP %6d\n'%(1.0, ib*2+1, iLine))

        # Node coordinates
        xyz, kLine = read_block_plot3d(lines, kLine, ni, nj, nk)
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    iLine += 1
                    f.write('%19.10e,%19.10e,%19.10e,%12dP %6d\n'%(
                        xyz[i,j,k,0], xyz[i,j,k,1], xyz[i,j,k,2],
                        ib*2+1, iLine))

        # Ending
        iLine += 1
        f.write('%14.6e,%14.6e,%14.6e,%14.6e;%12dP %6d\n'%(
            ximin, ximax, xjmin, xjmax, ib*2+1, iLine))

    #* Ending section
    f.write('S %6dG %6dD %6dP %6d %40s %6d\n'%(1, 3, 2*num_block, iLine, 'T', 1))
    f.close()

def idataline(ni: int, nj: int):
    '''
    Function for `plot3d_to_igs`
    '''
    i1 = ni+4
    i2 = nj+4
    i3 = ni*nj
    i4 = ni*nj
    i5 = 1+1

    return i1+i2+i3+i4+i5

def knotx(ni: int) -> np.ndarray:
    '''
    Function for `plot3d_to_igs`.
    
    Returns [0, 0, 0, 0, ...(ni-3)..., 1.0, 1.0, 1.0, 1.0].
    ''' 

    xKnot = np.zeros(ni+4)

    for i in range(ni-3):
        xKnot[i+4] = (i+1.0)/(ni-3.0)

    for i in range(4):
        xKnot[ni+i] = 1.0

    return xKnot


