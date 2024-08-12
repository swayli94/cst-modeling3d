'''
Basic classes for sections and surfaces, and fundamental functions
'''
import copy
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import CubicSpline

from .math import rotate, transform, stretch_fixed_point, toCylinder


class BasicSection():
    '''
    Coordinates of the 2D unit curve (profile) and the 3D curve (section).
    
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
    - chord:  chord length (m)
    - thick:  relative maximum thickness
    - tail:   absolute tail thickness (m)
    - twist:  rotation angle about the z-axis, +z direction (deg)
    - rot_x:  rotation angle about the x-axis, +x direction (deg)
    - rot_y:  rotation angle about the y-axis, +y direction (deg)
    
    Attributes
    ------------
    xLE, yLE, zLE : float
        coordinates of the leading edge.
        
    thick : float
        actual maximum relative thickness.
        
    chord : float
        chord length.
        
    twist : float
        twist angle (degree), i.e., rotation angle about z axis (rot_z).
        
    rot_x, rot_y : float
        rotation angle (degree) about x and y axis.
        
    specified_thickness : float or None
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
        
        self.xLE = 0.0
        self.yLE = 0.0
        self.zLE = 0.0
        self.chord = chord
        self.twist = twist
        self.thick = 0.0
        self.specified_thickness = thick
        
        #* Rotation angles
        # The rotation angles are applied after translation and scale
        # The rotation center is the leading edge (LE) by default
        self.rot_x = 0.0
        self.rot_y = 0.0
        
        self.lTwistAroundLE = lTwistAroundLE

        #* 2D unit curve
        self.xx : np.ndarray = None
        self.yy : np.ndarray = None     # open curve
        self.yu : np.ndarray = None     # upper surface of closed curve
        self.yl : np.ndarray = None     # lower surface of closed curve

        #* 3D section
        self.x : np.ndarray = None
        self.y : np.ndarray = None
        self.z : np.ndarray = None

    @property
    def has_profile(self) -> bool:
        '''
        Whether the 2D profile (unit curve) is constructed
        '''
        return isinstance(self.xx, np.ndarray)
    
    @property
    def has_section(self) -> bool:
        '''
        Whether the 3D section is constructed
        '''
        return isinstance(self.x, np.ndarray)
    
    @property
    def is_open_curve(self) -> bool:
        '''
        Whether the section is an open curve
        '''
        return isinstance(self.yy, np.ndarray)

    @property
    def n_point_profile(self) -> int:
        '''
        Number of points in the 2D profile (unit curve)
        '''
        if isinstance(self.xx, np.ndarray):
            return self.xx.shape[0]
        else:
            return 0
        
    @property
    def n_point_section(self) -> int:
        '''
        Number of points in the 3D section (3D curve)
        '''
        if isinstance(self.x, np.ndarray):
            return self.x.shape[0]
        else:
            return 0

    @property
    def scale(self) -> float:
        '''
        Scale factor, e.g., chord length (m).
        '''
        return self.chord

    @property
    def rot_z(self) -> float:
        '''
        Rotation angle (degree) about z axis, i.e., the twist angle.
        '''
        return self.twist
    
    @rot_z.setter
    def rot_z(self, value: float) -> None:
        self.twist = value

    def get_profile(self) -> List[np.ndarray]:
        '''
        Get the 2D profile (unit curve) coordinates.
        If it is a closed curve, the upper and lower surface curves are combined.
        
        Returns
        ---------
        profile : List[np.ndarray] [2][n_point_section]
            the 2D profile coordinates, [profile_x, profile_y].
        '''
        profile : List[np.ndarray] = [None, None]
        
        if self.xx is None:
            raise Exception('The 2D profile (xx, yy, yu, yl) has not been constructed')
        
        if not self.is_open_curve:
        
            profile[0] = np.concatenate((np.flip(self.xx), self.xx[1:]), axis=0)
            profile[1] = np.concatenate((np.flip(self.yl), self.yu[1:]), axis=0)
            
        else:
            
            profile[0] = self.xx
            profile[1] = self.yy
            
        return profile

    def section(self, flip_x=False, projection=True, nn=None) -> None:
        '''
        Calculate the 3D curve coordinates from the known 2D curve.

        Parameters
        ------------
        flip_x : bool
            whether flip `xx` in the reverse order, by default False.
            
        projection : bool
            whether keeps the projection length the same when rotating the section, by default True.
        
        nn : int
            number of points in `xx`, `yy`, `yu`, and `yl`. 
            It's here for the consistency with `Section.section` and `BasicSurface.update_sections`.
        
        Examples
        ------------
        >>> sec.section(flip_x=False, projection=True)

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
                            self.secs[i].specified_thickness = float(line[5])

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
        x_range = [self.secs[0].xLE, self.secs[0].xLE+self.secs[0].chord]
        y_range = [self.secs[0].yLE, self.secs[0].yLE]
        z_range = [self.secs[0].zLE, self.secs[0].zLE+self.l2d]
        
        if not self.l2d:
            
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
        
        if self.secs[0].xx is None:
            self.update_sections()

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

                    if i_sec>=n_piece-1:
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
        self.layout_center()
        
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

        if show:
            plt.show()
            
        return ax


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

