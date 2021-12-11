'''
Basic classes for sections and surfaces, and fundamental functions
'''
import copy
import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import CubicSpline
from scipy import spatial
from scipy.interpolate import interp1d


class BasicSection():
    '''
    Section: 3D curve and 2D unit curve
    '''
    def __init__(self, thick=None, chord=1.0, twist=0.0):
        self.xLE = 0.0
        self.yLE = 0.0
        self.zLE = 0.0
        self.chord = chord
        self.twist = twist
        self.thick = 0.0
        self.thick_set = thick

        #* 2D unit curve
        self.xx  = None
        self.yy  = None     # open curve
        self.yu  = None     # upper surface of closed curve
        self.yl  = None     # lower surface of closed curve

        #* 3D section
        self.x = np.zeros(1)
        self.y = np.zeros(1)
        self.z = np.zeros(1)

    def set_params(self, init=False, **kwargs):
        '''
        Set parameters of the section

        ### Inputs:
        ```text
        init:   True, set to default values
        ```

        ### kwargs:
        ```text
        xLE, yLE, zLE, chord, twist, thick (None)
        ```
        '''
        if init:
            self.xLE = 0.0
            self.yLE = 0.0
            self.zLE = 0.0
            self.chord = 1.0
            self.twist = 0.0
            self.thick = 0.0
            self.thick_set = None

            return

        if 'xLE' in kwargs.keys():
            self.xLE = kwargs['xLE']

        if 'yLE' in kwargs.keys():
            self.yLE = kwargs['yLE']

        if 'zLE' in kwargs.keys():
            self.zLE = kwargs['zLE']

        if 'chord' in kwargs.keys():
            self.chord = kwargs['chord']

        if 'twist' in kwargs.keys():
            self.twist = kwargs['twist']

        if 'thick' in kwargs.keys():
            self.thick_set = kwargs['thick']

    def section(self, nn=1001, flip_x=False, proj=True):
        '''
        ### Functions:
        ```text
        1. Construct 2D unit curve (null in the BasicSection)
        2. Transform to 3D curve
        ```

        ### Inputs:
        ```text
        nn:     total amount of points (it's here for function BasicSurface.geo_secs)
        flip_x: True ~ flip section.xx in reverse order
        proj:   True => for unit airfoil, the rotation keeps the projection length the same
        ```
        '''
        if not isinstance(self.xx, np.ndarray):
            raise Exception('The 2D curve has not been constructed')

        #* Flip xx
        if flip_x:
            self.xx = np.flip(self.xx)

        #* Transform to 3D for open section
        if isinstance(self.yy, np.ndarray):
            self.x, _, self.y, _ = transform(self.xx, self.xx, self.yy, self.yy, 
                scale=self.chord, rot=self.twist, dx=self.xLE, dy=self.yLE, proj=proj)

            self.z = np.ones_like(self.x)*self.zLE

        #* Transform to 3D for closed section
        if isinstance(self.yu, np.ndarray):
            xu_, xl_, yu_, yl_ = transform(self.xx, self.xx, self.yu, self.yl, 
                scale=self.chord, rot=self.twist, dx=self.xLE, dy=self.yLE, proj=proj)

            self.x = np.concatenate((np.flip(xl_),xu_[1:]), axis=0)
            self.y = np.concatenate((np.flip(yl_),yu_[1:]), axis=0)
            self.z = np.ones_like(self.x)*self.zLE

    def copyfrom(self, other):
        '''
        Copy from anthor BasicSection object
        '''
        if not isinstance(other, BasicSection):
            raise Exception('Must copy from another BasicSection object')
        
        self.xLE = other.xLE
        self.yLE = other.yLE
        self.zLE = other.zLE
        self.chord = other.chord
        self.twist = other.twist

        self.xx = copy.deepcopy(other.xx)
        self.yy = copy.deepcopy(other.yy)
        self.yu = copy.deepcopy(other.yu)
        self.yl = copy.deepcopy(other.yl)

        self.x = other.x.copy()
        self.y = other.y.copy()
        self.z = other.z.copy()


class BasicSurface():
    '''
    Construct multi-section surface with BasicSection objects.

    >>> BasicSurface(n_sec=0, name='Surf', nn=1001, ns=101, project=True)
    '''

    def __init__(self, n_sec=0, name='Surf', nn=1001, ns=101, project=True):

        n_ = max(1, n_sec)
        self.l2d   = n_ == 1    # type: bool
        self.name  = name       # type: str
        self.nn    = nn         # type: int
        self.ns    = ns         # type: int
        self.secs  = [ BasicSection() for _ in range(n_) ]
        self.surfs = []         # type: list[list]
        self.project = project  # type: bool

        # Parameters for plot
        self.half_s = 0.5       # type: float
        self.center = np.array([0.5, 0.5, 0.5])

    @property
    def n_sec(self):
        return len(self.secs)

    @property
    def zLE_secs(self):
        '''
        List of section zLE
        '''
        return [round(sec.zLE,5) for sec in self.secs]

    def read_setting(self, fname: str):
        '''
        Read in Surface layout parameters from file

        ### Inputs:
        ```text
        fname:  control file name
        ```
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

    def layout_center(self):
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
        self.half_s = span.max()/2.0
        self.center[0] = 0.5*(x_range[1]+x_range[0])
        self.center[1] = 0.5*(y_range[1]+y_range[0])
        self.center[2] = 0.5*(z_range[1]+z_range[0])

    def copyfrom(self, other):
        '''
        Copy from another BasicSurface object
        '''
        if not isinstance(other, BasicSurface):
            raise Exception('Must copy from a BasicSurface object')

        self.l2d   = other.l2d
        self.name  = other.name
        self.nn    = other.nn
        self.ns    = other.ns
        self.secs  = copy.deepcopy(other.secs)
        self.surfs = copy.deepcopy(other.surfs)

        self.half_s = other.half_s
        self.center = other.center.copy()

    def linear_interpolate_z(self, z: float, key='x'):
        '''
        Linear interpolation of key by given z

        >>> key_value = linear_interpolate_z(z: float, key='x')

        ### Inputs:
        ```text
        z:      location of the value
        key:    The value to be interpolated
                'x' or 'X'
                'y' or 'Y'
                'c' or 'C' or 'chord'
                't' or 'thick' or 'thickness'
                'twist'
        ```
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


    def geo_secs(self, flip_x=False):
        '''
        Update surface sections

        ### Functions:
        ```text
        1. Construct 2D unit curve (null in the BasicSection)
        2. Transform to 3D curve
        ```

        ### Inputs:
        ```text
        flip_x:     True ~ flip section.xx in reverse order
        ```
        '''
        for i in range(self.n_sec):
            self.secs[i].section(nn=self.nn, flip_x=flip_x, proj=self.project)

    def geo(self, flip_x=False, update_sec=True):
        '''
        Generate surface geometry

        ### Inputs:
        ```text
        flip_x:     True ~ flip section.xx in reverse order
        update_sec: True ~ update sections
        ```
        '''
        if update_sec:
            self.geo_secs(flip_x=flip_x)

        self.surfs = []

        if self.l2d:
            sec_ = copy.deepcopy(self.secs[0])
            sec_.zLE = 1.0
            surf = self.section_surf(self.secs[0], sec_, ns=self.ns)
            self.surfs.append(surf)

        else:
            for i in range(self.n_sec-1):
                surf = self.section_surf(self.secs[i], self.secs[i+1], ns=self.ns)
                self.surfs.append(surf)

    def geo_axisymmetric(self, phi, flip_x=False, update_sec=True):
        '''
        Generate axisymmetric surface geometry

        ### Inputs:
        ```text
        phi:        list or ndarray, position angle of control sections
        flip_x:     True ~ flip section.xx in reverse order
        update_sec: True ~ update sections
        ```
        '''
        if update_sec:
            self.geo_secs(flip_x=flip_x)

        self.surfs = []

        if self.l2d:
            raise Exception('Axisymmetric geometry can not be 2D surface')

        else:
            for i in range(self.n_sec-1):
                surf = self.section_surf_axisymmetric(self.secs[i], self.secs[i+1], phi[i], phi[i+1], ns=self.ns)
                self.surfs.append(surf)


    @staticmethod
    def section_surf(sec0, sec1, ns=101):
        '''
        Interplot surface section between curves

        >>> surf = section_surf(sec0, sec1, ns)

        ### Inputs:
        ```text
        sec0, sec1:     Section object
        ns:             number of spanwise points
        ```

        ### Return: 
        ```text
        surf:   [surf_x, surf_y, surf_z]
                list of ndarray [ns, nn]
        ```
        '''

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
    def section_surf_axisymmetric(sec0, sec1, phi0: float, phi1: float, ns=101):
        '''
        Interplot axisymmetric surface section between curves

        >>> surf = section_surf_axisymmetric(sec0, sec1, ns)

        ### Inputs:
        ```text
        sec0, sec1:     Section object
        phi0, phi1:     angle (degree) about X-axis (X-Y plane is 0 degree)
        ns:             number of spanwise points
        ```

        ### Return: 
        ```text
        surf:   [surf_x, surf_y, surf_z]
                list of ndarray [ns, nn]
        ```
        '''
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


    def flip(self, axis='None', plane='None'):
        '''
        For surfs, and center. (This should be the last action)

        The axis and plane can be a single string, 
        or a string contains multiple actions to take in order, e.g., '+X  +Y'.

        ### Inputs:
        ```text
        axis:  turn 90 degrees about axis: +X, -X, +Y, -Y, +Z, -Z
        plane: get symmetry about plane: 'XY', 'YZ', 'ZX'
        ```
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

    def translate(self, dX=0.0, dY=0.0, dZ=0.0):
        '''
        Translate surface coordinates

        >>> translate(dX=0.0, dY=0.0, dZ=0.0)
        '''
        for surf in self.surfs:
            surf[0] += dX
            surf[1] += dY
            surf[2] += dZ

        self.center[0] += dX
        self.center[1] += dY
        self.center[2] += dZ

    def scale(self, scale=1.0, X0=0.0, Y0=0.0, Z0=0.0):
        '''
        Scale surface coordinates about (X0, Y0, Z0)

        >>> scale(scale=1.0, X0=0.0, Y0=0.0, Z0=0.0)
        '''
        for surf in self.surfs:
            surf[0] = (surf[0]-X0)*scale + X0
            surf[1] = (surf[1]-Y0)*scale + Y0
            surf[2] = (surf[2]-Z0)*scale + Z0

        self.center[0] = (self.center[0]-X0)*scale + X0
        self.center[1] = (self.center[1]-Y0)*scale + Y0
        self.center[2] = (self.center[2]-Z0)*scale + Z0


    def smooth(self, i_sec0: int, i_sec1: int, smooth0=False, smooth1=False, dyn0=None):
        '''
        Smooth the spanwise curve between i_sec0 and i_sec1

        ### Inputs:
        ```text
        i_sec0, i_sec1:     the starting and ending section index of the smooth region
        smooth0, smooth1:   bool, whether have smooth transition to the neighboring surfaces
        dyn0:               (dy/dz)|n, set the slope of y-z curve at the end of section 0
        ```
        '''
        #* Do not have neighboring surfaces
        if i_sec0 == 0:
            smooth0 = False
        if i_sec1 == self.n_sec-1:
            smooth1 = False

        #* For each point in the section curve (n_point)
        n_point = self.surfs[0][0].shape[1]
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
                        _rr = _x2/_x1
                        _yz = _y1/_z2 * np.clip(_x2/_x1, -1, 1)
                        bcy0 = (1,_yz)
                    else:
                        bcy0 = (1,_yz)
                else:
                    
                    bcy0 = (1,dyn0)

            curve_y = CubicSpline(zz, yy, bc_type=(bcy0, bcy1))

            #* Use the spanwise spline to update the spanwise geometry
            for i_surf in range(i_sec0, i_sec1):
                nn = self.surfs[i_surf][0].shape[0]
                for j in range(nn):
                    zi = self.surfs[i_surf][2][j,ip]
                    self.surfs[i_surf][0][j,ip] = curve_x(zi)
                    self.surfs[i_surf][1][j,ip] = curve_y(zi)
    
    def smooth_axisymmetric(self, i_sec0: int, i_sec1: int, phi, linear_TEx=True, RTE=None, RTE_=None, func_trans=None):
        '''
        Smooth the axisymmetric curve between i_sec0 and i_sec1

        ### Inputs:
        ```text
        i_sec0, i_sec1:   the starting and ending section index of the smooth region
        phi:            list or ndarray, position angle of control sections: i_sec0 ~ i_sec1

        linear_TEx:     if True, the x coordinates of trailing edge curve are piece-wise 
                        linear distribution. Otherwise, they can be nonlinear distribution 
                        due to the leading edge curve

        RTE:            default None, then the trailing edge curve in YZ plane is generated
                        by the layout parameters. If provided a float, then the trailing 
                        edge curve in YZ plane is set to a circle. Its origin is (0,0), 
                        radius is RTE
        RTE_:           if RTE_ is provided, it means the control section is close sections
                        i.e., both upper and lower surfaces of the control section exist
                        Then, RTE_ is the inner circle radius

        func_trans:     optional function: ratio = func_trans(tx)
                        ratio is a float (0~1), representing the extent of the YZ-plane 
                        curve being similar to a circle. When ratio is 1, the curve is the 
                        specified circle of which the radius is RTE.
                        tx is a float  (0~1), representing the relative x-axis location of 
                        the YZ-plane curve
                        default None, means ratio = tx
        ```
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
    

    def bend(self, i_sec0: int, i_sec1: int, leader=None, kx=None, ky=None, kc=None, rot_x=False):
        '''
        Bend surfaces by a guide curve, i.e., leader.

        >>> bend(i_sec0: int, i_sec1: int, leader=None, 
        >>>         kx=None, ky=None, kc=None, rot_x=False)

        ### Inputs:
        ```text
        i_sec0:      the index of start section
        i_sec1:      the index of end section
        leader:     list of points (and chord length) in the guide curve. 
                    [[x,y,z(,c)], [x,y,z(,c)], ...]
        axis:       Z-axis, spanwise direction
        kx:         X-axis slope (dx/dz) at both ends [kx0, kx1]
        ky:         Y-axis slope (dy/dz) at both ends [ky0, ky1]
        kc:         Chord  slope (dc/dz) at both ends [kc0, kc1]
        rot_x:      if True, rotate sections in x-axis to 
                    make the section vertical to the leader
        ```

        ### Note:
        ```text
        The leader is a list of points to define the spline curve that 
        describes the leading edge curve. 
        Regenerate the surface between section i_sec0 and i_sec1
        X is the flow direction (chord direction)
        ```
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
                    leader_points.append(point)

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


    def Surf2Cylinder(self, flip=True, origin=None):
        '''
        Bend the surface (surfs) to cylinder (turbomachinery).
        The original surface is constructed by 2D sections.

        ### Inputs:
        ```text
        flip:   if True, flip X
        origin: default None, i.e., the cylinder origin axis is Z-axis for all sections
                otherwise, provide a list of actual cylinder origins, [O0, O1, ...]
                list length is the number of sections
                each element is the cylinder origin of that section, i.e., [xO, yO]
                can be ndarray or list
        ```
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

                    #! This linear interplotation of origins
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

    def read_cylinder_origins(self, fname):
        '''
        Read in orgins of each section from file

        >>> origins = read_cylinder_origins(fname)

        ### Inputs:
        ```text
        fname:  settings file name
        ```
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


    def output_tecplot(self, fname=None, one_piece=False):
        '''
        Output the surface to *.dat in Tecplot format

        ### Inputs:
        ```text
        fname:      the name of the file
        one_piece:  True ~ combine the spanwise sections into one piece
        ```
        '''
        # surf_x[ns,nt], ns => spanwise

        if fname is None:
            fname = self.name + '.dat'

        n_sec   = 1 if self.l2d else self.n_sec-1
        n_piece = len(self.surfs)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            nt = self.surfs[0][0].shape[1]
            ns = self.ns

            if not one_piece:

                for i_sec in range(n_piece):
                    surf_x = self.surfs[i_sec][0]
                    surf_y = self.surfs[i_sec][1]
                    surf_z = self.surfs[i_sec][2]

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

                    if i_sec>=n_piece-2:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))

    def output_plot3d(self, fname=None):
        '''
        Output the surface to *.grd in plot3d format

        ### Inputs:
        ```text
        fname: the name of the file
        ```
        '''
        if fname is None:
            fname = self.name + '.grd'

        n_piece = len(self.surfs)

        # X[ns][nn], ns => spanwise
        X = self.surfs[0][0]
        ns = X.shape[0]
        nn = X.shape[1]
        
        with open(fname, 'w') as f:
            f.write('%d \n '%(n_piece))     # Number of surfaces
            for i_sec in range(n_piece):
                f.write('%d %d 1\n '%(nn, ns))

            for i_sec in range(n_piece):
                X = self.surfs[i_sec][0]
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

    def output_section(self, fname=None, TwoD=True):
        '''
        Output the control sections

        ### Inputs:
        ```text
        fname:  file name of the output file
        TwoD:   if True, output the 2D unit curves
                otherwise, output the 3D control sections
        ```
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

    def plot(self, fig_id=1, type='wireframe'):
        '''
        Plot surface

        ### Inputs:
        ```text
        fig_id: ID of the figure
        type:   wireframe, surface
        ```
        '''
        fig = plt.figure(fig_id)
        ax = Axes3D(fig)

        for surf in self.surfs:
            if type in 'wireframe':
                ax.plot_wireframe(surf[0], surf[1], surf[2])
            else:
                ax.plot_surface(surf[0], surf[1], surf[2])

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(self.center[0]-self.half_s, self.center[0]+self.half_s)
        ax.set_ylim3d(self.center[1]-self.half_s, self.center[1]+self.half_s)
        ax.set_zlim3d(self.center[2]-self.half_s, self.center[2]+self.half_s)
        plt.show()


#* ===========================================
#* Supportive functions
#* ===========================================

def transform(xu, xl, yu, yl, scale=1.0, rot=None, x0=None, y0=None, dx=0.0, dy=0.0, proj=False):
    '''
    Apply chord length, twist angle(deg) and leading edge position to unit airfoil

    >>> xu_new, xl_new, yu_new, yl_new = transform()

    ### Inputs:
    ```text
    xu, xl, yu, yl:  current curve or unit airfoil (ndarray)
    scale:      scale factor, e.g., chord length
    rot:        rotate angle (deg), +z direction for x-y plane, 
                e.g., twist angle
    x0, y0:     rotation and scale center
    dx, dy:     translation, e.g., leading edge location
    proj:       if True, for unit airfoil, the rotation keeps 
                the projection length the same
    ```

    ### Return: 
    ```text
    xu_new, xl_new, yu_new, yl_new (ndarray)
    ```
    '''
    #* Translation
    xu_new = dx + xu
    xl_new = dx + xl
    yu_new = dy + yu
    yl_new = dy + yl

    #* Rotation center
    if x0 is None:
        x0 = xu_new[0]
    if y0 is None:
        y0 = 0.5*(yu_new[0]+yl_new[0])
    
    #* Scale (keeps the same projection length)
    rr = 1.0
    if proj and not rot is None:
        angle = rot/180.0*np.pi  # rad
        rr = np.cos(angle)

    xu_new = x0 + (xu_new-x0)*scale/rr
    xl_new = x0 + (xl_new-x0)*scale/rr
    yu_new = y0 + (yu_new-y0)*scale/rr
    yl_new = y0 + (yl_new-y0)*scale/rr

    #* Rotation
    if not rot is None:
        xu_new, yu_new, _ = rotate(xu_new, yu_new, None, angle=rot, origin=[x0, y0, 0.0], axis='Z')
        xl_new, yl_new, _ = rotate(xl_new, yl_new, None, angle=rot, origin=[x0, y0, 0.0], axis='Z')

    return xu_new, xl_new, yu_new, yl_new

def rotate(x, y, z, angle=0.0, origin=[0.0, 0.0, 0.0], axis='X'):
    '''
    Rotate the 3D curve according to origin

    >>> x_, y_, z_ = rotate(x, y, z, angle, origin, axis)

    ### Inputs:
    ```text
    x,y,z:  curve ndarray
    angle:  rotation angle (deg)
    origin: rotation origin
    axis:   rotation axis (use positive direction to define angle)
    ```

    ### Return:
    x_, y_, z_ (ndarray)
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

def stretch_fixed_point(x, y, dx=0.0, dy=0.0, xm=None, ym=None, xf=None, yf=None):
    '''
    Linearly stretch a curve when certain point is fixed

    >>> x_, y_ = stretch_fixed_point(x, y, dx, dy, xm, ym, xf, yf)

    ### Inputs:
    ```text
    x, y:   curve (ndarray)
    dx, dy: movement of the first element (scaler)
    xm, ym: The point that moves dx, dy (e.g., the first element of the curve)
    xf, yf: The fixed point (e.g., the last element of the curve)
    ```

    ### Returns:
    x_, y_ (ndarray)
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

def interplot_basic_sec(sec0: BasicSection, sec1: BasicSection, ratio: float):
    '''
    Interplot a basic section by ratio.

    >>> sec = interplot_basic_sec(sec0, sec1, ratio)
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


def fromCylinder(x, y, z, flip=True, origin=None):
    '''
    Bend the cylinder curve to a 2D plane curve.

    ### Inputs:
    ```text
    x, y ,z:    ndarray, point coordinates of curves on the cylinder
    flip:       if True, flip the X of plane curve
    origin:     default None.
                if provided a list [x0, y0], then the cylinder origin is [x0, y0]
    ```

    ### Return:
    X, Y, Z:    ndarray, point coordinates of curves bent to 2D X-Y planes

    ### Note:
    ```text
    Cylinder: origin (0,0,0), axis is z-axis
    x and y must not be 0 at the same time

    The origin of cylinder and plane curves is the same (0,0,0).
    
        Cylinder: x, y, z ~~ r, theta, z
        Plane:    X, Y, Z

        theta = arctan(y/x)
        r = sqrt(x^2+y^2)
        z = z

        X = r*theta
        Y = z
        Z = r
    ```
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

def toCylinder(X, Y, Z, flip=True, origin=None):
    '''
    Bend the plane sections to curves on a cylinder.

    ### Inputs:
    ```text
    X, Y, Z:    ndarray, point coordinates of curves on 2D X-Y planes
                Z must not be 0
    flip:       if True, flip the X of plane curve
    origin:     default None.
                if provided a list [x0, y0], then the cylinder origin is [x0, y0]
    ```

    ### Return:
    x, y ,z:    ndarray, point coordinate of curves bent to a cylinder

    ### Note:
    ```text
    The origin of cylinder and plane curves is the same (0,0,0).
    
        Plane:    X, Y, Z
        Cylinder: x, y, z ~~ r, theta, z
        
        theta = arctan(y/x)
        r = sqrt(x^2+y^2)
        z = z

        X = r*theta
        Y = z
        Z = r
    ```
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


def output_curve(x, y, fname='curve.dat', ID=0):
    '''
    Output airfoil data to tecplot ASCII format file

    ### Inputs:
    ```text
    x, y:   current curve (ndarray)
    ID:     >0 append to existed file. 0: write header
    ```
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

def read_curves(fname='curve.dat'):
    '''
    Read curves from a tecplot format file.

    >>> xs, ys = read_curves(fname='curve.dat')

    ### Return:
    ```text
    xs, ys: list [list]
            len(xs) = len(ys) = number of curves
            len(xs[i]) = number of points on curve i
    ```
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


#* ===========================================
#* Intersection and interplotation
#* ===========================================

def interplot_from_curve(x0, x, y) -> np.ndarray:
    '''
    Interplot points from curve represented points [x, y]

    >>> y0 = interplot_from_curve(x0, x, y)

    ### Inputs:
    ```text
    x0  : ndarray/value of x locations to be interploted
    x, y: points of curve (ndarray)
    ```

    ### Return: 
    y0: ndarray/float
    '''
    f  = interp1d(x, y, kind='cubic')
    y0 = f(x0)

    return y0

def curve_intersect(x1, y1, x2, y2):
    '''
    Find the intersect index between two curves.

    >>> i1, i2, points = curve_intersect(x1, y1, x2, y2)

    ### Inputs:
    ```text
    x1, y1: curve 1 coordinates, list or ndarray
    x2, y2: curve 2 coordinates, list or ndarray
    ```

    ### Return:
    ```text
    i1, i2: index of the closest points in curve 1 & 2
    points: tuple of two closest points in curve 1 & 2
    ```
    '''

    arr1 = np.vstack((np.array(x1),np.array(y1))).T
    arr2 = np.vstack((np.array(x2),np.array(y2))).T

    tree = spatial.KDTree(arr2)
    distance, arr2_index = tree.query(arr1)
    i1 = distance.argmin()  # type: int
    i2 = arr2_index[i1]     # type: int
    points = (arr1[i1], arr2[i2])

    return i1, i2, points

def intersect_point(p1, p2, p3, p4):
    '''
    Calculate intersection point of two segments p1p2 & p3p4
    
    ### Inputs:
    ```text
    px: ndarray [2] or [:,2]
    ```
    
    ### Return:
    ```text
    pi: ndarray [2] or [:,2]
    ```
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

def intersect_vec_plane(V0, V1, P0, P1, P3):
    '''
    Calculate the intersection point of a vector and a plane
    
    >>> xi, t1, t3, rv = intersect_vec_plane(V0, V1, P0, P1, P3)
    
    ### Inputs:
    ```text
    V0, V1:     ndarray [3], coordinates of vector: V01
    P0, P1, P3: ndarray [3], coordinates of three points of plane P0123
    ```
    
    ### Return:
    ```text
    xi:     ndarray [3], intersection point
    t1, t3: ratio of xi in P01, P03 direction
    ss:     ratio of xi in V01 direction
    ```
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


#* ===========================================
#* Format transfer
#* ===========================================

def read_block_plot3d(lines, iLine0, ni, nj, nk):
    '''
    Read block data from lines

    >>> xyz, iLine0_new = read_block_plot3d(lines, iLine0, ni, nj, nk)

    ### Inputs:
    ```text
    lines:      f.readlines() of the entire plot3d formate file
    iLine0:     the first line of this block is lines[iLine0]
    ni, nj, nk: size of this block
    ```

    ### Return:
    ```text
    xyz:    ndarray [ni,nj,nk,3]
    ```
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

def output_plot3d(X: list, Y: list, Z: list, fname: str, scale=1.0):
    '''
    Output surface to fname in plot3d format

    ### Inputs:
    ```text
    X, Y, Z:    list of ndarray [ns,nn]
                n0: number of surfaces
                ns: number of spanwise points
                nn: number of curve points
    fname:      the name of the file (*.grd)
    ```
    '''

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
    
    Original Fortran version by Prof. Zhang Yufei: zhangyufei@tsinghua.edu.cn
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

    i1 = ni+4
    i2 = nj+4
    i3 = ni*nj
    i4 = ni*nj
    i5 = 1+1

    return i1+i2+i3+i4+i5

def knotx(ni: int):
    '''
    [0, 0, 0, 0, ...(ni-3)..., 1.0, 1.0, 1.0, 1.0]
    ''' 

    xKnot = np.zeros(ni+4)

    for i in range(ni-3):
        xKnot[i+4] = (i+1.0)/(ni-3.0)

    for i in range(4):
        xKnot[ni+i] = 1.0

    return xKnot


