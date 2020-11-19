'''
This is a module containing functions to construct a surface.
The surface is interploted by sections, e.g., airfoils
'''
import copy
import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import CubicSpline

from cst_modeling.foil import (OpenSection, Section, cst_foil_fit, output_foil,
                               rotate, stretch_fixed_point, toCylinder,
                               transform)


class BasicSurface():
    '''
    Basic functions of Surface classes
    '''

    def __init__(self, n_sec=0, name='Wing', nn=1001, ns=101, project=True):

        n_ = max(1, n_sec)
        self.l2d   = n_ == 1    # type: bool
        self.name  = name       # type: str
        self.nn    = nn         # type: int
        self.ns    = ns         # type: int
        self.secs  = [ Section() for _ in range(n_) ]
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
        raise NotImplementedError

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


    def geo_secs(self, flip_x=False):
        '''
        Update surface sections

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
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][2]
                    self.surfs[isec][2] = copy.deepcopy(self.surfs[isec][1])
                    self.surfs[isec][1] = copy.deepcopy(temp)

                temp = self.center[2]*1.0
                self.center[2] = self.center[1]*1.0
                self.center[1] = -temp

            if '-X' in axis_:
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][1]
                    self.surfs[isec][1] = copy.deepcopy(self.surfs[isec][2])
                    self.surfs[isec][2] = copy.deepcopy(temp)

                temp = self.center[1]*1.0
                self.center[1] = self.center[2]
                self.center[2] = -temp

            if '+Y' in axis_:
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][0]
                    self.surfs[isec][0] = copy.deepcopy(self.surfs[isec][2])
                    self.surfs[isec][2] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[2]
                self.center[2] = -temp

            if '-Y' in axis_:
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][2]
                    self.surfs[isec][2] = copy.deepcopy(self.surfs[isec][0])
                    self.surfs[isec][0] = copy.deepcopy(temp)

                temp = self.center[2]
                self.center[2] = self.center[0]
                self.center[0] = -temp

            if '+Z' in axis_:
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][1]
                    self.surfs[isec][1] = copy.deepcopy(self.surfs[isec][0])
                    self.surfs[isec][0] = copy.deepcopy(temp)

                temp = self.center[1]
                self.center[1] = self.center[0]
                self.center[0] = -temp

            if '-Z' in axis_:
                for isec in range(len(self.surfs)):
                    temp = -self.surfs[isec][0]
                    self.surfs[isec][0] = copy.deepcopy(self.surfs[isec][1])
                    self.surfs[isec][1] = copy.deepcopy(temp)

                temp = self.center[0]
                self.center[0] = self.center[1]
                self.center[1] = -temp

        if 'XY' in plane:
            for isec in range(len(self.surfs)):
                self.surfs[isec][2] = -self.surfs[isec][2]
            self.center[2] = - self.center[2]

        if 'YZ' in plane:
            for isec in range(len(self.surfs)):
                self.surfs[isec][0] = -self.surfs[isec][0]
            self.center[0] = - self.center[0]

        if 'ZX' in plane:
            for isec in range(len(self.surfs)):
                self.surfs[isec][1] = -self.surfs[isec][1]
            self.center[1] = - self.center[1]

    def smooth(self, isec0: int, isec1: int, smooth0=False, smooth1=False):
        '''
        Smooth the spanwise curve between isec0 and isec1

        ### Inputs:
        ```text
        isec0, isec1:       the starting and ending section index of the smooth region
        smooth0, smooth1:   bool, whether have smooth transition to the neighboring surfaces
        ```
        '''
        #* Do not have neighboring surfaces
        if isec0 == 0:
            smooth0 = False
        if isec1 == self.n_sec-1:
            smooth1 = False

        #* For each point in the section curve (npoint)
        npoint = self.surfs[0][0].shape[1]
        for ip in range(npoint):

            #* Collect the spanwise control points
            xx = []
            yy = []
            zz = []
            for i_surf in range(isec0, isec1):
                xx.append(self.surfs[i_surf][0][0,ip])
                yy.append(self.surfs[i_surf][1][0,ip])
                zz.append(self.surfs[i_surf][2][0,ip])
            xx.append(self.surfs[isec1-1][0][-1,ip])
            yy.append(self.surfs[isec1-1][1][-1,ip])
            zz.append(self.surfs[isec1-1][2][-1,ip])

            #* Construct spanwise spline curve
            bcx0 = (2,0.0)
            bcx1 = (2,0.0)
            bcy0 = (2,0.0)
            bcy1 = (2,0.0)
            if smooth0:
                ii = isec0-1
                dz = self.surfs[ii][2][-1,ip] - self.surfs[ii][2][-2,ip]
                dxz0 = (self.surfs[ii][0][-1,ip] - self.surfs[ii][0][-2,ip])/dz
                dyz0 = (self.surfs[ii][1][-1,ip] - self.surfs[ii][1][-2,ip])/dz
                bcx0 = (1,dxz0)
                bcy0 = (1,dyz0)

            if smooth1:
                ii = isec1+1
                dz = self.surfs[ii][2][1,ip] - self.surfs[ii][2][0,ip]
                dxz1 = (self.surfs[ii][0][1,ip] - self.surfs[ii][0][0,ip])/dz
                dyz1 = (self.surfs[ii][1][1,ip] - self.surfs[ii][1][0,ip])/dz
                bcx1 = (1,dxz1)
                bcy1 = (1,dyz1)

            curve_x = CubicSpline(zz, xx, bc_type=(bcx0, bcx1))
            curve_y = CubicSpline(zz, yy, bc_type=(bcy0, bcy1))

            #* Use the spanwise spline to update the spanwise geometry
            for i_surf in range(isec0, isec1):
                nn = self.surfs[i_surf][0].shape[0]
                for j in range(nn):
                    zi = self.surfs[i_surf][2][j,ip]
                    self.surfs[i_surf][0][j,ip] = curve_x(zi)
                    self.surfs[i_surf][1][j,ip] = curve_y(zi)
    
    def smooth_axisymmetric(self, isec0: int, isec1: int, phi, linear_TEx=True, RTE=None, RTE_=None, func_trans=None):
        '''
        Smooth the axisymmetric curve between isec0 and isec1

        ### Inputs:
        ```text
        isec0, isec1:   the starting and ending section index of the smooth region
        phi:            list or ndarray, position angle of control sections: isec0 ~ isec1

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
        for i in range(isec0, isec1+1):
            xx.append(self.secs[i].xLE)

        if periodic:
            curve_x = CubicSpline(phi, xx, bc_type='periodic')
        else:
            curve_x = CubicSpline(phi, xx)

        for i_surf in range(isec0, isec1):

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
            for i_surf in range(isec0, isec1):
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
            for i_surf in range(isec0, isec1):
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
    

    def bend(self, isec0: int, isec1: int, leader=None, kx=None, ky=None, rot_x=False):
        '''
        Bend surfaces by a guide curve, i.e., leader.

        ### Inputs:
        ```text
        isec0:      the index of start section
        isec1:      the index of end section
        leader:     list of points (and chord length) in the guide curve. [[x,y,z(,c)], [x,y,z(,c)]]
        axis:       Z-axis, spanwise direction
        kx:         X-axis slope (dx/dz) at both ends [kx0, kx1]
        ky:         Y-axis slope (dy/dz) at both ends [ky0, ky1]
        rot_x:      True ~ rotate sections in x-axis to make the section vertical to the leader
        ```

        ### Note:
        ```text
        The leader is a list of points to define a spline curve describing the leading edge curve.
        Regenerate the surface between section isec0 and isec1
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
        if not leader is None:
            if len(leader[0]) == 4:
                # chord length specified
                spline_chord = True
                for i in range(isec0, isec1+1):
                    leader_points.append([self.secs[i].xLE, self.secs[i].yLE, self.secs[i].zLE, self.secs[i].chord])
            else:
                for i in range(isec0, isec1+1):
                    leader_points.append([self.secs[i].xLE, self.secs[i].yLE, self.secs[i].zLE])
            for point in leader:
                leader_points.append(point)
        else:
            for i in range(isec0, isec1+1):
                leader_points.append([self.secs[i].xLE, self.secs[i].yLE, self.secs[i].zLE])

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

        if spline_chord:
            leader_c = CubicSpline(u, c)

        #* Bend surfaces
        i0 = isec0
        i1 = isec1

        for i_surf in range(i0, i1):

            sec0 = self.secs[i_surf]
            sec1 = self.secs[i_surf+1]

            ns = self.surfs[i_surf][0].shape[0]
            for j in range(ns):

                # Transition of inner sections
                if isec0!=0 and j==0:
                    if i_surf==i0:
                        continue

                if isec1!=self.n_sec-1 and j==ns-1:
                    if i_surf==i1-1:
                        continue

                # Start bending
                xx  = self.surfs[i_surf][0][j,:]
                yy  = self.surfs[i_surf][1][j,:]
                zz  = self.surfs[i_surf][2][j,:]

                zLE = zz[0]
                xLE = leader_x(zLE)
                yLE = leader_y(zLE)

                tt  = 1.0*j/(ns-1.0)
                x0  = (1-tt)*sec0.xLE + tt*sec1.xLE
                y0  = (1-tt)*sec0.yLE + tt*sec1.yLE

                # Translation
                c0  = (1-tt)*sec0.chord + tt*sec1.chord
                if spline_chord:
                    xx, _, yy, _ = transform(xx, xx, yy, yy, dx=xLE-x0, dy=yLE-y0, x0=xLE, y0=yLE, scale=leader_c(zLE)/c0)
                else:
                    # The location of trailing edge (xTE, yTE) is fixed
                    xTE = xx[-1]
                    yTE = yy[-1]
                    xx, yy = stretch_fixed_point(xx, yy, dx=xLE-x0, dy=yLE-y0, xm=x0, ym=y0, xf=xTE, yf=yTE )

                # Rotation of x-axis (dy/dz)
                if rot_x:
                    angle = -np.arctan(leader_y(zLE, 1))/np.pi*180.0
                    xx, yy, zz = rotate(xx, yy, zz, angle=angle, origin=[xLE, yLE, zLE])

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

            for i in range(self.n_sec-1):

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

                for isec in range(n_piece):
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    f.write('zone T="sec %d" i= %d j= %d \n'%(isec, nt, ns))

                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j], surf_y[i,j], surf_z[i,j]))
                            
            else:
                
                npoint = n_sec*(self.ns-1) + 1
                
                f.write('zone T="sec" i= %d j= %d \n'%(nt, npoint))

                for isec in range(n_piece):
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    if isec>=n_piece-2:
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

        n_sec   = 1 if self.l2d else self.n_sec-1
        n_piece = len(self.surfs)

        # X[ns][nn], ns => spanwise
        X = self.surfs[0][0]
        ns = X.shape[0]
        nn = X.shape[1]
        
        with open(fname, 'w') as f:
            f.write('%d \n '%(n_piece))     # Number of surfaces
            for isec in range(n_piece):
                f.write('%d %d 1\n '%(nn, ns))

            for isec in range(n_piece):
                X = self.surfs[isec][0]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(X[i,j]))
                        ii += 1
                        if ii%3 == 0:
                            f.write(' \n ')

                Y = self.surfs[isec][1]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Y[i,j]))
                        ii += 1
                        if ii%3 == 0:
                            f.write(' \n ')

                Z = self.surfs[isec][2]
                ii = 0
                for i in range(ns):
                    for j in range(nn):
                        f.write(' %.9f '%(Z[i,j]))
                        ii += 1
                        if ii%3 == 0:
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


class OpenSurface(BasicSurface):
    '''
    Open surface defined by multiple OpenSection objects

    >>> OpenSurface(n_sec=0, name='Patch', nn=1001, ns=101, project=True)
    '''
    def __init__(self, n_sec=0, name='Patch', nn=1001, ns=101, project=True):

        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, project=project)

        n_ = max(1, n_sec)
        self.secs = [ OpenSection() for _ in range(n_) ]

    def read_setting(self, fname):
        '''
        Read in Surface layout and CST parameters from file

        ### Inputs:
        ```text
        fname:  settings file name
        ```
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3, 'CST_flip:': 4}

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
                            self.secs[i].thick = float(line[5])

                        if self.l2d:
                            self.secs[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_sec):
                        iL += 2
                        line = lines[iL].split()
                        self.secs[i].cst = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()
                        cst_r = np.zeros(n_cst_refine)

                        i1 = 0

                        for j in range(n_cst_refine):
                            if j>=i_cst_start-1 and i1<len(line1):
                                cst_r[j] = float(line1[i1])
                                i1 += 1

                        self.secs[i].set_params(refine=cst_r)

                    found_key = 0

                elif found_surf and found_key == 4:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()
                        cst_r = np.zeros(n_cst_refine)

                        i1 = 0
                        for j in range(n_cst_refine):
                            if i1<len(line1):
                                cst_r[j] = float(line1[i1])
                                i1 += 1

                        self.secs[i].set_params(cst_flip=cst_r)

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1
        
        print('Read surface [%s] settings'%(self.name))

        self.layout_center()


class Surface(BasicSurface):
    '''
    Surface defined by multiple Section objects, i.e., foils

    >>> Surface(n_sec=0, name='Wing', nn=1001, ns=101, project=True)

    ### Inputs:
    ```text
    n_sec:   number of control sections (2D if set to 0 or 1)
    name:    name of the surface
    nn:      number of points of upper/lower section
    ns:      number of spanwise points
    project: True ~ projected chord length does not change when twisted
    ```

    ### Note:
    ```text
    +x:     flow direction (m)
    +y:     upside (m)
    +z:     spanwise (m)
    twist:  +z direction (deg)
    chord:  chord length (m)
    thick:  relative maximum thickness
    tail:   absolute tail thickness (m)
    ```

    ### Attributes:
    ```text
    secs:   list of [Section] class
    surfs:  list of [surf_x, surf_y, surf_z], they are [ns, nn] ndarray
    ```
    '''
    def __init__(self, n_sec=0, name='Wing', nn=1001, ns=101, project=True):
        '''
        Initialize the CST surface (upper & lower)
        '''
        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, project=project)

        n_ = max(1, n_sec)
        self.secs = [ Section() for _ in range(n_) ]

    def read_setting(self, fname, tail=0.0):
        '''
        Read in Surface layout and CST parameters from file

        ### Inputs:
        ```text
        fname:  settings file name
        tail:   float or list, tail thickness (m) of each section
        ```
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3, 'CST_flip:': 4}

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
                            self.secs[i].thick = float(line[5])

                        if isinstance(tail, float):
                            self.secs[i].tail  = tail/self.secs[i].chord
                        elif len(tail)==self.n_sec:
                            self.secs[i].tail  = tail[i]/self.secs[i].chord
                        else:
                            raise Exception('tail must be a float or a list with length = section number')
                        
                        if self.secs[i].thick <= 0.0:
                            self.secs[i].thick = None

                        if self.l2d:
                            self.secs[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_sec):
                        iL += 2
                        line = lines[iL].split()
                        self.secs[i].cst_u = np.array([float(aa) for aa in line])

                        iL += 1
                        line = lines[iL].split()
                        self.secs[i].cst_l = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()

                        iL += 1
                        line2 = lines[iL].split()

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

                        self.secs[i].set_params(
                            refine_fixed_t=True,
                            refine_u=cst_ur,
                            refine_l=cst_lr,
                        )

                    found_key = 0

                elif found_surf and found_key == 4:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()

                        iL += 1
                        line2 = lines[iL].split()

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

                        self.secs[i].set_params(
                            cst_flip_u=cst_ur,
                            cst_flip_l=cst_lr,
                        )

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1
        
        print('Read surface [%s] settings'%(self.name))

        # Locate layout center for plot
        self.layout_center()

    def add_sec(self, location: list, axis='Z'):
        '''
        Add sections to the surface, the new sections are interpolated from current ones

        ### Inputs:
        ```text
        location: list of spanwise location (must within current sections)
        axis:     the direction for interplotation Y,Z
        ```

        ### Note:   
        ```text
        Must run before geo(), geo_secs() and flip()
        This will automatically update the curves of all sections
        X is the flow direction (chord direction)
        ```
        '''
        if self.l2d:
            print('Can not add sections in 2D case')
            return

        if len(location) == 0:
            print('Must specify locations when adding sections')
            return

        #* First update current sections
        self.geo_secs()

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
                    sec_add = interplot_sec(self.secs[j], self.secs[j+1], ratio=abs(rr))
                    self.secs.insert(j+1, sec_add)
                    break

    def output_tecplot(self, fname=None, one_piece=False, split=False):
        '''
        Output the surface to *.dat in Tecplot format

        ### Inputs:
        ```text
        fname:      the name of the file
        one_piece:  if True, combine the spanwise sections into one piece
        split:      if True, split to upper and lower surfaces
        ```
        '''
        if not split:
            super().output_tecplot(fname=fname, one_piece=one_piece)
            return

        if fname is None:
            fname = self.name + '.dat'

        n_sec   = 1 if self.l2d else self.n_sec-1
        n_piece = len(self.surfs)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            if not one_piece:

                for isec in range(n_piece):
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    # surf_x[ns,nt], ns => spanwise
                    ns = self.ns
                    nt = self.nn

                    f.write('zone T="sec-u %d" i= %d j= %d \n'%(isec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                    f.write('zone T="sec-l %d" i= %d j= %d \n'%(isec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))

            else:
                
                npoint = n_sec*(self.ns-1) + 1
                ns = self.ns
                nt = self.nn
            
                f.write('zone T="sec-u" i= %d j= %d \n'%(nt, npoint))

                for isec in range(n_piece):
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    if isec>=n_piece-2:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                f.write('zone T="sec-l" i= %d j= %d \n'%(nt, npoint))

                for isec in range(n_piece):
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    if isec>=n_piece-2:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))


#* ===========================================
#* Static functions
#* ===========================================
def interplot_sec(sec0: Section, sec1: Section, ratio: float):
    '''
    Interplot a section by ratio. CST coefficients are gained by cst_foil_fit.

    >>> sec = interplot_sec(sec0, sec1, ratio)
    '''
    
    sec = Section()
    sec.copyfrom(sec0)

    sec.xLE   = (1-ratio)*sec0.xLE   + ratio*sec1.xLE
    sec.yLE   = (1-ratio)*sec0.yLE   + ratio*sec1.yLE
    sec.zLE   = (1-ratio)*sec0.zLE   + ratio*sec1.zLE
    sec.chord = (1-ratio)*sec0.chord + ratio*sec1.chord
    sec.twist = (1-ratio)*sec0.twist + ratio*sec1.twist
    sec.thick = (1-ratio)*sec0.thick + ratio*sec1.thick
    sec.tail  = (1-ratio)*sec0.tail  + ratio*sec1.tail
    sec.RLE   = (1-ratio)*sec0.RLE   + ratio*sec1.RLE

    sec.xx = (1-ratio)*sec0.xx + ratio*sec1.xx
    sec.yu = (1-ratio)*sec0.yu + ratio*sec1.yu
    sec.yl = (1-ratio)*sec0.yl + ratio*sec1.yl

    sec.x  = (1-ratio)*sec0.x + ratio*sec1.x
    sec.y  = (1-ratio)*sec0.y + ratio*sec1.y
    sec.z  = (1-ratio)*sec0.z + ratio*sec1.z

    sec.cst_u, sec.cst_l = cst_foil_fit(sec.xx, sec.yu, sec.xx, sec.yl, n_order=sec0.cst_u.shape[0])

    return sec

def list_mul(list_, coef=1.0):
    '''
    Multiply each element in the list by coef
    '''
    if not isinstance(list_, list):
        print(str(list_))
        raise Exception('Can not use list_mul for a non-list object')
    
    temp = np.array(list_) * coef
    return temp.tolist()
