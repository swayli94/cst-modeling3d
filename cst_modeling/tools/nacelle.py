'''
Nacelle profiles for lofting to a surface of revolution.

- flow through nacelle (FTN);
- powered engine nacelle (PEN), or with powered nacelle (WPN);


References:

    [1] "Comparison of Overwing and Underwing Nacelle Aeropropulsion Optimization for Subsonic Transport Aircraft", Journal of Aircraft, 2024, Vol. 61, No. 2.
        https://doi.org/10.2514/1.C037508
        
    [2] "Non-axisymmetric aero-engine nacelle design by surrogate-based methods", Aerospace Science and Technology, 2021, Vol. 117, 106890.
        https://doi.org/10.1016/j.ast.2021.106890
        
    [3] "Impact of Droop and Scarf on the Aerodynamic Performance of Compact Aero-Engine Nacelles", AIAA SciTech, 2020.
        https://doi.org/10.2514/6.2020-1522

'''
import copy
import numpy as np

from typing import Dict, Tuple, List, Union

from ..math import rotate_vector, transform, cst_foil, dist_clustcos, interp_from_curve

import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve


class NacelleIntakeHighlight():
    '''
    The three-dimensional highlight point curve of the intake face.
    
    Parameters
    ----------
    offset_y_intake : float
        Offset of the intake face center from the axis of rotation.
        
    l_intake : float
        Offset of the intake face center from the fan face, i.e., the intake length.
        
    theta_droop : float
        Centre-line droop angle (degree), i.e., arctan(offset_y_intake/l_intake).
        After the revolution around the engine axis, the contoured intake was morphed by the droop angle, 
        measured as the offset between the engine and the inlet centre-line.
        
    theta_scarf : float
        Highlight scarf angle (degree).
        The geometry was scarfed by rotating the highlight plane around the intake face center.
        
    n_circum : int
        The number of circumferential points for the intake highlight curve.
        
    circum_control_psi : List[float]
        Circumferential location of control sections.
        
    circum_control_r_highlight : Union[float, List[float]]
        Radius of the intake highlight point to the intake face center.
        If float, the radius is constant for all sections.
    
    '''
    def __init__(self, l_intake: float, 
                    theta_droop: float, theta_scarf: float,
                    n_circum=101, 
                    circum_control_psi = [0.0, 90.0, 180.0, 270.0],
                    circum_control_r_highlight : Union[float, List[float]] = 1.0,
                ) -> None:
        
        if theta_droop + theta_scarf > 45.0:
            raise ValueError('The sum of theta_droop and theta_scarf should be less than 45.0 degree.')
        
        self.l_intake = l_intake
        self.theta_droop = theta_droop
        self.theta_scarf = theta_scarf
        self.n_circum = n_circum
        self.circum_control_psi = circum_control_psi + [360.0]
        
        if isinstance(circum_control_r_highlight, float):
            self.circum_control_r_highlight = [circum_control_r_highlight] * len(self.circum_control_psi)
        else:
            self.circum_control_r_highlight = circum_control_r_highlight + circum_control_r_highlight[0]
        
        self.intake_face_center = np.zeros(3)
        self.intake_face_center[0] = -l_intake
        self.intake_face_center[1] =  l_intake * np.tan(np.deg2rad(theta_droop))
        
    def calculate(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Calculate highlight point curve of the intake.
        
        Returns
        -------
        curve : np.ndarray [n_n_circum, 3]
            Highlight point coordinates in 3D space, i.e., (x, y, z).
            
        psi_curve: np.ndarray [n_n_circum, 4]
            Highlight point coordinates in 3D space, described in parametric coordinates,
            i.e., (psi, x, y, z).
        '''
        
        self.func_radius = CubicSpline(self.circum_control_psi, self.circum_control_r_highlight, bc_type='periodic')
        
        curve_psi = np.linspace(0.0, 360.0, self.n_circum, endpoint=True)
        
        curve_radius = self.func_radius(curve_psi)
        
        #* The intake highlight curve in the yz-plane before rotation
        
        curve = np.zeros((self.n_circum, 3))
        curve[:, 0] = self.intake_face_center[0] / np.cos(np.deg2rad(self.theta_droop))
        curve[:, 1] = curve_radius * np.cos(np.deg2rad(curve_psi))
        curve[:, 2] = curve_radius * np.sin(np.deg2rad(curve_psi))
        
        
        #* Rotate yz-plane about z-axis, the origin is the fan face center
        #  After the revolution around the engine axis, the contoured intake was morphed by the droop angle, 
        #  measured as the offset between the engine and the inlet centre-line.
        
        curve = rotate_vector(curve[:, 0], curve[:, 1], curve[:, 2], angle=self.theta_droop,
                       origin=[0.0, 0.0, 0.0], axis_vector=[0.0, 0.0, 1.0])
        
        #* Rotate xyz about z-axis, the origin is the intake face center
        #  The geometry was scarfed by rotating the highlight plane around the intake face center.
        
        curve = rotate_vector(curve[:, 0], curve[:, 1], curve[:, 2], angle=self.theta_scarf,
                       origin=self.intake_face_center, axis_vector=[0.0, 0.0, 1.0])
                
        avg = 0.5*(curve[0] + curve[-1])
        for i in range(avg.shape[0]):
            if abs(avg[i]) < 1e-6:
                avg[i] = 0.0
        
        curve[ 0,:] = avg
        curve[-1,:] = avg
           
        #* Project the highlight point to the profile plane
        
        new_curve_psi = [np.arctan2(z, y) for y, z in zip(curve[:self.n_circum, 1], curve[:self.n_circum, 2])]

        for i in range(len(new_curve_psi)):
            if new_curve_psi[i] < 0.0:
                new_curve_psi[i] += 2.0 * np.pi

        new_curve_psi[-1] = 2.0 * np.pi

        self.func_x = CubicSpline(new_curve_psi, curve[:,0], bc_type='periodic')
        self.func_y = CubicSpline(new_curve_psi, curve[:,1], bc_type='periodic')
        self.func_z = CubicSpline(new_curve_psi, curve[:,2], bc_type='periodic')
        
        psi_curve = np.zeros((self.n_circum, 4))
        psi_curve[:, 0] = curve_psi
        psi_curve[:, 1] = self.func_x(new_curve_psi)
        psi_curve[:, 2] = self.func_y(new_curve_psi)
        psi_curve[:, 3] = self.func_z(new_curve_psi)
                
        return curve, psi_curve

    def get_coordinate_3d(self, psi: float) -> np.ndarray:
        '''
        Get the highlight point coordinate at a given circumferential angle.
        
        Parameters
        ----------
        psi : float
            Circumferential angle (degree).
            
        Returns
        -------
        point : np.ndarray [3]
            Highlight point coordinate in 3D space, i.e., (x, y, z).
        '''
        psi = np.deg2rad(psi % 360.0)
        
        x = self.func_x(psi)
        y = self.func_y(psi)
        z = self.func_z(psi)
        
        return np.array([x, y, z])
    
    def get_coordinate_2d(self, psi: float) -> np.ndarray:
        '''
        Get the highlight point coordinate at a given circumferential angle.
        
        Parameters
        ----------
        psi : float
            Circumferential angle (degree).
            
        Returns
        -------
        point : np.ndarray [3]
            Highlight point coordinate in 3D space, i.e., (x, y, z).
        '''
        point = self.get_coordinate_3d(psi)
        
        y = np.sqrt(point[1]**2 + point[2]**2)
        
        return np.array([point[0], y])


class PoweredNacelleProfile():
    '''
    Two-dimensional profile of a powered engine nacelle.
    
    The reference point (0,0) is located at the intersection of rotation axis and the fan inlet face. 
    The definition of the profile parameters can be found in Figure 1 of [2].
    
    The intake face center is the same for all profiles in different meridional angles,
    which locates at (-l_intake, offset_y_intake, 0.0).
    
    
    Nacelle profile consists of the following curves:
    
    - Outer cowl surface: (3,4);
    - Inlet surface: (3,2);
    - Bypass duct outer and inner surfaces: (5,4), (6,7);
    - Core duct outer and inner surfaces: (9,8), (10,11);
    
    - Spinner (conical surface): (0,1);
    - Core cowl (conical surface): (7, 8);
    - Core plug (conical surface): (11,12);
    
    - Fan face intake (planar surface): (1,2);
    - Bypass duct inlet (planar surface): (5,6);
    - Core duct inlet (planar surface): (9,10);
    
    
    Nacelle profile has the following features:
    
    - l_fore: fore-body length;
    - l_lip: lip length;
    
    - r_if: initial fore-body radius;
    - r_throat: inlet throat radius;
    - r_max: maximum radius;

    - theta_airfoil: airfoil geometric incidence angle (degree);
    - theta_bp = bypass exit angle (degree);
    - theta_cr = core exit angle (degree);
    - theta_nac = nacelle boat tail angle (degree);
    
    - x_max: Axial location of maximum radius;
    - x_thr: Axial location of inlet throat;
    
    Parameters
    ----------
    psi : float
        Profile's meridional angle (degree), [0, 360].
        Zero angle corresponds to the half x-y plane with positive y.
    
    n_point_segment : int, optional
        Number of points for each segment of the profile. The default is 201.
    
    Attributes
    ----------
    params : Dict[str, float]
        Nacelle profile parameters.
        
    segment_points : Dict[int, np.ndarray]
        Nacelle profile segment points.
        
    profile_segments : Dict[str, np.ndarray]
        Nacelle profile segments.
    
    features : Dict[str, float]
        Nacelle profile features.
        
    profile_x, profile_y : np.ndarray [n_point_profile]
        Nacelle profile coordinates.
    
    '''

    def __init__(self, psi=0.0, n_point_segment=201) -> None:
        
        self.psi = psi
        self.n_point_segment = n_point_segment
        
        self.params : Dict[str, np.ndarray] = {}
        self.features : Dict[str, np.ndarray] = {}
        self.segment_points : Dict[int, np.ndarray] = {}
        self.profile_segments : Dict[int, np.ndarray] = {}
        
        self.profile_x : np.ndarray = None
        self.profile_y : np.ndarray = None
        
    def set_parameters(self, theta_spinner: float, r_spinner: float, r_fan: Union[None, float], 
                        highlight_x: float, highlight_y: float, intake_face_center: np.ndarray, 
                        l_nacelle: float, r_te: float,
                        l_fan: float, r_bypass_outer: float, r_bypass_inner: float,
                        x_core_cowl_0: float, y_core_cowl_0: float,
                        x_core_cowl_1: float, y_core_cowl_1: float,
                        x_core_duct: float, r_core_outer: float, r_core_inner: float,
                        x_core_plug_0: float, y_core_plug_0: float, x_core_plug_1: float,
                        cst_u: List[float], cst_l: List[float],
                        bypass_inner_angle: float = 0.0,
                        bypass_inner_control_points: List[Tuple[float, float]] = None,
                        core_outer_control_points: List[Tuple[float, float]] = None,
                        core_inner_control_points: List[Tuple[float, float]] = None,
                       ) -> None:
        '''
        Set profile parameters.
        
        Parameters
        ----------
        cst_u, cst_l : List[float]
            Upper and lower CST coefficients of the airfoil profile.
        
        bypass_inner_angle : float
            Bypass duct inner surface half-angle (degree) at point (6).
        
        bypass_inner_control_points : List[Tuple[float, float]], optional
            Control points for the bypass duct inner surface. The default is None.
            
        core_outer_control_points : List[Tuple[float, float]], optional
            Control points for the core duct outer surface. The default is None.
            
        core_inner_control_points : List[Tuple[float, float]], optional
            Control points for the core duct inner surface. The default is None.
        
        Profile parameters:
        
        - theta_spinner: half-angle of the spinner cone (degree) (0);
        
        - r_spinner: spinner radius (1);
        - r_fan: fan radius (optional) (2);
        
        - intake_face_center: intake face center, ndarray [3];
        - highlight_x: highlight point x-coordinate (3);
        - highlight_y: highlight point y-coordinate (3);
        
        - l_nacelle: nacelle length in the x-axis direction (4);
        - r_te: nacelle trailing edge radius (4);
        
        - l_fan: fan length in the x-axis direction (5, 6);
        - r_bypass_outer: bypass duct outer radius (5);
        - r_bypass_inner: bypass duct inner radius (6);

        - x_core_cowl_0: Axial location of core cowl start (7);
        - y_core_cowl_0: Radial location of core cowl start (7);
        - x_core_cowl_1: Axial location of core cowl end (8);
        - y_core_cowl_1: Radial location of core cowl end (8);
        
        - x_core_duct: Axial location of core duct inlet (9);
        - r_core_outer: core duct outer radius (9);
        - r_core_inner: core duct inner radius (10);
        
        - x_core_plug_0: Axial location of core plug start (11);
        - y_core_plug_0: Radial location of core plug start (11);
        - x_core_plug_1: Axial location of core plug end (12);
        
        '''
        self.params['r_spinner'] = r_spinner
        self.params['theta_spinner'] = theta_spinner
        self.params['r_fan'] = r_fan
        self.params['highlight_x'] = highlight_x
        self.params['highlight_y'] = highlight_y
        self.params['intake_face_center'] = intake_face_center
        self.params['l_nacelle'] = l_nacelle
        self.params['r_te'] = r_te
        self.params['l_fan'] = l_fan
        self.params['r_bypass_outer'] = r_bypass_outer
        self.params['r_bypass_inner'] = r_bypass_inner
        self.params['x_core_cowl_0'] = x_core_cowl_0
        self.params['y_core_cowl_0'] = y_core_cowl_0
        self.params['x_core_cowl_1'] = x_core_cowl_1
        self.params['y_core_cowl_1'] = y_core_cowl_1
        self.params['x_core_duct'] = x_core_duct
        self.params['r_core_outer'] = r_core_outer
        self.params['r_core_inner'] = r_core_inner
        self.params['x_core_plug_0'] = x_core_plug_0
        self.params['y_core_plug_0'] = y_core_plug_0
        self.params['x_core_plug_1'] = x_core_plug_1
        
        self.params['cst_u'] = cst_u
        self.params['cst_l'] = cst_l
        
        self.params['bypass_inner_angle'] = bypass_inner_angle
        self.params['bypass_inner_control_points'] = bypass_inner_control_points
        self.params['core_outer_control_points'] = core_outer_control_points
        self.params['core_inner_control_points'] = core_inner_control_points
    
    def get_profile(self) -> Dict[np.ndarray, np.ndarray]:
        '''
        Get nacelle profile segments.
        
        Returns
        -------
        profile_x, profile_y : np.ndarray 
            Nacelle profile coordinates.
        '''
        self.calculate_segment_points()
        self.calculate_profile_segments()
        self.calculate_profile_features()
        
        profile_x = [self.profile_segments[0][:, 0]] + [self.profile_segments[i+1][1:, 0] for i in range(self.n_segment-1)]
        profile_y = [self.profile_segments[0][:, 1]] + [self.profile_segments[i+1][1:, 1] for i in range(self.n_segment-1)]
        
        self.profile_x = np.concatenate(profile_x)
        self.profile_y = np.concatenate(profile_y)

        return self.profile_x, self.profile_y
    
    @property
    def n_point_profile(self) -> int:
        '''
        Number of points for the nacelle profile.
        '''
        return self.profile_x.shape[0]
    
    @property
    def n_segment(self) -> int:
        '''
        Number of nacelle profile segments.
        '''
        return len(self.profile_segments)
    
    def calculate_segment_points(self) -> None:
        '''
        Calculate nacelle profile segment end points based on the profile parameters.
        '''
        self.segment_points[0] = np.array([-self.params['r_spinner'] / np.tan(np.deg2rad(self.params['theta_spinner'])), 0.0])
        
        self.segment_points[1] = np.array([0.0, self.params['r_spinner']])
        
        self.segment_points[2] = np.array([0.0, self.params['r_fan']])
        
        self.segment_points[3] = np.array([self.params['highlight_x'], self.params['highlight_y']])
        
        self.segment_points[4] = np.array([self.params['intake_face_center'][0] + self.params['l_nacelle'], self.params['r_te']])

        self.segment_points[5] = np.array([self.params['l_fan'], self.params['r_bypass_outer']])
        
        self.segment_points[6] = np.array([self.params['l_fan'], self.params['r_bypass_inner']])
        
        self.segment_points[7] = np.array([self.params['x_core_cowl_0'], self.params['y_core_cowl_0']])
        
        self.segment_points[8] = np.array([self.params['x_core_cowl_1'], self.params['y_core_cowl_1']])
        
        self.segment_points[9] = np.array([self.params['x_core_duct'], self.params['r_core_outer']])
        
        self.segment_points[10] = np.array([self.params['x_core_duct'], self.params['r_core_inner']])
        
        self.segment_points[11] = np.array([self.params['x_core_plug_0'], self.params['y_core_plug_0']])
        
        self.segment_points[12] = np.array([self.params['x_core_plug_1'], 0.0])

    def calculate_profile_segments(self) -> None:
        '''
        Calculate nacelle profile segments:
        
        - Spinner (conical surface): (0,1);
        - Fan face intake (planar surface): (1,2);
        - Inlet surface: (3,2);
        - Outer cowl surface: (3,4);
        - Bypass duct outer surface: (5,4);
        - Bypass duct inlet (planar surface): (5,6);
        - Bypass duct inner surface: (6,7);
        - Core cowl (conical surface): (7, 8);
        - Core duct outer surface: (9,8);
        - Core duct inlet (planar surface): (9,10);
        - Core duct inner surface: (10,11);
        - Core plug (conical surface): (11,12);
        '''
        
        #* (0,1) Spinner (conical surface)
        self.profile_segments[0] = self.create_straight_line(   self.segment_points[0][0], self.segment_points[0][1],
                                                                self.segment_points[1][0], self.segment_points[1][1])
                
        #* (3,4) Outer cowl surface
        self.outer_cowl_profile()

        #* (2,3) Inlet surface
        self.inlet_profile()

        #* (4,5) Bypass duct outer surface
        self.bypass_duct_outer_profile()
        
        #* (1,2) Fan face intake (planar surface)
        self.profile_segments[1] = self.create_straight_line(   self.segment_points[1][0], self.segment_points[1][1],
                                                                self.segment_points[2][0], self.segment_points[2][1])
        
        #* (5,6) Bypass duct inlet (planar surface)
        self.profile_segments[5] = self.create_straight_line(   self.segment_points[5][0], self.segment_points[5][1],
                                                                self.segment_points[6][0], self.segment_points[6][1])
        
        #* (6,7) Bypass duct inner surface
        self.bypass_duct_inner_profile()
        
        #* (7, 8) Core cowl (conical surface)
        self.profile_segments[7] = self.create_straight_line(   self.segment_points[7][0], self.segment_points[7][1],
                                                                self.segment_points[8][0], self.segment_points[8][1])
        
        #* (8, 9) Core duct outer surface
        self.core_duct_outer_profile()
        
        #* (9,10) Core duct inlet (planar surface)
        self.profile_segments[9] = self.create_straight_line(   self.segment_points[9][0], self.segment_points[9][1],
                                                                self.segment_points[10][0], self.segment_points[10][1])
        
        #* (10,11) Core duct inner surface
        self.core_duct_inner_profile()
        
        #* (11,12) Core plug (conical surface)
        self.profile_segments[11] = self.create_straight_line(  self.segment_points[11][0], self.segment_points[11][1],
                                                                self.segment_points[12][0], self.segment_points[12][1])

    def calculate_profile_features(self) -> None:
        '''
        Calculate nacelle profile features:
        
        - l_fore: fore-body length;
        - l_lip: lip length;
        
        - r_if: initial fore-body radius;
        - r_throat: inlet throat radius;
        - r_max: maximum radius;

        - theta_airfoil: airfoil geometric incidence angle (degree);
        - theta_bp = bypass exit angle (degree);
        - theta_cr = core exit angle (degree);
        - theta_nac = nacelle boat tail angle (degree);
        
        - x_max: Axial location of maximum radius;
        - x_thr: Axial location of inlet throat;
        '''
        pass

    def create_straight_line(self, x0: float, y0: float, x1: float, y1: float) -> np.ndarray:
        '''
        Generate a straight line
        
        Parameters
        ----------
        x0, y0 : float
            Starting point coordinates.
            
        x1, y1 : float
            Ending point coordinates.
            
        Returns
        -------
        line : np.ndarray [n_point_segment, 2]
            Straight line coordinates.
        '''
        line = np.zeros((self.n_point_segment, 2))
        line[:, 0] = np.linspace(x0, x1, self.n_point_segment, endpoint=True)
        
        if abs(x0 - x1) < 1e-6:
            line[:, 1] = np.linspace(y0, y1, self.n_point_segment, endpoint=True)
        else:
            line[:, 1] = y0 + (y1 - y0) / (x1 - x0) * (line[:, 0] - x0)
        
        return line

    def outer_cowl_profile(self) -> None:
        '''
        Generate outer cowl curve, i.e., profile segment (3) between segment points (3,4).
        '''
        
        #* Airfoil profile
        self._xx, self._yu, self._yl, _, _ = cst_foil(
                self.n_point_segment, self.params['cst_u'], self.params['cst_l'])
        
        scale = np.linalg.norm(self.segment_points[4] - self.segment_points[3])
        angle = np.arctan2( self.segment_points[4][1] - self.segment_points[3][1],
                            self.segment_points[4][0] - self.segment_points[3][0])
        angle = np.rad2deg(angle)
        
        self._foil_xu, self._foil_xl, self._foil_yu, self._foil_yl = transform(
                    self._xx, self._xx, self._yu, self._yl, scale=scale, rot=angle, 
                    dx=self.segment_points[3][0], dy=self.segment_points[3][1])
        
        self.features['theta_airfoil'] = angle
        self.features['scale_airfoil'] = scale
        
        #* Outer cowl curve
        self.profile_segments[3] = np.concatenate(
            (self._foil_xu[:, None], self._foil_yu[:, None]), axis=1)

    def inlet_profile(self) -> None:
        '''
        Generate inlet curve, i.e., profile segment (2) between segment points (2,3).
        
        Must run after `outer_cowl_profile`.
        '''
        xx = dist_clustcos(self.n_point_segment, a0=0.5, a1=0.99, beta=1.0)
        
        curve = np.zeros((self.n_point_segment, 2))
        
        curve[:, 0] = xx * (self.segment_points[3][0] - self.segment_points[2][0]) + self.segment_points[2][0]
        
        scale = 1.0
        
        #* Solve scale to match the fan radius
        
        if self.params['r_fan'] is not None:

            def func(scale):

                _, _foil_xl, _, _foil_yl = transform(
                            self._xx, self._xx, scale*self._yl, scale*self._yl, 
                            scale=self.features['scale_airfoil'], 
                            rot=self.features['theta_airfoil'], 
                            dx=self.segment_points[3][0], dy=self.segment_points[3][1])
                
                y = interp_from_curve(self.segment_points[2][0], _foil_xl, _foil_yl, extrapolate=True)
                
                return abs(y-self.params['r_fan'])
                
            root = fsolve(func, x0=1.0)
            scale = root[0]
            
            print('> Inlet surface scale: %.2f'%(scale))
        
        _, _foil_xl, _, _foil_yl = transform(
                    self._xx, self._xx, scale*self._yl, scale*self._yl, 
                    scale=self.features['scale_airfoil'], 
                    rot=self.features['theta_airfoil'], 
                    dx=self.segment_points[3][0], dy=self.segment_points[3][1])
        
        curve[:, 1] = interp_from_curve(curve[:, 0], _foil_xl, _foil_yl, extrapolate=True)
        
        self.profile_segments[2] = curve
        
        #* Update segment point and parameters
        self.segment_points[2][1] = curve[0, 1]
        self.params['r_fan'] = curve[0, 1]

    def bypass_duct_outer_profile(self) -> None:
        '''
        Generate bypass duct outer surface curve, i.e., profile segment (4) between segment points (4,5).
        
        Must run after `outer_cowl_profile`.
        '''
        xx = dist_clustcos(self.n_point_segment, a0=0.1, a1=0.9, beta=1.0)
        
        curve = np.zeros((self.n_point_segment, 2))
        
        curve[:, 0] = xx * (self.segment_points[5][0] - self.segment_points[4][0]) + self.segment_points[4][0]
        
        scale = 1.0
        
        #* Solve scale to match the fan radius
        
        if self.params['r_bypass_outer'] is not None:

            def func(scale):
                
                _, _foil_xl, _, _foil_yl = transform(
                            self._xx, self._xx, scale*self._yl, scale*self._yl, 
                            scale=self.features['scale_airfoil'], 
                            rot=self.features['theta_airfoil'], 
                            dx=self.segment_points[3][0], dy=self.segment_points[3][1])
                
                y = interp_from_curve(self.segment_points[5][0], _foil_xl, _foil_yl, extrapolate=True)
                
                return abs(y-self.params['r_bypass_outer'])
                
            root = fsolve(func, x0=1.0)
            scale = root[0]
            
            print('> Bypass duct outer surface scale: %.2f'%(scale))
        
        _, _foil_xl, _, _foil_yl = transform(
                    self._xx, self._xx, scale*self._yl, scale*self._yl, 
                    scale=self.features['scale_airfoil'], 
                    rot=self.features['theta_airfoil'], 
                    dx=self.segment_points[3][0], dy=self.segment_points[3][1])
        
        curve[:, 1] = interp_from_curve(curve[:, 0], _foil_xl, _foil_yl, extrapolate=True)
        
        self.profile_segments[4] = curve
        
        #* Update segment point and parameters
        self.segment_points[5][1] = curve[-1, 1]
        self.params['r_bypass_outer'] = curve[0, 1]

    def bypass_duct_inner_profile(self) -> None:
        '''
        Generate bypass duct inner surface curve, i.e., profile segment (6) between segment points (6,7).
        '''       
        
        if isinstance(self.params['bypass_inner_angle'], float):
            bcx0 = (1, np.tan(np.deg2rad(self.params['bypass_inner_angle'])))
        else:
            bcx0 = (2, 0.0)
            
        bcx1 = (1, (self.segment_points[8][1] - self.segment_points[7][1]) / (self.segment_points[8][0] - self.segment_points[7][0]))
        
        bc_type = (bcx0, bcx1)
        
        control_points = copy.deepcopy(self.params['bypass_inner_control_points'])
        
        xs = [self.segment_points[6][0]] + [x for x, _ in control_points] + [self.segment_points[7][0]]
        ys = [self.segment_points[6][1]] + [y for _, y in control_points] + [self.segment_points[7][1]]
        
        func = CubicSpline(xs, ys, bc_type=bc_type)

        curve = np.zeros((self.n_point_segment, 2))
        curve[:, 0] = np.linspace(self.segment_points[6][0], self.segment_points[7][0], self.n_point_segment, endpoint=True)
        curve[:, 1] = func(curve[:, 0])
        
        self.profile_segments[6] = curve

    def core_duct_outer_profile(self) -> None:
        '''
        Generate core duct outer surface curve, i.e., profile segment (8) between segment points (8,9).
        '''       
        
        bcx0 = (1, 0.0)
        bcx1 = (1, (self.segment_points[12][1] - self.segment_points[11][1]) / (self.segment_points[12][0] - self.segment_points[11][0]))
        
        bc_type = (bcx0, bcx1)
        
        control_points = copy.deepcopy(self.params['core_outer_control_points'])
        
        xs = [self.segment_points[9][0]] + [x for x, _ in control_points] + [self.segment_points[8][0]]
        ys = [self.segment_points[9][1]] + [y for _, y in control_points] + [self.segment_points[8][1]]
        
        func = CubicSpline(xs, ys, bc_type=bc_type)

        curve = np.zeros((self.n_point_segment, 2))
        curve[:, 0] = np.linspace(self.segment_points[8][0], self.segment_points[9][0], self.n_point_segment, endpoint=True)
        curve[:, 1] = func(curve[:, 0])
        
        self.profile_segments[8] = curve

    def core_duct_inner_profile(self) -> None:
        '''
        Generate core duct inner surface curve, i.e., profile segment (10) between segment points (10,11).
        '''       
        
        bcx0 = (1, 0.0)
        bcx1 = (1, (self.segment_points[12][1] - self.segment_points[11][1]) / (self.segment_points[12][0] - self.segment_points[11][0]))
        
        bc_type = (bcx0, bcx1)
        
        control_points = copy.deepcopy(self.params['core_outer_control_points'])
        
        xs = [self.segment_points[10][0]] + [x for x, _ in control_points] + [self.segment_points[11][0]]
        ys = [self.segment_points[10][1]] + [y for _, y in control_points] + [self.segment_points[11][1]]
        
        func = CubicSpline(xs, ys, bc_type=bc_type)

        curve = np.zeros((self.n_point_segment, 2))
        curve[:, 0] = np.linspace(self.segment_points[10][0], self.segment_points[11][0], self.n_point_segment, endpoint=True)
        curve[:, 1] = func(curve[:, 0])
        
        self.profile_segments[10] = curve


    def plot(self, show=True) -> None:
        '''
        Plot nacelle profile.
        '''
        fig, ax = plt.subplots(figsize=(16, 8))
        
        plt.plot(self.profile_x, self.profile_y, 'k-', lw=2, label='Nacelle profile')
        
        plt.xlim((-0.2, 1.2))
        plt.ylim((-0.2, 0.2))
        plt.axis('equal')
        plt.legend()
        
        def add_label(x, y, label, dx=0.01, dy=0.01, color='b'):
            ax.plot(x, y, color+'*')
            ax.text(x+dx, y+dy, label, color=color)
        
        for i_point in range(len(self.segment_points)):
            
            add_label(self.segment_points[i_point][0], self.segment_points[i_point][1], str(i_point))
        
        if show:
            plt.show()
