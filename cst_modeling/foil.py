'''
Airfoil geometric features and modification functions
'''
from typing import Tuple, List

import numpy as np
from scipy.interpolate import interp1d

from .math import curve_curvature, find_circle_3p, CoordinateTransformation
from .section import cst_foil, cst_foil_fit, bump_function, foil_increment


def check_validity(x: np.ndarray, yu: np.ndarray, yl: np.ndarray, eps=1e-6) -> bool:
    '''
    Check if the airfoil is valid
    
    Returns
    -------
    is_valid : bool
        True if the airfoil is valid, False otherwise
    '''
    if len(x) != len(yu) or len(x) != len(yl):
        print('>>> Error [check_validity]')
        print('    Airfoil x, yu, and yl have different lengths.')
        return False
    
    if abs(x[0]) > eps or abs(x[-1]-1.0) > eps:
        print('>>> Error [check_validity]')
        print('    Airfoil x coordinates are not normalized.')
        print('    x[0] = ', x[0], ', x[-1] = ', x[-1])
        return False
    else:
        x[0] = 0.0
        x[-1] = 1.0
    
    if abs(yu[-1] + yl[-1]) > eps:
        print('>>> Error [check_validity]')
        print('    Airfoil trailing edge is not symmetrical about the x-axis.')
        return False        
    
    if abs(yu[0] - yl[0]) > eps:
        print('>>> Error [check_validity]')
        print('    Airfoil leading edge is not connected.')
        return False
    
    if not np.all(yu >= yl):
        print('>>> Error [check_validity]')
        print('    Airfoil upper surface is below lower surface.')
        return False
    
    return True


class FoilGeoFeatures():
    '''
    Airfoil geometric features (unit chord length and sharp trailing edge)
    
    Parameters
    ----------
    x : np.ndarray
        Airfoil x-coordinates
        
    yu, yl : np.ndarray
        Airfoil surface y-coordinates
    '''
    def __init__(self, x: np.ndarray, yu: np.ndarray, yl: np.ndarray) -> None:
        
        self.x = x
        self.yu = yu
        self.yl = yl
        self.nn = x.shape[0]
        
        if not check_validity(x, yu, yl):
            raise ValueError('Invalid airfoil geometry')
        
        self.thickness = self.yu - self.yl
        self.camber = 0.5*(self.yu + self.yl)
        
        self.tail = self.thickness[-1]
        
        self.features ={'x': self.x, 'yu': self.yu, 'yl': self.yl, 
                        'thickness': self.thickness, 'camber': self.camber,
                        'tail': self.tail}
        
        self.interp_func_yu = None
        self.interp_func_yl = None
        
        #* Constants
        self.X_RLE = 0.005
        self.X_LE = 0.02
        self.X_TE = 0.98
            
    @property
    def feature_names(self) -> Tuple[str]:
        return ('x', 'yu', 'yl', 'thickness', 'camber', 'curvature_upper', 'curvature_lower',
                't_max', 'x_t', 'i_t', 'c_max', 'x_c', 'i_c', 
                'x_uc', 'y_uc', 'i_uc', 'x_lc', 'y_lc', 'i_lc',
                'volume', 'c_mean', 'c_weight_mean', 'x_c_weight_mean',
                't_20', 't_70', 'c_f60', 'c_r40',
                'r_le', 'wedge_angle', 'slope_angle_le', 'slope_angle_te'
                )
    
    def get_feature(self, feature_name: str):
        
        if feature_name not in self.feature_names:
            
            print('>>> Error [FoilFeatures.get_feature()]')
            print('    Available features: ', self.feature_names)
            print('    Feature not available: ', feature_name)
            raise ValueError()
        
        elif feature_name not in self.features.keys():
            
            print('>>> Error [FoilFeatures.get_feature()]')
            print('    Feature not computed yet: ', feature_name)
            raise ValueError()
        
        return self.features[feature_name]

    def interp_y(self, x0: float, side='upper') -> float:
        '''
        Interpolate value from the airfoil surfaces
        
        Parameters
        ----------
        x0 : float or ndarray
            ndarray/value of x locations to be interpolated.

        Returns
        ----------
        y0 : float or ndarray
            interpolated coordinates
        '''
        if self.interp_func_yu is None:
            self.interp_func_yu = interp1d(self.x, self.yu, kind='cubic')
            
        if self.interp_func_yl is None:
            self.interp_func_yl = interp1d(self.x, self.yl, kind='cubic')
        
        if side == 'upper':
            y0 = self.interp_func_yu(x0)

        elif side == 'lower':
            y0 = self.interp_func_yl(x0)
            
        else:
            print('>>> Error [FoilFeatures.interp_y()]')
            print('    Invalid side: ', side)
            raise ValueError('    Side must be either "upper" or "lower"')

        return y0

    def get_maximum_thickness(self) -> Tuple[float, float, int]:
        '''
        Get airfoil maximum thickness and its location
        
        Returns
        -------
        t_max : float
            Maximum thickness
            
        x_t : float
            Maximum thickness location
            
        i_t : int
            Maximum thickness location index
        '''
        i_t = np.argmax(self.thickness)
        x_t = self.x[i_t]
        
        self.features['t_max'] = self.thickness[i_t]
        self.features['x_t'] = x_t
        self.features['i_t'] = i_t
        
        return self.thickness[i_t], x_t, i_t

    def get_thickness_at(self, x0: float) -> float:
        '''
        Get airfoil thickness at a specific location
        
        Parameters
        ----------
        x0 : float
            X-coordinate location
        
        Returns
        -------
        t0 : float
            Airfoil thickness at x0
        '''
        yu0 = self.interp_y(x0, side='upper')
        yl0 = self.interp_y(x0, side='lower')
        
        t0 = yu0 - yl0
        
        name = 't_%d'%int(x0*100)
        
        if name in self.feature_names:
            self.features[name] = t0
        
        return t0

    def get_maximum_camber(self) -> Tuple[float, float, int]:
        '''
        Get airfoil maximum camber
        
        Returns
        -------
        c_max : float
            Maximum camber
            
        x_c : float
            Maximum camber location
            
        i_c : int
            Maximum camber location index
        '''
        i_c = np.argmax(self.camber)
        x_c = self.x[i_c]
        
        self.features['c_max'] = self.camber[i_c]
        self.features['x_c'] = x_c
        self.features['i_c'] = i_c
        
        return self.camber[i_c], x_c, i_c

    def get_volume(self) -> float:
        '''
        Get airfoil volume.
        Volume is also a measure of the airfoil thickness.
        
        Returns
        -------
        volume : float
            Airfoil volume
        '''
        volume = np.trapz(self.thickness, self.x)
        
        self.features['volume'] = volume
        
        return volume

    def get_average_camber(self) -> float:
        '''
        Get airfoil average camber, i.e., the area of camber line.
        
        c_mean = (the area of camber line) / (the chord length)
        
        Returns
        -------
        c_mean : float
            Average camber
        '''
        c_mean = np.trapz(self.camber, self.x)
        
        self.features['c_mean'] = c_mean
        
        return c_mean

    def get_weighted_average_camber(self) -> Tuple[float, float]:
        '''
        Get airfoil weighted average camber, i.e.,
        the camber is averaged by using thickness as weights.
        
        Returns
        -------
        c_weight_mean : float
            Weighted average camber
            
        x_c_weight_mean : float
            Weighted average camber location
        '''
        c_weight_mean = 0.0
        x_c_weight_mean = 0.0
        
        sum_cambers = 0.0
        
        for i in range(self.nn-1):
            
            mc = (self.camber[i] + self.camber[i+1])*0.5
            mt = (self.thickness[i] + self.thickness[i+1])*0.5
            dx = self.x[i+1] - self.x[i]
            
            c_weight_mean += mc * mt * dx
            x_c_weight_mean += mc * dx * 0.5 * (self.x[i] + self.x[i+1])
            
            sum_cambers += mc * dx
            
        if 'volume' not in self.features.keys():
            self.get_volume()
            
        c_weight_mean = c_weight_mean / self.features['volume']
        x_c_weight_mean = x_c_weight_mean / sum_cambers
        
        self.features['c_weight_mean'] = c_weight_mean
        self.features['x_c_weight_mean'] = x_c_weight_mean
        
        return c_weight_mean, x_c_weight_mean

    def get_average_camber_front_60p(self) -> float:
        '''
        Get airfoil average camber in the front 60% of the chord length,
        i.e., the area of camber line in the front 60% of the chord length.
        
        Returns
        -------
        c_f60 : float
            Average camber in the front 60% of the chord length
        '''
        c_f60 = 0.0

        for i in range(self.nn-1):
                
            mc = (self.camber[i] + self.camber[i+1])*0.5
            dx = self.x[i+1] - self.x[i]
            
            if self.x[i] < 0.6:
                c_f60 += mc * dx
        
        self.features['c_f60'] = c_f60
        
        return c_f60
    
    def get_average_camber_rear_40p(self) -> float:
        '''
        Get airfoil average camber in the rear 40% of the chord length,
        i.e., the area of camber line in the rear 40% of the chord length.
        
        Returns
        -------
        c_r40 : float
            Average camber in the rear 40% of the chord length
        '''
        c_r40 = 0.0

        for i in range(self.nn-1):
                
            mc = (self.camber[i] + self.camber[i+1])*0.5
            dx = self.x[i+1] - self.x[i]
            
            if self.x[i] > 0.6:
                c_r40 += mc * dx
        
        self.features['c_r40'] = c_r40
        
        return c_r40

    def get_curvature(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Get airfoil curvature.
        Curvature is the second derivative of the camber line.
        
        Returns
        -------
        curvature_upper, curvature_lower : np.ndarray
            Curvature of airfoil upper surface and lower surface
        '''
        curvature_upper = curve_curvature(self.x, self.yu)
        curvature_lower = curve_curvature(self.x, self.yl)
        
        self.features['curvature_upper'] = curvature_upper
        self.features['curvature_lower'] = curvature_lower
        
        return curvature_upper, curvature_lower

    def get_leading_edge_radius(self) -> float:
        '''
        Get airfoil leading edge radius.
        
        Returns
        -------
        r_le : float
            Leading edge radius
        '''
        yu_RLE = self.interp_y(self.X_RLE, side='upper')
        yl_RLE = self.interp_y(self.X_RLE, side='lower')
        r_le, _ = find_circle_3p([0.0, 0.0], [self.X_RLE, yu_RLE], [self.X_RLE, yl_RLE])
        
        self.features['r_le'] = r_le
        
        return r_le

    def get_leading_edge_slope_angle(self) -> float:
        '''
        Get airfoil leading edge slope angle.
        
        Returns
        -------
        slope_angle : float
            Leading edge slope angle (degree)
        '''
        yu_LE = self.interp_y(self.X_LE, side='upper')
        yl_LE = self.interp_y(self.X_LE, side='lower')
        
        slope_angle = np.arctan(0.5*(yu_LE+yl_LE)/self.X_LE)
        slope_angle = np.rad2deg(slope_angle)
        
        self.features['slope_angle_le'] = slope_angle
        
        return slope_angle

    def get_trailing_edge_wedge_angle(self) -> float:
        '''
        Get airfoil trailing edge wedge angle.
        
        Returns
        -------
        wedge_angle : float
            Trailing edge wedge angle (degree)
        '''
        yu_TE = self.interp_y(self.X_TE, side='upper') - 0.5*self.tail
        yl_TE = self.interp_y(self.X_TE, side='lower') + 0.5*self.tail
        
        wedge_angle = np.arctan(yu_TE/(1.0-self.X_TE)) - np.arctan(yl_TE/(1.0-self.X_TE))
        wedge_angle = np.rad2deg(wedge_angle)
        
        self.features['wedge_angle'] = wedge_angle
        
        return wedge_angle

    def get_trailing_edge_slope_angle(self) -> float:
        '''
        Get airfoil trailing edge slope angle.
        
        Returns
        -------
        slope_angle : float
            Trailing edge slope angle (degree)
        '''
        yu_TE = self.interp_y(self.X_TE, side='upper')
        yl_TE = self.interp_y(self.X_TE, side='lower')
        
        slope_angle = np.arctan(-0.5*(yu_TE+yl_TE)/(1.0-self.X_TE))
        slope_angle = np.rad2deg(slope_angle)
        
        self.features['slope_angle_te'] = slope_angle
        
        return slope_angle

    def get_upper_crest_point(self) -> Tuple[float, float, int]:
        '''
        Get airfoil upper crest position and curvature.
        
        Returns
        -------
        x_uc : float
            X-coordinate of the upper crest point
            
        y_uc : float
            Y-coordinate of the upper crest point
            
        i_uc : int
            Index of the upper crest point
        '''
        i_uc = np.argmax(self.yu)
        x_uc = self.x[i_uc]
        y_uc = self.yu[i_uc]
        
        self.features['x_uc'] = x_uc
        self.features['y_uc'] = y_uc
        self.features['i_uc'] = i_uc

        return x_uc, y_uc, i_uc
    
    def get_lower_crest_point(self) -> Tuple[float, float, int]:
        '''
        Get airfoil lower crest position and curvature.
        
        Returns
        -------
        x_lc : float
            X-coordinate of the lower crest point
            
        y_lc : float
            Y-coordinate of the lower crest point
            
        i_lc : int
            Index of the lower crest point
        '''
        i_lc = np.argmin(self.yl)
        x_lc = self.x[i_lc]
        y_lc = self.yl[i_lc]
        
        self.features['x_lc'] = x_lc
        self.features['y_lc'] = y_lc
        self.features['i_lc'] = i_lc

        return x_lc, y_lc, i_lc


class FoilModification():
    '''
    Modify an airfoil by specifying geometric features and adding bumps
    
    Parameters
    ----------
    x : np.ndarray
        Airfoil x-coordinates
        
    yu, yl : np.ndarray
        Airfoil surface y-coordinates
        
    n_cst : int
        Number of CST coefficients for the upper and lower surfaces
    '''
    def __init__(self, x: np.ndarray, yu: np.ndarray, yl: np.ndarray, n_cst=10) -> None:
        
        if not check_validity(x, yu, yl):
            raise ValueError('Invalid airfoil geometry')
        
        self.x = x.copy()
        self.yu = yu.copy()
        self.yl = yl.copy()
        
        self.n_cst = n_cst
        
        self.nn = x.shape[0]
        
        #* Constants
        self.X_LE = 0.05
        self.X_TE = 0.95
        
        self.MAX_TRY = 5
        
    @property
    def tail(self) -> float:
        '''
        Airfoil trailing edge thickness
        '''
        return self.yu[-1] - self.yl[-1]
        
    def get_cst_coefficients(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Get CST coefficients for the current airfoil
        '''
        cst_u, cst_l = cst_foil_fit(self.x, self.yu, self.yl, n_cst=self.n_cst)
        
        return cst_u, cst_l
    
    def set_thickness(self, t_new: float) -> None:
        '''
        Change airfoil thickness to a new value
        '''
        t_old = np.max(self.yu - self.yl)

        self.yu = self.yu * t_new / t_old
        self.yl = self.yl * t_new / t_old

    def set_thickness_at(self, x0: float, t_new: float, width_bump=0.4
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change the airfoil thickness at location `x0` to a new value `t_new`.
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        t_old = geo.get_thickness_at(x0)
        
        n_try = 0
        
        while abs(t_new-t_old) > 1E-4 and n_try < self.MAX_TRY:
        
            cst_u, cst_l, _, _ = self.add_bump_to_thickness(
                x0, t_new-t_old, width_bump, kind='H', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            t_old = geo.get_thickness_at(x0)

            n_try += 1
            
        return cst_u, cst_l, t_old, t_new

    def set_maximum_thickness_location(self, x_t_new: float, slope0=1.0, slope1=1.0) -> None:
        '''
        Change airfoil maximum thickness location
        '''
        _, x_t, _ = FoilGeoFeatures(self.x, self.yu, self.yl).get_maximum_thickness()
        
        xShift = CoordinateTransformation()
        xShift.set_function_by_interpolation(x=[x_t], xp=[x_t_new], slope0=slope0, slope1=slope1)
        
        xx = self.x.copy()
        
        self.x = xShift.transform(self.x)
        
        # import matplotlib.pyplot as plt
        # plt.plot(xx, self.x, 'k')
        # plt.xlabel('x')
        # plt.ylabel('x\'')
        # plt.axis('equal')
        # plt.show()

    def set_camber_front(self, c_new: float, width_bump=0.9
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change the airfoil average camber of the front 60% part.
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        cf_old = geo.get_average_camber_front_60p()
        
        x0 = 0.3
        
        n_try = 0
        
        while abs(c_new-cf_old) > 1E-4 and n_try < self.MAX_TRY:
        
            cst_u, cst_l, _, _ = self.add_bump_to_camber(
                x0, 2*(c_new-cf_old), width_bump, kind='H', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            cf_old = geo.get_average_camber_front_60p()

            n_try += 1
            
        return cst_u, cst_l, cf_old, c_new

    def set_camber_rear(self, c_new: float, width_bump=0.6
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change the airfoil average camber of the rear 40% part.
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        cf_old = geo.get_average_camber_rear_40p()
        
        x0 = 0.8
        
        n_try = 0
        
        while abs(c_new-cf_old) > 1E-4 and n_try < self.MAX_TRY:
        
            cst_u, cst_l, _, _ = self.add_bump_to_camber(
                x0, 2*(c_new-cf_old), width_bump, kind='H', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            cf_old = geo.get_average_camber_rear_40p()

            n_try += 1
            
        return cst_u, cst_l, cf_old, c_new

    def add_bump(self, bumps: List[Tuple[float, float, float, str, str]],
                    keep_tmax=False) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Add bumps to the airfoil
        
        Parameters
        ----------
        bumps : list of tuple
            Bump parameters: (xc, h, w, side, kind)
            
        keep_tmax : bool
            Keep the maximum thickness of the airfoil
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        rLE : float
            Leading edge radius
        '''
        if keep_tmax:
            tmax = np.max(self.yu - self.yl)
        else:
            tmax = None

        for xc, h, w, side, kind in bumps:
            
            if side == 'upper':

                self.yu = self.yu + bump_function(self.x, xc, h, w, kind)
                
            elif side == 'lower':
                
                self.yl = self.yl + bump_function(self.x, xc, h, w, kind)
                
            else:
                
                print('>>> Error [FoilModification.add_bump()]')
                print('    Invalid side: ', side)
                print('    Side must be either "upper" or "lower')
                raise ValueError()
            
        cst_u, cst_l = cst_foil_fit(self.x, self.yu, self.x, self.yl, n_cst=self.n_cst)
        
        _, self.yu, self.yl, tmax, rLE = cst_foil(self.nn, cst_u, cst_l, x=self.x, t=tmax, tail=self.tail)

        return cst_u, cst_l, tmax, rLE

    def add_bump_to_thickness(self, xc: float, h: float, w: float, kind='H', 
                                keep_tmax=False) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Add a bump to the airfoil thickness.
        
        Parameters
        ----------
        xc : float
            Bump location
        
        h : float
            Bump height
        
        w : float
            Bump width
        
        kind : str
            Bump function type, i.e., 'G', 'H'.
            
        keep_tmax : bool
            Keep the maximum thickness of the airfoil
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        rLE : float
            Leading edge radius
        '''
        thickness = self.yu - self.yl
        
        if keep_tmax:
            tmax = np.max(thickness)
        else:
            tmax = None
        
        dy = bump_function(self.x, xc, h, w, kind)
        
        new_thickness = np.clip(thickness + dy, 0.0, None)
        dy = new_thickness - thickness
        
        self.yu = self.yu + 0.5*dy
        self.yl = self.yl - 0.5*dy
        
        cst_u, cst_l = cst_foil_fit(self.x, self.yu, self.x, self.yl, n_cst=self.n_cst)
        
        _, self.yu, self.yl, tmax, rLE = cst_foil(self.nn, cst_u, cst_l, x=self.x, t=tmax, tail=self.tail)

        return cst_u, cst_l, tmax, rLE

    def add_bump_to_camber(self, xc: float, h: float, w: float, kind='H',
                            keep_tmax=False) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Add a bump to the airfoil camber.
        
        Parameters
        ----------
        xc : float
            Bump location
        
        h : float
            Bump height
        
        w : float
            Bump width
        
        kind : str
            Bump function type, i.e., 'G', 'H'.
            
        keep_tmax : bool
            Keep the maximum thickness of the airfoil
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        rLE : float
            Leading edge radius
        '''
        if keep_tmax:
            tmax = np.max(self.yu - self.yl)
        else:
            tmax = None
            
        dy = bump_function(self.x, xc, h, w, kind)
        
        self.yu = self.yu + dy
        self.yl = self.yl + dy

        cst_u, cst_l = cst_foil_fit(self.x, self.yu, self.x, self.yl, n_cst=self.n_cst)
        
        _, self.yu, self.yl, tmax, rLE = cst_foil(self.nn, cst_u, cst_l, x=self.x, t=tmax, tail=self.tail)

        return cst_u, cst_l, tmax, rLE

    def add_cst_incremental_curves(self, increment_cst_u=None, increment_cst_l=None, 
                keep_tmax=False) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Add incremental CST curves to the airfoil.
        
        Parameters
        ----------
        increment_cst_u, increment_cst_l : ndarray or None
            Incremental CST coefficients for the upper or lower surface
            
        keep_tmax : bool
            Keep the maximum thickness of the airfoil
        '''
        if keep_tmax:
            tmax = np.max(self.yu - self.yl)
        else:
            tmax = None
        
        self.yu, self.yl = foil_increment(self.x, self.yu, self.yl, increment_cst_u, increment_cst_l, t=tmax)
        
        cst_u, cst_l = cst_foil_fit(self.x, self.yu, self.x, self.yl, n_cst=self.n_cst)
        
        _, self.yu, self.yl, tmax, rLE = cst_foil(self.nn, cst_u, cst_l, x=self.x, t=tmax, tail=self.tail)

        return cst_u, cst_l, tmax, rLE
        
    def set_leading_edge_radius(self, rLE_new: float, width_bump=0.8
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change airfoil leading edge radius to a new value (the maximum thickness is kept).
        
        Parameters
        ----------
        rLE_new : float
            New leading edge radius
            
        width_bump : float
            Width of the bump added to the airfoil thickness
            
        Returns
        -------
        cst_u, cst_l : ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        rLE : float
            Leading edge radius
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        rLE = geo.get_leading_edge_radius()
        
        n_try = 0
        
        while abs(rLE-rLE_new) > 1E-4 and n_try < self.MAX_TRY:
        
            cst_u, cst_l, tmax, rLE = self.add_bump_to_thickness(
                geo.X_RLE, 2.0*(rLE_new-rLE), width_bump, kind='H', keep_tmax=True)

            n_try += 1
            
        return cst_u, cst_l, tmax, rLE

    def set_leading_edge_slope_angle(self, slope_angle_new: float, width_bump=0.2
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change airfoil leading edge slope angle to a new value.
        
        Parameters
        ----------
        slope_angle_new : float
            New leading edge slope angle
            
        width_bump : float
            Width of the bump added to the airfoil thickness
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        slope_angle_new : float
            Leading edge slope angle (degree)
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        slope_angle = geo.get_leading_edge_slope_angle()
        
        n_try = 0
        
        while abs(slope_angle-slope_angle_new) > 1E-1 and n_try < self.MAX_TRY:
        
            dt = 0.5*np.deg2rad(slope_angle_new-slope_angle) * self.X_LE
        
            cst_u, cst_l, tmax, _ = self.add_bump_to_camber(
                self.X_LE, dt, width_bump, kind='G', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            slope_angle = geo.get_leading_edge_slope_angle()

            n_try += 1
            
        return cst_u, cst_l, tmax, slope_angle_new

    def set_trailing_edge_wedge_angle(self, wedge_angle_new: float, width_bump=0.2
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change airfoil trailing edge wedge angle to a new value.
        
        Parameters
        ----------
        wedge_angle_new : float
            New trailing edge wedge angle
            
        width_bump : float
            Width of the bump added to the airfoil thickness
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        wedge_angle_new : float
            Trailing edge wedge angle (degree)
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        wedge_angle = geo.get_trailing_edge_wedge_angle()
        
        n_try = 0
        
        while abs(wedge_angle-wedge_angle_new) > 1E-1 and n_try < self.MAX_TRY:
        
            dt = 0.5*np.deg2rad(wedge_angle_new-wedge_angle) * (1-self.X_TE)
        
            cst_u, cst_l, tmax, _ = self.add_bump_to_thickness(
                self.X_TE, dt, width_bump, kind='G', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            wedge_angle = geo.get_trailing_edge_wedge_angle()

            n_try += 1
            
        return cst_u, cst_l, tmax, wedge_angle_new

    def set_trailing_edge_slope_angle(self, slope_angle_new: float, width_bump=0.2
                ) -> Tuple[np.ndarray, np.ndarray, float, float]:
        '''
        Change airfoil trailing edge slope angle to a new value.
        
        Parameters
        ----------
        slope_angle_new : float
            New trailing edge slope angle
            
        width_bump : float
            Width of the bump added to the airfoil thickness
            
        Returns
        -------
        cst_u, cst_l : np.ndarray
            CST coefficients of the upper and lower surfaces
            
        tmax : float
            Maximum thickness of the airfoil
            
        slope_angle_new : float
            Trailing edge slope angle (degree)
        '''
        geo = FoilGeoFeatures(self.x, self.yu, self.yl)
        slope_angle = geo.get_trailing_edge_slope_angle()
        
        n_try = 0
        
        while abs(slope_angle-slope_angle_new) > 1E-1 and n_try < self.MAX_TRY:
        
            dt = 0.5*np.deg2rad(slope_angle_new-slope_angle) * (1-self.X_TE)
        
            cst_u, cst_l, tmax, _ = self.add_bump_to_camber(
                self.X_TE, dt, width_bump, kind='G', keep_tmax=True)
            
            geo = FoilGeoFeatures(self.x, self.yu, self.yl)
            slope_angle = geo.get_trailing_edge_slope_angle()

            n_try += 1
            
        return cst_u, cst_l, tmax, slope_angle_new


