'''
Math functions for the CST modeling package.
'''
from typing import Tuple, List, Union, Callable

import copy
import numpy as np

from scipy import spatial
from scipy.special import factorial
from scipy.interpolate import interp1d, CubicSpline
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

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

def find_circle_3p(p1, p2, p3) -> Tuple[float, np.ndarray]:
    '''
    Determine the radius and origin of a circle by 3 points (2D)
    
    Parameters
    -----------
    p1, p2, p3: list or ndarray [2]
        coordinates of points, [x, y]
        
    Returns
    ----------
    R: float
        radius
    XC: ndarray [2]
        circle center

    Examples
    ----------
    >>> R, XC = find_circle_3p(p1, p2, p3)

    '''

    # http://ambrsoft.com/TrigoCalc/Circle3D.htm

    A = p1[0]*(p2[1]-p3[1]) - p1[1]*(p2[0]-p3[0]) + p2[0]*p3[1] - p3[0]*p2[1]
    if np.abs(A) <= 1E-20:
        raise Exception('Finding circle: 3 points in one line')
    
    p1s = p1[0]**2 + p1[1]**2
    p2s = p2[0]**2 + p2[1]**2
    p3s = p3[0]**2 + p3[1]**2

    B = p1s*(p3[1]-p2[1]) + p2s*(p1[1]-p3[1]) + p3s*(p2[1]-p1[1])
    C = p1s*(p2[0]-p3[0]) + p2s*(p3[0]-p1[0]) + p3s*(p1[0]-p2[0])
    D = p1s*(p3[0]*p2[1]-p2[0]*p3[1]) + p2s*(p1[0]*p3[1]-p3[0]*p1[1]) + p3s*(p2[0]*p1[1]-p1[0]*p2[1])

    x0 = -B/2/A
    y0 = -C/2/A
    R  = np.sqrt(B**2+C**2-4*A*D)/2/np.abs(A)

    '''
    x21 = p2[0] - p1[0]
    y21 = p2[1] - p1[1]
    x32 = p3[0] - p2[0]
    y32 = p3[1] - p2[1]

    if x21 * y32 - x32 * y21 == 0:
        raise Exception('Finding circle: 3 points in one line')

    xy21 = p2[0]*p2[0] - p1[0]*p1[0] + p2[1]*p2[1] - p1[1]*p1[1]
    xy32 = p3[0]*p3[0] - p2[0]*p2[0] + p3[1]*p3[1] - p2[1]*p2[1]
    
    y0 = (x32 * xy21 - x21 * xy32) / 2 * (y21 * x32 - y32 * x21)
    x0 = (xy21 - 2 * y0 * y21) / (2.0 * x21)
    R = np.sqrt(np.power(p1[0]-x0,2) + np.power(p1[1]-y0,2))
    '''

    return R, np.array([x0, y0])

def angle_between_vectors(a: np.ndarray, b: np.ndarray, n=None, in_degree=True) -> float:
    '''
    Calculate the angle between two vectors.
    
    Parameters
    ------------
    a, b : ndarray [3]
        vectors
        
    n : ndarray [3]
        positive normal vector, by default None. 
        If None, the angle is in [0, 180] or [0, pi].
        
    in_degree : bool
        if True, return the angle in degree.
        
    Returns
    ---------
    angle : float
        angle between two vectors, in [-180, 180] or [-pi, pi].
    '''
    # Calculate the dot product of vectors a and b
    dot_product = np.dot(a, b)
    
    # Calculate the magnitudes of vectors a and b
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    
    # Calculate the cosine of the angle
    cos_theta = dot_product / (norm_a * norm_b)
    
    # Calculate the angle in radians in [0, pi]
    angle_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Clip to avoid numerical issues
    
    if in_degree:
        angle = np.rad2deg(angle_radians)
    else:
        angle = angle_radians
    
    # Calculate the cross product
    cross = np.cross(a, b)
    
    # Determine the direction of the angle using the sign of the z-component of the cross product
    if n is not None:
        if np.dot(cross, n) < 0:
            angle = - angle
            
    return angle

def project_vector_to_plane(v: np.ndarray, n: np.ndarray) -> np.ndarray:
    '''
    Project a vector `v` to a plane with normal vector `n`.
    
    Parameters
    ------------
    v : ndarray [3]
        vector to be projected.
        
    n : ndarray [3]
        normal vector of the plane.
    
    Returns
    ---------
    vp : ndarray [3]
        projected vector.
    '''
    n = n/np.linalg.norm(n)
    vp = v - np.dot(v, n)*n
    return vp


#* ===========================================
#* CST foils
#* ===========================================

def cst_foil(nn: int, cst_u, cst_l, x=None, t=None, tail=0.0, xn1=0.5, xn2=1.0, a0=0.0079, a1=0.96):
    '''
    Constructing upper and lower curves of an airfoil based on CST method

    CST: class shape transformation method (Kulfan, 2008)
    
    Parameters
    -----------
    nn: int
        total amount of points
    cst_u, cst_l: list or ndarray
        CST coefficients of the upper and lower surfaces
    x: ndarray [nn]
        x coordinates in [0,1] (optional)
    t: float
        specified relative maximum thickness (optional)
    tail: float
        relative tail thickness (optional)
    xn1, xn12: float
        CST parameters
        
    Returns
    --------
    x, yu, yl: ndarray
        coordinates
    t0: float
        actual relative maximum thickness
    R0: float
        leading edge radius
    
    Examples
    ---------
    >>> x_, yu, yl, t0, R0 = cst_foil(nn, cst_u, cst_l, x, t, tail)

    '''
    cst_u = np.array(cst_u)
    cst_l = np.array(cst_l)
    x_, yu = cst_curve(nn, cst_u, x=x, xn1=xn1, xn2=xn2, a0=a0, a1=a1)
    x_, yl = cst_curve(nn, cst_l, x=x, xn1=xn1, xn2=xn2, a0=a0, a1=a1)
    
    thick = yu-yl
    it = np.argmax(thick)
    t0 = thick[it]

    # Apply thickness constraint
    if t is not None:
        r  = (t-tail*x_[it])/t0
        t0 = t
        yu = yu * r
        yl = yl * r

    # Add tail
    for i in range(nn):
        yu[i] += 0.5*tail*x_[i]
        yl[i] -= 0.5*tail*x_[i]
        
    # Update t0 after adding tail
    if t is None:
        thick = yu-yl
        it = np.argmax(thick)
        t0 = thick[it]

    # Calculate leading edge radius
    x_RLE = 0.005
    yu_RLE = interp_from_curve(x_RLE, x_, yu)
    yl_RLE = interp_from_curve(x_RLE, x_, yl)
    R0, _ = find_circle_3p([0.0,0.0], [x_RLE,yu_RLE], [x_RLE,yl_RLE])

    return x_, yu, yl, t0, R0

def clustcos(i: int, nn: int, a0=0.0079, a1=0.96, beta=1.0) -> float:
    '''
    Point distribution on x-axis [0, 1]. (More points at both ends)
    
    Parameters
    ----------
    i: int
        index of current point (start from 0)
        
    nn: int
        total amount of points
        
    a0: float
        Parameter for distributing points near x=0.
        Smaller a0, more points near x=0.
        
    a1: float
        Parameter for distributing points near x=1.
        Larger a1, more points near x=1.
        
    beta: float
        Parameter for distribution points.

    Returns
    ---------
    float

    Examples
    ---------
    >>> c = clustcos(i, n, a0, a1, beta)

    '''
    aa = np.power((1-np.cos(a0*np.pi))/2.0, beta)
    dd = np.power((1-np.cos(a1*np.pi))/2.0, beta) - aa
    yt = i/(nn-1.0)
    a  = np.pi*(a0*(1-yt)+a1*yt)
    c  = (np.power((1-np.cos(a))/2.0,beta)-aa)/dd

    return c

def dist_clustcos(nn: int, a0=0.0079, a1=0.96, beta=1.0) -> np.ndarray:
    '''
    Point distribution on x-axis [0, 1]. (More points at both ends)

    Parameters
    ----------
    nn: int
        total amount of points
        
    a0: float
        Parameter for distributing points near x=0.
        Smaller a0, more points near x=0.
        
    a1: float
        Parameter for distributing points near x=1.
        Larger a1, more points near x=1.
        
    beta: float
        Parameter for distribution points.
    
    Examples
    ---------
    >>> xx = dist_clustcos(n, a0, a1, beta)

    '''
    aa = np.power((1-np.cos(a0*np.pi))/2.0, beta)
    dd = np.power((1-np.cos(a1*np.pi))/2.0, beta) - aa
    yt = np.linspace(0.0, 1.0, num=nn)
    a  = np.pi*(a0*(1-yt)+a1*yt)
    xx = (np.power((1-np.cos(a))/2.0,beta)-aa)/dd

    return xx

def cst_curve(nn: int, coef: np.array, x=None, xn1=0.5, xn2=1.0, a0=0.0079, a1=0.96) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Generating single curve based on CST method.

    CST: class shape transformation method (Kulfan, 2008)

    Parameters
    ----------
    nn: int
        total amount of points
    coef: ndarray
        CST coefficients
    x: ndarray [nn]
        coordinates of x distribution in [0,1] (optional)
    xn1, xn12: float
        CST parameters
    
    Returns
    --------
    x, y: ndarray
        coordinates
    
    Examples
    ---------
    >>> x, y = cst_curve(nn, coef, x, xn1, xn2)

    '''
    if x is None:
        x = dist_clustcos(nn, a0, a1)
    elif x.shape[0] != nn:
        raise Exception('Specified point distribution has different size %d as input nn %d'%(x.shape[0], nn))
    
    n_cst = coef.shape[0]

    s_psi = np.zeros(nn)
    for i in range(n_cst):
        xk_i_n = factorial(n_cst-1)/factorial(i)/factorial(n_cst-1-i)
        s_psi += coef[i] * xk_i_n * np.power(x, i) * np.power(1 - x, n_cst - 1 - i)

    C_n1n2 = np.power(x, xn1) * np.power(1 - x, xn2)
    y = C_n1n2 * s_psi
    y[0] = 0.0
    y[-1] = 0.0

    return x, y


#* ===========================================
#* Transformation
#* ===========================================

def transform(xu: np.ndarray, xl: np.ndarray, yu: np.ndarray, yl: np.ndarray, 
              scale=1.0, rot=None, x0=None, y0=None, xr=None, yr=None, dx=0.0, dy=0.0, 
              projection=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
    Apply chord length, twist angle(deg) and leading edge position to a 2D curve.
    
    The transformation is applied in the following order:
    
    1. Translation
    2. Scaling
    3. Rotation

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
        xu_new, yu_new, _ = rotate(xu_new, yu_new, np.zeros_like(xu_new), angle=rot, origin=[xr, yr, 0.0], axis='Z')
        xl_new, yl_new, _ = rotate(xl_new, yl_new, np.zeros_like(xu_new), angle=rot, origin=[xr, yr, 0.0], axis='Z')

    return xu_new, xl_new, yu_new, yl_new

def transform_curve(xx: np.ndarray, yy: np.ndarray, dx=0.0, dy=0.0, dz=0.0,
                    scale=1.0, x0=None, y0=None, 
                    rot_z=0.0, rot_x=0.0, rot_y=0.0, rot_axis=0.0,
                    xr=None, yr=None, zr=None
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Transform a 2D (unit) curve to a 3D curve by translation, scaling and rotation.
    
    The transformation is applied in the following order:
    
    1. Translation
    2. Scaling
    
        - Scale center: (x0, y0), the first point of the curve by default.
        - Scale factor.
    
    3. Rotation
    
        - Rotate about the z axis by `rot_z` degree.
        - Rotate about the x axis by `rot_x` degree.
        - Rotate about the y axis by `rot_y` degree.
        - Rotate about the main axis of the curve by `rot_axis` degree.
        - Rotate center: (xr, yr, zr), the scale center by default.

    Parameters
    -------------
    xx, yy : ndarray
        a 2D (unit) curve.
        
    dx, dy : float
        translation vector, e.g., leading edge location.
        
    scale : bool
        scale factor.
        
    x0, y0 : float
        the scale center for the 2D curve in the x-y plane.
        
    rot_x, rot_y, rot_z : float
        rotate angle (degree) about the x, y, z axis.
        
    rot_axis : float
        rotate angle (degree) about the main axis of the curve,
        e.g., the chord line of an airfoil.

    xr, yr, zr : float
        the rotation center for the 2D curve

    Returns
    ---------
    x, y, z : ndarray
        coordinates of the 3D curve.
    '''
    #* Translation
    x = dx + xx
    y = dy + yy
    z = dz + np.zeros_like(x)

    #* Scale center
    x0 = x0 if x0 is not None else x[0]
    y0 = y0 if y0 is not None else y[0]
        
    #* Scaling
    x = x0 + (x-x0)*scale
    y = y0 + (y-y0)*scale

    #* Rotation center
    xr = xr if xr is not None else x0
    yr = yr if yr is not None else y0
    zr = zr if zr is not None else dz
    
    #* Rotation
    xv = [1.0]; yv = [0.0]; zv = [0.0]
    
    if abs(rot_z) > 1.0E-12:
        x,  y,  z  = rotate(x,  y,  z,  angle=rot_z, origin=[xr, yr, zr], axis='Z')
        xv, yv, zv = rotate(xv, yv, zv, angle=rot_z, axis='Z')

    if abs(rot_x) > 1.0E-12:
        x,  y,  z  = rotate(x,  y,  z,  angle=rot_x, origin=[xr, yr, zr], axis='X')
        xv, yv, zv = rotate(xv, yv, zv, angle=rot_z, axis='X')

    if abs(rot_y) > 1.0E-12:
        x,  y,  z  = rotate(x,  y,  z,  angle=rot_y, origin=[xr, yr, zr], axis='Y')
        xv, yv, zv = rotate(xv, yv, zv, angle=rot_z, axis='Y')
        
    if abs(rot_axis) > 1.0E-12:
        points = rotate_vector(x, y, z, angle=rot_axis, origin=[xr, yr, zr], axis_vector=[xv[0], yv[0], zv[0]])
        x = points[:,0]
        y = points[:,1]
        z = points[:,2]
        
    return x, y, z

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
    if axis in 'X':
        axis_vector=[1,0,0]
    if axis in 'Y':
        axis_vector=[0,1,0]
    if axis in 'Z':
        axis_vector=[0,0,1]

    points = rotate_vector(x, y, z, angle=angle, origin=origin, axis_vector=axis_vector)
    
    x_ = points[:,0]
    y_ = points[:,1]
    z_ = points[:,2]

    return x_, y_, z_

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

    rot = Rotation.from_rotvec(angle*rotation_vector/180.0*np.pi)
    
    # In terms of rotation matrices, this application is the same as rot.as_matrix().dot(vector).
    points = rot.apply(vector) + origin
    
    return points

def rotation_3d(pp: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float):
    '''
    The rotation_3d is derived from Chenyu Wu. 2022. 11. 5

    ### Description
    This function rotate a set of points based on the origin and the axis given by the inputs

    ### Inputs
    `pp`: The point set that is going to be rotated. `pp.shape = (n_points, 3)`

    `origin`: The numpy array that defines the origin of the rotation axis. The shape must be `(3,0)`

    `axis`: The direction of the rotation axis. This axis does not need to be normalized. The shape must be `(3,0)`

    `angle`: The rotation angle in degree

    ### Outputs
    `xnew, ynew, znew`: The rotated points.
    '''
    # Translate the points to a coordinate system that has the origin defined by the input
    # The points have to be translated back to the original frame before return.
    
    nn = pp.shape[0]
    for i in range(nn):
        pp[i, :] = pp[i, :] - origin
    xnew, ynew, znew = np.zeros(nn), np.zeros(nn), np.zeros(nn)
    
    norm = np.sqrt(axis @ axis)
    if norm < 1e-8:
        raise Exception("The length of the axis is too short!")
    e3 = axis / norm

    angle_rad = np.pi * angle / 180.0

    for i in range(nn):
        vec = pp[i, :].copy()

        # compute the parallel component
        vec_p = (vec @ e3) * e3
        # compute the normal component
        vec_n = vec - vec_p

        # define the local coordinate system
        e1 = vec_n
        e2 = np.cross(e3, e1)

        # rotate
        vec_n_rot = e1 * np.cos(angle_rad) + e2 * np.sin(angle_rad)

        # assemble the vector
        vec_new = vec_n_rot + vec_p
        xnew[i], ynew[i], znew[i] = vec_new[0], vec_new[1], vec_new[2]

    # transform back to the original frame
    xnew, ynew, znew = xnew + origin[0], ynew + origin[1], znew + origin[2]

    pp_new = np.hstack((xnew.reshape(-1,1), ynew.reshape(-1,1), znew.reshape(-1,1)))
    # print(pp_new.shape)
    return pp_new

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
                 flip=True, origin=[0, 0]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Bend the cylinder curve to a 2D plane curve.

    Parameters
    ----------
    x, y ,z : ndarray
        coordinates of the curve on a cylinder. `x` and `y` must not be 0 at the same time.
    flip : bool
        if True, flip the X of the extracted plane curve.
    origin: array_like
        the cylinder origin, [x0, y0] (or [x0, y0, 0]).

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

    x = x - origin[0]
    y = y - origin[1]

    rr = np.sqrt(x*x+y*y)
    tt = np.arctan2(y, x) * coef

    X = rr*tt
    Y = z.copy()
    Z = rr

    return X, Y, Z

def toCylinder(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, 
               flip=True, origin=[0, 0]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Bend the plane sections to curves on a cylinder.

    Parameters
    ----------
    X, Y, Z : ndarray
        coordinates of the curve on a plane. `Z` must not be 0.
    flip : bool
        if True, flip the X of the extracted plane curve.
    origin: array_like
        the cylinder origin, [x0, y0] (or [x0, y0, 0]).

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

    x = x + origin[0]
    y = y + origin[1]

    return x, y, z


#* ===========================================
#* Interpolation
#* ===========================================

def interp_from_curve(x0: Union[float, np.ndarray], x: np.ndarray, y: np.ndarray, 
                        extrapolate=False) -> Union[float, np.ndarray]:
    '''
    Interpolate points from curve represented points [x, y].
    
    Parameters
    ----------
    x0 : Union[float, np.ndarray]
        ndarray/value of x locations to be interpolated.
        
    x, y : ndarray
        coordinates of the curve.

    Returns
    ----------
    y0 : Union[float, np.ndarray]
        interpolated coordinates

    Examples
    ---------
    >>> y0 = interp_from_curve(x0, x, y)
    '''
    if extrapolate:
        f  = interp1d(x, y, kind='cubic', fill_value='extrapolate')
    else:
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

def plane_3points(P0: np.ndarray, P1: np.ndarray, P3: np.ndarray, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    '''
    Calculate the plane function and the coordinates of given points (`xs`, `ys`).
    
    The plane function is `a*x+b*y+c*z+d=0`.
    
    Parameters
    -------------
    P0, P1, P3 : ndarray [3]
        coordinates of three points of plane P0123.
    xs, ys : ndarray [:] or [:,:]
        X and Y coordinates of plane points.
    
    Returns
    -------
    zs : ndarray [:] or [:,:]
    
    Examples
    ---------
    >>> xs = plane_3points(P0, P1, P3, xs, ys)
    '''
    a1 = P1[0] - P0[0]
    b1 = P1[1] - P0[1]
    c1 = P1[2] - P0[2]
    a2 = P3[0] - P0[0]
    b2 = P3[1] - P0[1]
    c2 = P3[2] - P0[2]
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * P0[0] - b * P0[1] - c * P0[2])
    
    if c == 0:
        return np.zeros_like(xs)
    else:
        return (a*xs+b*ys+d)/(-c)
    
    
#* ===========================================
#* Intersection
#* ===========================================
    
def intersect_index(x1, y1, x2, y2):
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
    >>> i1, i2, points = intersect_index(x1, y1, x2, y2)

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

def join_curves(curves: List[np.ndarray], cri_dup=1e-6) -> np.ndarray:
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
    
def extract_slice(data: List[np.ndarray], locations: List[float], Pref: np.ndarray, dir_norm: np.ndarray, dir_ref=np.array([1.,0.,0.]),
                    zone_id=[], index_xyz=[0,1,2], index_var=None, arrange_method='join'):
    '''
    Extract data sliced by planes.
    
    Parameters
    --------------
    data : List[ndarray]  [n_zone][ni,nj,nk,nv]
        data of all surfaces. 
        
    locations : List[float]
        list of distances to the reference point in the given direction.
        
    Pref : ndarray [3]
        coordinates of the reference point.
        
    dir_norm : ndarray [3]
        direction vector normal to the slice plane (will be normalized).
        
    dir_ref : ndarray [3]
        direction vector that roughly sets the xi-axis in the slice plane.
        
    fname : str
        file name.
        
    zone_id : List[int]
        index of zones in the tecplot format file, start from 0.
        
    index_xyz : List[int]
        index of variables in the data for X, Y and Z.
        
    index_var : List[int] or None
        index of variables of interest in the data.
        
    arrange_method : str
        if 'join', keeps the original order of points (suitable for surface with a few blocks).
        If 'rearrange', rearrange points by minimal distance.
    
    Returns
    ------------
    sections : list of ndarray [:,3+nv]
        coordinates and data on the slice.

    
    Examples
    ----------
    >>> sections = extract_slice(data, locations, Pref, dir_norm, dir_ref=np.array([1.,0.,0.]),
                    zone_id=[], index_xyz=[0,1,2], arrange_method='join')
    '''
    if index_var is None:
        index_var = []
        for i in range(data[0].shape[3]):
            if i not in index_xyz:
                index_var.append(i)
    
    if len(zone_id)>0:
        data = [data[i] for i in zone_id]
    
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

    return sections

    
#* ===========================================
#* Functions
#* ===========================================

def smooth_omega_shape_function(x: np.ndarray, c0=0.1, c1=0.9, b0=50, b1=50) -> np.ndarray:
    '''
    Smooth Omega-shape function, the output equals 0 at both ends, 1 at the middle.
    There can be plateau at the middle and both ends. The function is smooth.
    
    Parameters
    ----------
    x : ndarray
        input values, range in [0,1].
        
    c0, c1 : float
        sharp transition location at both ends.
        
    b0, b1 : float
        sharpness of the sigmoid at both ends.
        `b` = 1  gives an almost A shape function, 
        `b` = 25 gives a sigmoid transition with width about 0.1,
        `b` = 50 gives a sigmoid transition with width about 0.05.     
        
    Returns
    ---------
    y : ndarray
        output values, scaled to [0,1].
        The output equals 0 at both ends, 1 at the middle.
    '''
    #* Scale x to [-1,1], and translate by c0, c1, then flip x1
    
    if c0 > 0:
        x0 = 2*(x-c0)
        r0 = scaled_sigmoid(x0, b0)
    else:
        r0 = 0.0
        
    if c1 < 1:
        x1 = 2*(x-c1); x1 = -x1
        r1 = scaled_sigmoid(x1, b1)
    else:
        r1 = 0.0
    
    y = r0 + r1

    return (y-np.min(y))/(np.max(y)-np.min(y))

def scaled_sigmoid(x: np.ndarray, b=1) -> np.ndarray:
    '''
    Scaled Sigmoid function. 
    
    Parameters
    ----------
    x : ndarray
        input values.
        The sharp transition is at `x=0`, therefore, `x` should be in [-1,1].
        
    b : float
        sharpness of the sigmoid. A larger `b` gives a steeper sigmoid.
        When `x` in [-1,1], `b` = 1 gives a almost linear function, 
        `b` = 10 gives a sigmoid transition with width about 1.0,
        `b` = 50 gives a sigmoid transition with width about 0.1.
        
    Returns
    ---------
    y : ndarray
        output values, scaled to [0,1], i.e., `y` at `x[0]` is 0, `y` at `x[-1]` is 1.    
    '''
    y = 1.0/(1.0+np.exp(-b*x))
    return (y-np.min(y))/(np.max(y)-np.min(y))


class CoordinateTransformation():
    '''
    Transform (x) coordinates to another (x') coordinates, i.e.,    
    `x' = f(x)`, `x` in `[0,1]`, `x'` in `[0,1]`.
    
    '''
    
    def __init__(self):
        
        self.func : Callable[[np.ndarray], np.ndarray] = None

    def set_function(self, func: Callable[[np.ndarray], np.ndarray]) -> None:
        '''
        Set the transformation function.
        
        Parameters
        ----------
        func : Callable[[np.ndarray], np.ndarray]
            transformation function.
        '''
        self.func = func
    
    def transform(self, x: np.ndarray) -> np.ndarray:
        '''
        Transform x to x'.
        
        Parameters
        ----------
        x : ndarray
            input values.
        
        Returns
        ---------
        x' : ndarray
            transformed values.
        '''
        if self.func is None:
            raise Exception('Transformation function is not set.')
        else:
            return self.func(x)
        
    def set_function_by_interpolation(self, x: List[float], xp: List[float],
            slope0=None, slope1=None) -> None:
        '''
        Set the transformation function by interpolation.
        
        Parameters
        ----------
        x, xp : list of float
            coordinates of the points to be interpolated.
            
        slope0, slope1 : float or None
            slope at the two ends, i.e., d(x')/dx at x=0 and x=1.
        '''

        if isinstance(slope0, float):
            bc0 = (1, slope0)
        else:
            bc0 = (2, 0.0)
            
        if isinstance(slope1, float):
            bc1 = (1, slope1)
        else:
            bc1 = (2, 0.0)

        x0 = [0.0] + x + [1.0]
        x1 = [0.0] + xp + [1.0]
        
        self.func = CubicSpline(x0, x1, bc_type=(bc0, bc1))


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    plt.figure(figsize=(16,4))
    
    xx = np.linspace(-1, 1, 1001, endpoint=True)
    
    plt.subplot(1,2,1)
    plt.title('Scaled Sigmoid function')
    
    yy = scaled_sigmoid(xx, b=1)
    plt.plot(xx, yy, 'k', label='b=1')
    
    yy = scaled_sigmoid(xx, b=10)
    plt.plot(xx, yy, 'g', label='b=10')
    
    yy = scaled_sigmoid(xx, b=50)
    plt.plot(xx, yy, 'r', label='b=50')
    
    
    plt.legend()
    
    xx = np.linspace(0, 1, 1001, endpoint=True)

    plt.subplot(1,2,2)
    plt.title('Smooth Ratio function')
    
    yy = smooth_omega_shape_function(xx, c0=0.1, c1=0.9, b0=50, b1=50)
    plt.plot(xx, yy, 'k', label='c0=0.1, c1=0.9, b0=50, b1=50')
    
    yy = smooth_omega_shape_function(xx, c0=0.0, c1=0.8, b0=25, b1=25)
    plt.plot(xx, yy, 'g', label='c0=0.0, c1=0.8, b0=25, b1=25')
    
    yy = smooth_omega_shape_function(xx, c0=0.1, c1=0.9, b0=1, b1=1)
    plt.plot(xx, yy, 'r', label='c0=0.1, c1=0.9, b0=1, b1=1')
    
    plt.legend()
    
    plt.show()
    
    
    xShift = CoordinateTransformation()
    xShift.set_function_by_interpolation(x=[0.6], xp=[0.4], slope0=1.0, slope1=None)
    
    xp = xShift.transform(xx)
    plt.plot(xx, xp, 'k')
    plt.xlabel('x')
    plt.ylabel('x\'')
    plt.axis('equal')
    plt.show()
    
    
