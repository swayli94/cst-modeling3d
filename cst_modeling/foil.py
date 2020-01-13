'''
This is a module containing functions to construct an airfoil
'''
import copy
import numpy as np 
from numpy.linalg import lstsq
from scipy.special import factorial
from scipy.interpolate import interp1d

class Section:
    '''
    Section curve (3D) generated by CST foil
    '''
    def __init__(self, thick=None, chord=1.0, twist=0.0, tail=0.0):
        self.xLE = 0.0
        self.yLE = 0.0
        self.zLE = 0.0
        self.chord = chord
        self.twist = twist
        self.thick = thick
        self.tail = tail
        self.RLE = 0.0

        #* 2D unit airfoil
        self.cst_u = []
        self.cst_l = []
        self.xx = []
        self.yu = []
        self.yl = []

        #* 3D section
        self.x = []
        self.y = []
        self.z = []

    def set_params(self, thick=None, chord=1.0, twist=0.0, tail=0.0, xLE=0.0, yLE=0.0, zLE=0.0):
        '''
        Set parameters of the section
        '''
        self.xLE = xLE
        self.yLE = yLE
        self.zLE = zLE
        self.chord = chord
        self.twist = twist
        self.thick = thick
        self.tail = tail

    def foil(self, cst_u=None, cst_l=None, nn=1001):
        '''
        Generating the section (3D) by cst_foil. 
            nn:     total amount of points
            cst_u:  CST coefficients of upper surface (list, optional)
            cst_l:  CST coefficients of lower surface (list, optional)
        '''
        if not cst_u is None and not cst_l is None:
            self.cst_u = copy.deepcopy(cst_u)
            self.cst_l = copy.deepcopy(cst_l)

        self.xx, self.yu, self.yl, self.thick, self.RLE = cst_foil(
            nn, self.cst_u, self.cst_l, t=self.thick, tail=self.tail)

        xx_, yu_, yl_ = transform(self.xx, self.yu, self.yl, 
            scale=self.chord, rotate=self.twist, dx=self.xLE, dy=self.yLE, proj=True)

        self.x = []
        self.y = []
        self.z = []
        for i in range(nn):
            self.x.append(xx_[-1-i])
            self.y.append(yl_[-1-i])
            self.z.append(self.zLE)

        for i in range(1,nn):
            self.x.append(xx_[i])
            self.y.append(yu_[i])
            self.z.append(self.zLE)

    def copyfrom(self, other):
        '''
        Copy from anthor section object
        '''
        if not isinstance(other, Section):
            raise Exception('Must copy from another section object')
        
        self.xLE = other.xLE
        self.yLE = other.yLE
        self.zLE = other.zLE
        self.chord = other.chord
        self.twist = other.twist
        self.thick = other.thick
        self.tail = other.tail
        self.RLE = other.RLE

        self.cst_u = copy.deepcopy(other.cst_u)
        self.cst_l = copy.deepcopy(other.cst_l)
        self.xx = copy.deepcopy(other.xx)
        self.yu = copy.deepcopy(other.yu)
        self.yl = copy.deepcopy(other.yl)

        self.x = copy.deepcopy(other.x)
        self.y = copy.deepcopy(other.y)
        self.z = copy.deepcopy(other.z)

#TODO: ===========================================
#TODO: Static functions
#TODO: ===========================================

def cst_foil(nn, coef_upp, coef_low, x=None, t=None, tail=0.0):
    '''
    Constructing upper and lower curves of an airfoil based on CST method
        nn:         total amount of points
        coef_upp:   CST coefficients of upper surface (list)
        coef_low:   CST coefficients of lower surface (list)
        x:          point x [0,1] (optional list, size is nn)
        t:          relative maximum thickness (optional)
        tail:       relative tail thickness (optional)

        CST:    class shape transfermation method (Kulfan, 2008)

    Return lists of x, y_upp, y_low
    '''
    x_, yu_ = cst_curve(nn, coef_upp, x=x)
    x_, yl_ = cst_curve(nn, coef_low, x=x)
    yu = np.array(yu_)
    yl = np.array(yl_)
    
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

    # Calculate leading edge radius
    x_RLE = 0.005
    yu_RLE = interplot_from_curve(x_RLE, x_, yu)
    yl_RLE = interplot_from_curve(x_RLE, x_, yl)
    R0, _ = find_circle_3p([0.0,0.0], [x_RLE,yu_RLE], [x_RLE,yl_RLE])

    return x_, yu.tolist(), yl.tolist(), t0, R0

def cst_foil_fit(xu, yu, xl, yl, n_order=7):
    '''
    Using CST method to fit an airfoil
        xu, yu:  upper surface points (list)
        xl, yl:  lower surface points (list)
        n_order: number of CST parameters

    Return: coef_upp, coef_low (list)
    '''
    coef_upp = fit_curve(xu, yu, n_order=n_order)
    coef_low = fit_curve(xl, yl, n_order=n_order)
    return coef_upp, coef_low

def foil_bump_modify(x, yu, yl, xc, h, s, side, n_order=None):
    '''
    Add bumps on the airfoil
        x, yu, yl: current airfoil (list)
        xc:        x of the bump center
        h:         relative height of the bump (to maximum thickness)
        s:         span of the bump
        side:      +1/-1 upper/lower side of the airfoil
        n_order:   if specified (>0), then use CST to fit the new foil

    Return: yu_new, yl_new
    '''
    yu_new = copy.deepcopy(yu)
    yl_new = copy.deepcopy(yl)

    yu_ = np.array(yu)
    yl_ = np.array(yl)
    t0 = np.max(yu_-yl_)

    if isinstance(xc, list):

        for i in range(len(xc)):
            if xc[i]<0.1 or xc[i]>0.9:
                kind = 'H'
            else:
                kind = 'G'

            if side[i] > 0:
                yu_new = add_bump(x, yu_new, xc[i], h[i]*t0, s[i], kind=kind)
            else:
                yl_new = add_bump(x, yl_new, xc[i], h[i]*t0, s[i], kind=kind)

            yu_ = np.array(yu_new)
            yl_ = np.array(yl_new)
            it = np.argmax(yu_-yl_)
            tu = np.abs(yu_new[it])
            tl = np.abs(yl_new[it])

            if side[i] > 0:
                rl = (t0-tu)/tl
                yl_new = (rl * np.array(yl_new)).tolist()
            else:
                ru = (t0-tl)/tu
                yl_new = (ru * np.array(yu_new)).tolist()

    else:
        if xc<0.1 or xc>0.9:
            kind = 'H'
        else:
            kind = 'G'

        if side > 0:
            yu_new = add_bump(x, yu_new, xc, h*t0, s, kind=kind)
        else:
            yl_new = add_bump(x, yl_new, xc, h*t0, s, kind=kind)

        yu_ = np.array(yu_new)
        yl_ = np.array(yl_new)
        it = np.argmax(yu_-yl_)
        tu = np.abs(yu_new[it])
        tl = np.abs(yl_new[it])

        if side > 0:
            rl = (t0-tu)/tl
            yl_new = (rl * np.array(yl_new)).tolist()
        else:
            ru = (t0-tl)/tu
            yl_new = (ru * np.array(yu_new)).tolist()

    if not n_order is None:
        # CST reverse
        tail = yu[-1] - yl[-1]

        c_upp, c_low = cst_foil_fit(x, yu_new, x, yl_new, n_order=n_order)
        _, yu_new, yl_new, _, _ = cst_foil(len(x),c_upp,c_low,x=x, t=t0,tail=tail)

    return yu_new, yl_new

def check_foil(x, yu, yl):
    '''
    Check if the airfoil is reasonable.
    Calculate curvature, thickness, camber distribution.
        x, yu, yl: current airfoil (list)

    Return: curv_u, curv_l, thickness, camber
    '''
    curv_u = curve_curvature(x, yu)
    curv_l = curve_curvature(x, yl)

    thickness = []
    camber = []
    for i in range(len(x)):
        tt = yu[i] - yl[i]
        thickness.append(tt)
        camber.append(0.5*(yu[i]+yl[i]))

        if tt < 0:
            print('Unreasonable Airfoil: negative thickness')

    return curv_u, curv_l, thickness, camber

def transform(x, yu, yl, scale=1.0, rotate=None, x0=None, y0=None, dx=0.0, dy=0.0, proj=False):
    '''
    Apply chord length, twist angle(deg) and leading edge position to unit airfoil
        x, yu, yl:  current curve or unit airfoil (list)
        scale:      scale factor, e.g., chord length
        rotate:     rotate angle (deg), +z direction for x-y plane, e.g., twist angle
        x0, y0:     rotation center (scaler)
        dx, dy:     translation, e.g., leading edge location
        proj:       True => for unit airfoil, the rotation keeps the projection length the same

    Return: x_new, yu_new, yl_new
    '''
    #* Translation
    x_new  = dx + np.array(copy.deepcopy(x)) 
    yu_new = dy + np.array(copy.deepcopy(yu))
    yl_new = dy + np.array(copy.deepcopy(yl))

    #* Rotation center
    if x0 is None:
        x0 = x_new[0]
    if y0 is None:
        y0 = 0.5*(yu_new[0]+yl_new[0])
    
    #* Scale (keeps the same projection length)
    rr = 1.0
    if proj and not rotate is None:
        angle = rotate/180.0*np.pi  # rad
        rr = np.cos(angle)

    x_new  = x0 + (x_new -x0)*scale/rr
    yu_new = y0 + (yu_new-y0)*scale/rr
    yl_new = y0 + (yl_new-y0)*scale/rr

    #* Rotation
    if not rotate is None:
        angle  = rotate/180.0*np.pi  # rad
        for i in range(1,len(x)):
            aa = x_new [i] - x0
            bu = yu_new[i] - y0
            bl = yl_new[i] - y0
            lu = np.sqrt(aa*aa+bu*bu)
            ll = np.sqrt(aa*aa+bl*bl)
            tu = np.arctan(bu/aa)
            tl = np.arctan(bl/aa)

            x_new [i] = x0 + 0.5*(lu*np.cos(tu-angle) + ll*np.cos(tl-angle))
            yu_new[i] = y0 + lu*np.sin(tu-angle)
            yl_new[i] = y0 + ll*np.sin(tl-angle)

    return x_new.tolist(), yu_new.tolist(), yl_new.tolist()

def rotate(x, y, z, angle=0.0, origin=[0.0, 0.0, 0.0], axis='X'):
    '''
    Rotate the 3D curve according to origin
        x,y,z:  curve lists
        angle:  rotation angle (deg)
        origin: rotation origin
        axis:   rotation axis (use positive direction to define angle)

    Return:
        x_, y_, z_
    '''
    cc = np.cos( angle/180.0*np.pi )
    ss = np.sin( angle/180.0*np.pi )
    x_ = copy.deepcopy(x)
    y_ = copy.deepcopy(y)
    z_ = copy.deepcopy(z)

    if axis in 'X':
        for i in range(len(x)):
            y_[i] = origin[1] + (y[i]-origin[1])*cc - (z[i]-origin[2])*ss
            z_[i] = origin[2] + (y[i]-origin[1])*ss + (z[i]-origin[2])*cc

    if axis in 'Y':
        for i in range(len(x)):
            z_[i] = origin[2] + (z[i]-origin[2])*cc - (x[i]-origin[0])*ss
            x_[i] = origin[0] + (z[i]-origin[2])*ss + (x[i]-origin[0])*cc

    if axis in 'Z':
        for i in range(len(x)):
            x_[i] = origin[0] + (x[i]-origin[0])*cc - (y[i]-origin[1])*ss
            y_[i] = origin[1] + (x[i]-origin[0])*ss + (y[i]-origin[1])*cc

    return x_, y_, z_

#TODO: ===========================================
#TODO: Supportive functions
#TODO: ===========================================

def clustcos(i, nn, a0=0.0079, a1=0.96, beta=1.0):
    '''
    Point distribution on x-axis [0, 1]. (More points at both ends)
        i:      index of current point (start from 0)
        nn:     total amount of points
        a0:     parameter for distributing points near x=0
        a1:     parameter for distributing points near x=1
        beta:   parameter for distribution points 
    '''
    aa = np.power((1-np.cos(a0*np.pi))/2.0, beta)
    dd = np.power((1-np.cos(a1*np.pi))/2.0, beta) - aa
    yt = i/(nn-1.0)
    a  = np.pi*(a0*(1-yt)+a1*yt)
    c  = (np.power((1-np.cos(a))/2.0,beta)-aa)/dd

    return c

def cst_curve(nn, coef, x=None, xn1=0.5, xn2=1.0):
    '''
    Generating single curve based on CST method.
        nn:     total amount of points
        coef:   CST coefficients (list)
        x:      points x [0,1] (optional list, size= nn)
        xn1,2:  CST parameters

        CST:    class shape transfermation method (Kulfan, 2008)

    Return lists of x, y
    '''
    if x is None:
        x = []
        for i in range(nn):
            x.append(clustcos(i, nn))
    elif len(x) != nn:
        raise Exception('Specified point distribution has different size %d as input nn %d'%(len(x), nn))
    
    n_order = len(coef)
    y = []
    for ip in range(nn):
        s_psi = 0.0
        for i in range(n_order):
            xk_i_n = factorial(n_order-1)/factorial(i)/factorial(n_order-1-i)
            s_psi += coef[i]*xk_i_n * np.power(x[ip],i) * np.power(1-x[ip],n_order-1-i)

        C_n1n2 = np.power(x[ip],xn1) * np.power(1-x[ip],xn2)
        y.append(C_n1n2*s_psi)

    y[0] = 0.0
    y[-1] = 0.0

    return x, y

def find_circle_3p(p1, p2, p3):
    '''
    Determine the radius and origin of a circle by 3 points (2D)
        p1, p2, p3: lists of [x, y]

    Return: R, [xc, yc]
    '''
    x21 = p2[0] - p1[0]
    y21 = p2[1] - p1[1]
    x32 = p3[0] - p2[0]
    y32 = p3[1] - p2[1]

    if x21 * y32 - x32 * y21 == 0:
        print('Finding circle: 3 points in one line')
        return None

    xy21 = p2[0]*p2[0] - p1[0]*p1[0] + p2[1]*p2[1] - p1[1]*p1[1]
    xy32 = p3[0]*p3[0] - p2[0]*p2[0] + p3[1]*p3[1] - p2[1]*p2[1]
    
    y0 = (x32 * xy21 - x21 * xy32) / 2 * (y21 * x32 - y32 * x21)
    x0 = (xy21 - 2 * y0 * y21) / (2.0 * x21)
    R = np.sqrt(np.power(p1[0]-x0,2) + np.power(p1[1]-y0,2))

    return R, [x0, y0]

def curve_curvature(x, y):
    '''
    Calculate curvature of points in the curve
        x, y: points of curve (list or ndarray)

    Return: curv (list)
    '''
    nn = len(x)
    if nn<3:
        raise Exception('curvature needs at least 3 points')
    
    curv = [0.0]
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

        curv.append(curv_)

    curv[0] = curv[1]
    curv.append(curv[-1])

    return curv

def interplot_from_curve(x0, x, y):
    '''
    Interplot points from curve represented points [x, y]
        x0  : list/number of x locations to be interploted
        x, y: points of curve (list or ndarray)

    Return: y0 list/number
    '''
    x_ = np.array(x)
    y_ = np.array(y)
    f = interp1d(x_, y_, kind='cubic')
    x0_ = np.array(x0)

    if isinstance(x0, list):
        y0 = f(x0_).tolist()
    elif np.size(x0_)==1:
        y0 = f(x0_)
    else:
        raise Exception('Interplot: x0 must be a list or number')

    return y0

def add_bump(x, y, xc, h, s, kind='G'):
    '''
    Add a bump on current curve [x, y]
        x, y:   current curve (list, x[0,1])
        xc:     x of the bump center
        h:      height of the bump
        s:      span of the bump
        kind:   bump function
            G: Gaussian, less cpu cost
            H: Hicks-Henne, better when near leading edge

    Return: y_new (list, new curve)
    '''
    if xc<=0 or xc>=1:
        print('Bump location not valid (0,1): xc = %.3f'%(xc))
        return y

    if 'G' in kind:
        y_new = []
        for i in range(len(x)):
            if xc-s<0.0 and x[i]<xc:
                sigma = xc/3.5
            elif  xc+s>1.0 and x[i]>xc:
                sigma = (1.0-xc)/3.5
            else:
                sigma = s/6.0
            aa = -np.power(x[i]-xc,2)/2.0/sigma**2
            y_new.append(y[i]+h*np.exp(aa))

    else:
        
        s0 = np.log(0.5)/np.log(xc) 

        Pow = 1
        span = 1.0
        hm = np.abs(h)
        while Pow<100 and span>s:
            x1  = -1.0
            x2  = -1.0
            for i in range(0, 201):
                xx = i*0.005
                rr = np.pi*np.power(xx,s0)
                yy = hm * np.power(np.sin(rr),Pow)
                if yy > 0.01*hm and x1<0.0 and xx<xc:
                    x1 = xx
                if yy < 0.01*hm and x2<0.0 and xx>xc:
                    x2 = xx
            if x2 < 0.0:
                x2 = 1.0
            
            span = x2 - x1
            Pow = Pow + 1

        y_new = []
        for i in range(len(x)):
            rr = np.pi*np.power(x[i],s0)
            dy = h*np.power(np.sin(rr),Pow)
            y_new.append(y[i]+dy)

    return y_new

def fit_curve(x, y, n_order=7, xn1=0.5, xn2=1.0):
    '''
    Using least square method to fit a CST curve
        x, y:    curve points (list)
        n_order: number of CST parameters

    Array A: A[nn, n_order], nn=len(x)
    Array b: b[nn]

    Return: coef (list)
    '''
    nn = len(x)
    L  = x[-1] - x[0]
    x_ = []
    b_ = []
    for ip in range(nn):
        x_.append((x[ip]-x[0])/L)    # scaling x to 0~1
        b_.append(y[ip]-x_[ip]*y[-1])  # removing tail
    b = np.array(b_)

    A = np.zeros((nn, n_order))
    for ip in range(nn):
        C_n1n2 = np.power(x_[ip],xn1) * np.power(1-x_[ip],xn2)
        for i in range(n_order):
            xk_i_n = factorial(n_order-1)/factorial(i)/factorial(n_order-1-i)
            A[ip][i] = xk_i_n * np.power(x_[ip],i) * np.power(1-x_[ip],n_order-1-i) * C_n1n2

    solution = lstsq(A, b, rcond=None)

    return solution[0].tolist()

def output_foil(x, yu, yl, fname='airfoil.dat', ID=0, info=False):
    '''
    Output airfoil data to tecplot ASCII format file
        x, yu, yl:  current airfoil (list)
        ID:         >0 append to existed file. 0: write header
        info:       True: include curvature, thickness and camber
    '''
    if ID == 0:
        # Write header
        with open(fname, 'w') as f:
            if info: 
                line = 'Variables= X  Y  Curvature Thickness Camber \n '
            else:
                line = 'Variables= X  Y  \n '
            f.write(line)

    if info:
        curv_u, curv_l, thickness, camber = check_foil(x, yu, yl)

    with open(fname, 'a') as f:
        f.write('zone T="Upp-%d" i= %d \n'%(ID, len(x)))
        for i in range(len(x)):
            line = '   %.9f  %.9f'%(x[i], yu[i])
            if info:
                line = line + '  %.9f  %.9f  %.9f'%(curv_u[i], thickness[i], camber[i])
            f.write(line+'\n')
            
        f.write('zone T="Low-%d" i= %d \n'%(ID, len(x)))
        for i in range(len(x)):
            line = '   %.9f  %.9f'%(x[i], yl[i])
            if info:
                line = line + '  %.9f  %.9f  %.9f'%(curv_l[i], thickness[i], camber[i])
            f.write(line+'\n')













