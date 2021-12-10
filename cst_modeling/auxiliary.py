'''
Auxiliary functions for surface modeling
'''
import numpy as np
import copy

import basic
import foil


from foil import Section, transform, interplot_from_curve, curve_intersect
from surface import Surface

def section_flap(sec: Section, ratio, angle, dy_axis=None):
    '''
    Deflect flap by angle (degree) of a section.
    Flap starts from chord-wise ratio = ratio.

    ### Inputs:
    ```text
    ratio:      chord-wise ratio of the flap rotation axis
    angle:      deflection angle (degree), +z direction for x-y plane
    dy_axis:    scaled y location of the rotation axis
    ```
    '''
    if angle == 0.0:
        return

    nn = len(sec.xx)
    xu_, xl_, yu_, yl_ = transform(sec.xx, sec.xx, sec.yu, sec.yl, rot=angle, x0=ratio, y0=dy_axis)

    iu1, iu2, _ = curve_intersect(sec.xx, sec.yu, xu_, yu_)
    il1, il2, _ = curve_intersect(sec.xx, sec.yl, xl_, yl_)
    nu_flap = nn - iu1
    nl_flap = nn - il1

    #* Adjust number of points on the flap
    xu_new2 = np.concatenate((sec.xx[:iu1], xu_[iu2:]), axis=0)
    yu_new2 = np.concatenate((sec.yu[:iu1], yu_[iu2:]), axis=0)
    xl_new2 = np.concatenate((sec.xx[:il1], xl_[il2:]), axis=0)
    yl_new2 = np.concatenate((sec.yl[:il1], yl_[il2:]), axis=0)

    xx_u = np.linspace(sec.xx[iu1], xu_[-1], nu_flap)
    yy_u = interplot_from_curve(xx_u, xu_new2, yu_new2)
    xu_new = np.concatenate((sec.xx[:iu1], xx_u), axis=0)
    yu_new = np.concatenate((sec.yu[:iu1], yy_u), axis=0)

    xx_l = np.linspace(sec.xx[il1], xl_[-1], nl_flap)
    yy_l = interplot_from_curve(xx_l, xl_new2, yl_new2)
    xl_new = np.concatenate((sec.xx[:il1], xx_l), axis=0)
    yl_new = np.concatenate((sec.yl[:il1], yy_l), axis=0)

    #* Update 3D section
    xu_, xl_, yu_, yl_ = transform(xu_new, xl_new, yu_new, yl_new, 
        scale=sec.chord, rot=sec.twist, dx=sec.xLE, dy=sec.yLE, proj=True)

    sec.x = np.concatenate((np.flip(xl_),xu_[1:]), axis=0)
    sec.y = np.concatenate((np.flip(yl_),yu_[1:]), axis=0)
    sec.z = np.ones(2*nn-1)*sec.zLE


class WingVariableCamber(Surface):
    '''
    Wing with variable camber, sub-class of surface class.

    ### Inputs:
    ```text
    n_sec:   number of control sections (2D if set to 0 or 1)
    tail:    tail thickness (m)
    name:    name of the surface
    fname:   name of control file
    nn:      number of points of upper/lower section
    ns:      number of spanwise
    project: True ~ projected chord length does not change when twisted

    flap_loc:   list [2*n_flap], z coordinates of the flap ends. 
                [z_flap1_1, z_flap1_2, z_flap2_1, z_flap2_2, ...]
    flap_trans: float, transition length of the flap deflection
    flap_angle: list [n_flap], deflect angle of flaps

    axis_xloc:  list [n_flap], chord-wise ratio of the flap rotation axis
    axis_dy:    empty list or list [n_flap], scaled y location of the rotation axis
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
    '''
    def __init__(self, n_sec=0, name='WingVC', fname='Wing.txt', **kwargs):

        tail = 0.0
        project = True
        nn = 1001
        ns = 101

        if 'tail' in kwargs.keys():
            tail = kwargs['tail']

        if 'project' in kwargs.keys():
            project = kwargs['project']

        if 'nn' in kwargs.keys():
            nn = kwargs['nn']

        if 'ns' in kwargs.keys():
            ns = kwargs['ns']

        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, project=project)

        self.read_setting(fname, tail=tail)

        self.flap_trans = 0.2
        self.flap_loc = []
        self.flap_angle = []

        self.axis_xloc = []
        self.axis_dy = []

        if 'flap_trans' in kwargs.keys():
            self.flap_trans = kwargs['flap_trans']

        if 'flap_loc' in kwargs.keys():
            self.flap_loc = kwargs['flap_loc']

        if 'flap_angle' in kwargs.keys():
            self.flap_angle = kwargs['flap_angle']

        if 'axis_xloc' in kwargs.keys():
            self.axis_xloc = kwargs['axis_xloc']

        if 'axis_dy' in kwargs.keys():
            self.axis_dy = kwargs['axis_dy']

        self.n_flap = len(self.flap_angle)
        if len(self.flap_loc) != 2*self.n_flap:
            raise Exception('Size of flap_loc %d does not match 2*flap_angle %d'%(len(self.flap_loc), 2*self.n_flap))

        if len(self.axis_xloc) != self.n_flap:
            raise Exception('Size of axis_xloc %d does not match flap_angle %d'%(len(self.axis_xloc), self.n_flap))

    def build(self, split=True, one_piece=False, f_tecplot='Wing.dat', f_plot3d='Wing.grd'):
        '''
        Build wing geometry.

        ### Inputs:
        ```text
        split:      True ~ generate [surfs] as upper and lower separately
        one_piece:  True ~ combine the spanwise sections into one piece (for tecplot format)

        f_tecplot:  file name of tecplot format file. If None, do not output.
        f_plot3d:   file name of tecplot format file. If None, do not output.
        ```
        '''

        z_secs = []
        if self.n_flap > 0:
            for i in range(self.n_flap):
                z_secs.append(self.flap_loc[2*i])
                z_secs.append(self.flap_loc[2*i]+self.flap_trans)
                z_secs.append(self.flap_loc[2*i+1]-self.flap_trans)
                z_secs.append(self.flap_loc[2*i+1])

            self.add_sec(z_secs)

        self.geo_secs()

        zLE_secs = self.zLE_secs

        for i in range(self.n_flap):
            i1 = zLE_secs.index(z_secs[4*i]) + 1
            i2 = zLE_secs.index(z_secs[4*i+1]) + 1
            ax = self.axis_xloc[i]
            aa = self.flap_angle[i]

            if len(self.axis_dy) == self.n_flap:
                ay = self.axis_dy[i]
            else:
                ay = None

            section_flap(self.secs[i1], ax, aa, dy_axis=ay)
            section_flap(self.secs[i2], ax, aa, dy_axis=ay)

        self.geo(update_sec=False)

        if f_tecplot is not None:
            self.output_tecplot(fname=f_tecplot, one_piece=one_piece, split=split)

        if f_plot3d is not None:
            self.output_plot3d(fname=f_plot3d)

class DeflectSurf():
    '''
    Deflecting surface by axis defined by (x0, y0, z0) and (z1, y1, z1).

    Movable region: chord-wise ratio from r0 to 1, span-wise from z0 to z1

    >>> DeflectSurf(surf: Surface, z0, z1, r0, r1, trans_len=0.5)

    ### Inputs:
    ```text
    z0, z1:     span-wise coordinates of the start and end sections
    r0, r1:     chord-wise ratio defining the deflection axis
    trans_len:  span-wise transition length (m)
    ```

    ### Parameters:
    ```text
    LE0, LE1, TE0, TE1: ndarray, [x, y, z], LE & TE coordinates of the two end sections of the deflection region
    AX0, AX1: ndarray, [x, y, z], axis ends coordinates. By default, AX = (1-r)*LE + r*TE

    i_sec0, i_sec1: section index to locate z0 and z1
    ```
    '''

    def __init__(self, surf: Surface, z0, z1, r0, r1, trans_len=0.5):
        
        self.surf = surf
        self.z0 = z0
        self.z1 = z1
        self.r0 = r0
        self.r1 = r1
        self.trans_len = trans_len

        self._locate_ends()

    def set_axis(self, x0, y0, x1, y1):
        '''
        Define the rotation axis.

        AX0, AX1: ndarray, [x, y, z], axis ends coordinates \n

        If by default, (x, y) = (1-r)*LE + r*TE.
        '''
        self.AX0 = np.array([self.x0, self.y0, self.z0])
        self.AX1 = np.array([self.x1, self.y1, self.z1])

    def deflect(self, angle: float):
        '''
        Deflecting surface

        Inputs:
        ---
        angle:  deflection angle (degree), positive means downwards deflection
        '''

        surf2 = copy.deepcopy(self.surfs)




        pass

    def _locate_ends(self):
        '''
        Locate two end sections of the deflection region. \n
        
        LE0, LE1, TE0, TE1, AX0, AX1 (ndarray, [x, y, z])
        '''

        secs = self.surf.secs

        for j in range(self.surf.n_sec-1):

            if (secs[j].zLE-self.z0)*(secs[j+1].zLE-self.z0)<0.0:
                s0 = (self.z0 - secs[j].zLE)/(secs[j+1].zLE-secs[j].zLE)

                X0 = np.array([secs[j  ].xLE, secs[j  ].yLE, secs[j  ].zLE])
                X1 = np.array([secs[j+1].xLE, secs[j+1].yLE, secs[j+1].zLE])
                self.LE0 = (1-s0)*X0 + s0*X1

                X0 = np.array([secs[j  ].x[0]+secs[j  ].x[-1], 
                               secs[j  ].y[0]+secs[j  ].y[-1], 
                               secs[j  ].z[0]+secs[j  ].z[-1]])*0.5
                X1 = np.array([secs[j+1].x[0]+secs[j+1].x[-1], 
                               secs[j+1].y[0]+secs[j+1].y[-1], 
                               secs[j+1].z[0]+secs[j+1].z[-1]])*0.5
                self.TE0 = (1-s0)*X0 + s0*X1
                self.AX0 = (1-self.r0)*self.LE0 + self.r0*self.TE0

                self.i_sec0 = j
                
            if (secs[j].zLE-self.z1)*(secs[j+1].zLE-self.z1)<0.0:
                s1 = (self.z1 - secs[j].zLE)/(secs[j+1].zLE-secs[j].zLE)
                i1 = j

                X0 = np.array([secs[j  ].xLE, secs[j  ].yLE, secs[j  ].zLE])
                X1 = np.array([secs[j+1].xLE, secs[j+1].yLE, secs[j+1].zLE])
                self.LE1 = (1-s1)*X0 + s1*X1

                X0 = np.array([secs[j  ].x[0]+secs[j  ].x[-1], 
                               secs[j  ].y[0]+secs[j  ].y[-1], 
                               secs[j  ].z[0]+secs[j  ].z[-1]])*0.5
                X1 = np.array([secs[j+1].x[0]+secs[j+1].x[-1], 
                               secs[j+1].y[0]+secs[j+1].y[-1], 
                               secs[j+1].z[0]+secs[j+1].z[-1]])*0.5
                self.TE1 = (1-s1)*X0 + s1*X1
                self.AX1 = (1-self.r1)*self.LE1 + self.r1*self.TE1

                self.i_sec1 = j



#! ===========================================
#! No use
#! ===========================================

def cst_foil_roundtail(nn, cst_u, cst_l, cst_tu, cst_tl, n_rev=10, x=None, t=None):
    '''
    Constructing upper and lower curves of an airfoil based on CST method

    CST:    class shape transfermation method (Kulfan, 2008)

    >>> x_, yu, yl, t0, R0, (yu_rev, yl_rev) = 
    >>>                 cst_foil_roundtail(nn, cst_u, cst_l, 
    >>>                 cst_tu, cst_tl, n_rev=10, x=None, t=None)

    ### Inputs:
    ```text
    nn:             total amount of points
    cst_u,  cst_l:  CST coefficients of upper/lower surface (ndarray)
    cst_tu, cst_tl: the (1st) coefficient of the CST curve describing the round tail
    n_rev:          order of CST for round tail
    x:      point x [0,1] (optional ndarray, size is nn)
    t:      relative maximum thickness (optional)
    ```

    ### Return
    ```text
    x (ndarray), y_upp (ndarray), y_low (ndarray), t0, R0
    ```
    '''
    x_, yu_ = foil.cst_curve(nn, cst_u, x=x)
    x_, yl_ = foil.cst_curve(nn, cst_l, x=x)

    rev_u = np.zeros(n_rev)
    rev_l = np.zeros(n_rev)
    if isinstance(cst_tu, float):
        rev_u[0] = cst_tu
    elif isinstance(cst_tu, np.ndarray):
        rev_u[:cst_tu.shape[0]] = cst_tu
    else:
        raise Exception('Wrong type of cst_tu')

    if isinstance(cst_tl, float):
        rev_l[0] = cst_tl
    elif isinstance(cst_tl, np.ndarray):
        rev_l[:cst_tl.shape[0]] = cst_tl
    else:
        raise Exception('Wrong type of cst_tl')

    _,  yu_rev = foil.cst_curve(nn, rev_u, x=1.0-x_)
    _,  yl_rev = foil.cst_curve(nn, rev_l, x=1.0-x_)

    yu = yu_ + yu_rev
    yl = yl_ + yl_rev

    thick = yu-yl
    t0 = np.max(thick)

    #* Apply thickness constraint
    if t is not None:
        r  = t/t0
        t0 = t
        yu = yu * r
        yl = yl * r

    # Calculate leading edge radius
    x_RLE = 0.005
    yu_RLE = basic.interplot_from_curve(x_RLE, x_, yu)
    yl_RLE = basic.interplot_from_curve(x_RLE, x_, yl)
    R0, _ = foil.find_circle_3p([0.0,0.0], [x_RLE,yu_RLE], [x_RLE,yl_RLE])

    return x_, yu, yl, t0, R0, (yu_rev, yl_rev)

def fit_curve_roundtail(x, y, n_cst=7, n_rev=10, n_tail=1, ratio_tail=0.1, xn1=0.5, xn2=1.0):
    '''
    Using CST method to fit a curve with a round tail

    This function allows the curve chord length not equals to one.
    Also allows the curve has a twist angle.

    >>> cst, coef_tail = fit_curve_roundtail(x, y, 
    >>>         n_cst=7, n_rev=10, n_tail=1, ratio_tail=0.1, xn1=0.5, xn2=1.0)

    ### Inputs:
    ```text
    x, y:       curve points, ndarray
    n_cst:      number of CST parameters
    n_rev:      size of cst_rev (ndarray), which is the CST parameters for 
                the reverse curve (y_tail) that is used to describe the round tail
    n_tail:     the first n_tail CST coefficients in cst_rev are not zero 
    ratio_tail: the reverse curve (y_tail) is used to mainly describe (1-ratio_tail~1)
                part of the entire curve
    ```

    ### Return: 
    ```text
    cst_u, cst_l:   ndarray [n_cst], the CST coefficients describing the remaining curve,
                    i.e., the original y substracts the round tail curve (y_tail)
    coef_tail:      ndarray [n_tail], the first n_tail CST coefficients
                    for the reverse curve (y_tail) that is used to describe the round tail
    ```
    '''
    n0 = x.shape[0]
    n_tail = max(1, min(n_tail, n_rev))

    #* Curve with twist and chord not 1.0
    chord = np.sqrt((x[0]-x[-1])**2+(y[0]-y[-1])**2)
    twist = np.arctan((y[-1]-y[0])/(x[-1]-x[0]))*180/np.pi

    x_ = (x - x[0])/chord
    y_ = (y - y[0])/chord
    x_, y_, _ = basic.rotate(x_, y_, None, angle=-twist, axis='Z')

    #* Fit the trailing part
    ip1 = max(2, int(ratio_tail*n0))
    coef_tail = foil.fit_curve_partial(x_, np.flip(y_), ip0=0, ip1=ip1,
                    n_cst=n_rev, ic0=0, ic1=n_tail, xn1=xn1, xn2=xn2)

    cst_rev = np.zeros(n_rev)
    cst_rev[:n_tail] = coef_tail
    _, y_tail = foil.cst_curve(n0, cst_rev, x=1.0-x_, xn1=xn1, xn2=xn2)
    
    #* Fit the remaining curve
    coef = foil.fit_curve(x_, y_-y_tail, n_cst=n_cst, xn1=xn1, xn2=xn2)

    return coef, coef_tail

def cst_foil_roundend(nn, cst_u, cst_l, cst_lu, cst_ll, n_lead_all,
                        cst_tu, cst_tl, n_tail_all, x=None, t=None):
    '''
    Constructing upper and lower curves of an airfoil based on CST method

    CST:    class shape transfermation method (Kulfan, 2008)

    >>> x_, yu, yl, t0, R0, (yu_lead, yl_lead), (yu_tail, yl_tail) = 
    >>>     cst_foil_roundend(nn, cst_u, cst_l, cst_lu, cst_ll, n_lead_all,
    >>>                         cst_tu, cst_tl, n_tail_all, x=None, t=None)

    ### Inputs:
    ```text
    nn:             total amount of points
    cst_u,  cst_l:  CST coefficients of upper/lower surface (ndarray)
    cst_lu, cst_ll: the (1st) coefficient of the CST curve describing the round leading edge
    cst_tu, cst_tl: the (1st) coefficient of the CST curve describing the round tail
    n_lead_all:     order of CST for round leading edge
    n_tail_all:     order of CST for round tail
    x:      point x [0,1] (optional ndarray, size is nn)
    t:      relative maximum thickness (optional)
    ```

    ### Return
    x (ndarray), y_upp (ndarray), y_low (ndarray), t0, R0
    '''
    #* Baseline curve
    x_, yu_ = foil.cst_curve(nn, cst_u, x=x)
    x_, yl_ = foil.cst_curve(nn, cst_l, x=x)

    #* Round leading edge
    lead_u = np.zeros(n_lead_all)
    lead_l = np.zeros(n_lead_all)
    if isinstance(cst_lu, float):
        lead_u[0] = cst_lu
    elif isinstance(cst_lu, np.ndarray):
        lead_u[:cst_lu.shape[0]] = cst_lu
    else:
        raise Exception('Wrong type of cst_lu')

    if isinstance(cst_ll, float):
        lead_l[0] = cst_ll
    elif isinstance(cst_ll, np.ndarray):
        lead_l[:cst_ll.shape[0]] = cst_ll
    else:
        raise Exception('Wrong type of cst_ll')

    _,  yu_lead = foil.cst_curve(nn, lead_u, x=x_)
    _,  yl_lead = foil.cst_curve(nn, lead_l, x=x_)

    #* Round tail
    rev_u = np.zeros(n_tail_all)
    rev_l = np.zeros(n_tail_all)
    if isinstance(cst_tu, float):
        rev_u[0] = cst_tu
    elif isinstance(cst_tu, np.ndarray):
        rev_u[:cst_tu.shape[0]] = cst_tu
    else:
        raise Exception('Wrong type of cst_tu')

    if isinstance(cst_tl, float):
        rev_l[0] = cst_tl
    elif isinstance(cst_tl, np.ndarray):
        rev_l[:cst_tl.shape[0]] = cst_tl
    else:
        raise Exception('Wrong type of cst_tl')

    _,  yu_tail = foil.cst_curve(nn, rev_u, x=1.0-x_)
    _,  yl_tail = foil.cst_curve(nn, rev_l, x=1.0-x_)

    #* Combine
    yu = yu_ + yu_lead + yu_tail
    yl = yl_ + yl_lead + yl_tail

    thick = yu-yl
    t0 = np.max(thick)

    #* Apply thickness constraint
    if t is not None:
        r  = t/t0
        t0 = t
        yu = yu * r
        yl = yl * r

    # Calculate leading edge radius
    x_RLE = 0.005
    yu_RLE = basic.interplot_from_curve(x_RLE, x_, yu)
    yl_RLE = basic.interplot_from_curve(x_RLE, x_, yl)
    R0, _ = foil.find_circle_3p([0.0,0.0], [x_RLE,yu_RLE], [x_RLE,yl_RLE])

    return x_, yu, yl, t0, R0, (yu_lead, yl_lead), (yu_tail, yl_tail)

def fit_curve_roundend(x, y, n_cst=7, n_lead_all=10, n_lead=1, ratio_lead=0.1, 
            n_tail_all=10, n_tail=1, ratio_tail=0.1, xn1=0.5, xn2=1.0):
    '''
    Using CST method to fit a curve with round ends

    This function allows the curve chord length not equals to one.
    Also allows the curve has a twist angle.

    >>> cst_u, cst_l, coef_tail = fit_curve_roundtail(x, y, 
    >>>         n_cst=7, n_rev=10, n_tail=1, ratio_tail=0.1, xn1=0.5, xn2=1.0)

    ### Inputs:
    ```text
    x, y:       curve points, ndarray
    n_cst:      number of CST parameters
    n_tail_all: size of cst_rev (ndarray), which is the CST parameters for 
                the reverse curve (y_tail) that is used to describe the round tail
    n_tail:     the first n_tail CST coefficients in cst_rev are not zero 
    ratio_tail: the reverse curve (y_tail) is used to mainly describe (1-ratio_tail~1)
                part of the entire curve
    ```

    ### Return: 
    ```text
    cst_u, cst_l:   ndarray [n_cst], the CST coefficients describing the remaining curve,
                    i.e., the original y substracts the round tail curve (y_tail)
    coef_tail:      ndarray [n_tail], the first n_tail CST coefficients
                    for the reverse curve (y_tail) that is used to describe the round tail
    ```
    '''
    n0 = x.shape[0]
    n_lead = max(1, min(n_lead, n_lead_all))
    n_tail = max(1, min(n_tail, n_tail_all))

    #* Curve with twist and chord not 1.0
    chord = np.sqrt((x[0]-x[-1])**2+(y[0]-y[-1])**2)
    twist = np.arctan((y[-1]-y[0])/(x[-1]-x[0]))*180/np.pi

    x_ = (x - x[0])/chord
    y_ = (y - y[0])/chord
    x_, y_, _ = basic.rotate(x_, y_, None, angle=-twist, axis='Z')

    #* Fit the leading part
    ip1 = max(2, int(ratio_lead*n0))
    coef_lead = foil.fit_curve_partial(x_, y_, ip0=0, ip1=ip1,
                    n_cst=n_lead_all, ic0=0, ic1=n_lead, xn1=xn1, xn2=xn2)

    cst_lead = np.zeros(n_lead_all)
    cst_lead[:n_lead] = coef_lead
    _, y_lead = foil.cst_curve(n0, cst_lead, x=x_, xn1=xn1, xn2=xn2)

    #* Fit the trailing part
    ip1 = max(2, int(ratio_tail*n0))
    coef_tail = foil.fit_curve_partial(x_, np.flip(y_), ip0=0, ip1=ip1,
                    n_cst=n_tail_all, ic0=0, ic1=n_tail, xn1=xn1, xn2=xn2)

    cst_rev = np.zeros(n_tail_all)
    cst_rev[:n_tail] = coef_tail
    _, y_tail = foil.cst_curve(n0, cst_rev, x=1.0-x_, xn1=xn1, xn2=xn2)
    
    #* Fit the remaining curve
    coef = foil.fit_curve(x_, y_-y_lead-y_tail, n_cst=n_cst, xn1=xn1, xn2=xn2)

    return coef, coef_lead, coef_tail




