'''
Auxiliary functions for surface modeling
'''
import numpy as np
import copy

from ..math import intersect_index, interp_from_curve
from ..basic import transform
from ..section import Section
from ..surface import Surface


def section_flap(sec: Section, ratio: float, angle: float, dy_axis=None) -> None:
    '''
    Deflect flap by angle (degree) of a section.
    Flap starts from chord-wise ratio = ratio.
    
    Parameters
    ----------
    sec: Section object
        section to be deflected
    
    ratio: float
        chord-wise ratio of the flap rotation axis
        
    angle: float
        deflection angle (degree), +z direction for x-y plane
        
    dy_axis: float, optional
        scaled y location of the rotation axis
    '''
    if angle == 0.0:
        return

    nn = len(sec.xx)
    xu_, xl_, yu_, yl_ = transform(sec.xx, sec.xx, sec.yu, sec.yl, rot=angle, x0=ratio, y0=dy_axis)

    iu1, iu2, _ = intersect_index(sec.xx, sec.yu, xu_, yu_)
    il1, il2, _ = intersect_index(sec.xx, sec.yl, xl_, yl_)
    nu_flap = nn - iu1
    nl_flap = nn - il1

    #* Adjust number of points on the flap
    xu_new2 = np.concatenate((sec.xx[:iu1], xu_[iu2:]), axis=0)
    yu_new2 = np.concatenate((sec.yu[:iu1], yu_[iu2:]), axis=0)
    xl_new2 = np.concatenate((sec.xx[:il1], xl_[il2:]), axis=0)
    yl_new2 = np.concatenate((sec.yl[:il1], yl_[il2:]), axis=0)

    xx_u = np.linspace(sec.xx[iu1], xu_[-1], nu_flap)
    yy_u = interp_from_curve(xx_u, xu_new2, yu_new2)
    xu_new = np.concatenate((sec.xx[:iu1], xx_u), axis=0)
    yu_new = np.concatenate((sec.yu[:iu1], yy_u), axis=0)

    xx_l = np.linspace(sec.xx[il1], xl_[-1], nl_flap)
    yy_l = interp_from_curve(xx_l, xl_new2, yl_new2)
    xl_new = np.concatenate((sec.xx[:il1], xx_l), axis=0)
    yl_new = np.concatenate((sec.yl[:il1], yy_l), axis=0)

    #* Update 3D section
    xu_, xl_, yu_, yl_ = transform(xu_new, xl_new, yu_new, yl_new, 
        scale=sec.chord, rot=sec.twist, dx=sec.xLE, dy=sec.yLE, projection=True)

    sec.x = np.concatenate((np.flip(xl_),xu_[1:]), axis=0)
    sec.y = np.concatenate((np.flip(yl_),yu_[1:]), axis=0)
    sec.z = np.ones(2*nn-1)*sec.zLE


class WingVariableCamber(Surface):
    '''
    Wing with variable camber, sub-class of surface class.

    Parameters
    ----------
    n_sec: int
        number of control sections (2D if set to 0 or 1)
        
    tail: float
        tail thickness (m)
        
    name: str
        name of the surface
        
    fname: str
        name of control file
        
    nn : int
        number of points in the unit 2D curve's `xx`, by default 1001.
        
    ns : int
        number of points in the sweep direction between sections, by default 101.
        
    projection : bool
        whether keeps the projection length the same when rotating the section, by default True.
        
    flap_loc : List[float] [2*n_flap]
        z coordinates of the flap ends, i.e., [z_flap1_1, z_flap1_2, z_flap2_1, z_flap2_2, ...]
    
    flap_trans : float
        transition length of the flap deflection
    
    flap_angle : List[float] [n_flap]
        deflection angle of the flaps

    axis_xloc : List[float] [n_flap]
        chord-wise ratio of the flap rotation axis
        
    axis_dy : List[float] [n_flap]
        scaled y location of the rotation axis

    Notes
    -------
    - +x:     flow direction (m)
    - +y:     upside (m)
    - +z:     span-wise (m)
    - twist:  +z direction (deg)
    - chord:  chord length (m)
    - thick:  relative maximum thickness
    - tail:   absolute tail thickness (m)

    '''
    def __init__(self, n_sec=0, name='WingVC', fname='Wing.txt', **kwargs):

        tail = 0.0
        projection = True
        nn = 1001
        ns = 101

        if 'tail' in kwargs.keys():
            tail = kwargs['tail']

        if 'projection' in kwargs.keys():
            projection = kwargs['projection']

        if 'nn' in kwargs.keys():
            nn = kwargs['nn']

        if 'ns' in kwargs.keys():
            ns = kwargs['ns']

        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, projection=projection)

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

    def build(self, split=True, one_piece=False, f_tecplot='Wing.dat', f_plot3d='Wing.xyz'):
        '''
        Build wing geometry.

        Parameters
        ----------
        split: bool
            generate [surfs] as upper and lower separately
            
        one_piece: bool
            combine the spanwise sections into one piece (for tecplot format)
            
        f_tecplot: str
            file name of tecplot format file. If None, do not output.
            
        f_plot3d: str
            file name of tecplot format file. If None, do not output.
        '''

        z_secs = []
        if self.n_flap > 0:
            for i in range(self.n_flap):
                z_secs.append(self.flap_loc[2*i])
                z_secs.append(self.flap_loc[2*i]+self.flap_trans)
                z_secs.append(self.flap_loc[2*i+1]-self.flap_trans)
                z_secs.append(self.flap_loc[2*i+1])

            self.add_sec(z_secs)

        self.update_sections()

        zLE_secs = self.zLEs

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

    Parameters
    ----------
    z0, z1: float
        span-wise coordinates of the start and end sections
        
    r0, r1: float
        chord-wise ratio defining the deflection axis
        
    trans_len: float
        span-wise transition length (m)


    Attributes
    ------------
    LE0, LE1, TE0, TE1: ndarray
        [x, y, z], LE & TE coordinates of the two end sections of the deflection region
        
    AX0, AX1: ndarray
        [x, y, z], axis ends coordinates. By default, AX = (1-r)*LE + r*TE

    i_sec0, i_sec1: int
        section index to locate z0 and z1

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
        self.AX0 = np.array([x0, y0, self.z0])
        self.AX1 = np.array([x1, y1, self.z1])


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



