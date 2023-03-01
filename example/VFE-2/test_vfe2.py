
import numpy as np
from cst_modeling.section import dist_clustcos, intersect_point
from cst_modeling.surface import surf_axisymmetric, output_plot3d, plot3d_to_igs
from numpy.core.numeric import zeros_like

'''
Experimental Surface Pressure Data Obtained on 65° Delta Wing Across Reynolds Number
and Mach Number Ranges (Volume 2—Small-Radius Leading Edges)

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960025648.pdf
'''


def general_eqn(x, x0: float, x1: float, rr: float, l: float, check=True):
    '''
    General equations to define the leading edge semithickness, 
    the flat plate semithickness, the trailing edge closure semithickness,
    and the tranverse radius of the sting fairing.

    >>> phi = general_eqn(x, x0, x1, rr, l, check=True) # phi >= 0

    ### Inputs:
    ```text
    x:      current location
    x0:     min x
    x1:     range of x
    rr:     relative radius
    l:      max phi
    ```

    from Appendix A
    '''
    r = rr*l
    a = np.sqrt(2*r/x1)
    b = -15/8.*a + 3*l/x1
    c = 5/4.*a - 3*l/x1
    d = -3/8.*a + l/x1
    
    xi  = (x-x0)/x1

    if isinstance(xi, float):
        if xi>1.0 or xi<0.0:
            if check:
                raise Exception('Xi >1.0 or < 0.0', xi)
            else:
                xi = max(0.0, min(1.0, xi))
    else:
        if np.any(xi>1.0) or np.any(xi<0.0):
            if check:
                raise Exception('Xi >1.0 or < 0.0', np.max(xi), np.min(xi))
            else:
                for i in range(xi.shape[0]):
                    xi[i] = max(0.0, min(1.0, xi[i]))

    phi = x1*(a*np.sqrt(xi)+b*xi+c*np.power(xi,2)+d*np.power(xi,3))

    return phi

def leading_edge(x, xle: float, rr: float, x1=0.15, l=0.0170008, check=True):
    '''
    xle <= x <= xle+0.15

    x1  = 0.15
    r/l = 0.00, 0.05, 0.15, 0.30
    l   = 0.0170008
    '''
    return general_eqn(x, xle, x1, rr, l, check=check)

def trailing_edge(x, x1=0.1, l=0.0170008, rr=0.0, check=True):
    '''
    The streamwise trailing-edge closure is designed to
    produce a sharp trailing edge and to match the flat plate
    wing at the 90-percent root chord station with continuity
    through the second derivative. 

    0.9 <= x <= 1.0

    x0  = 1
    x1  = 0.10
    r/l = 0
    l   = 0.0170008
    '''
    return general_eqn(1.0-x, 0.0, x1, rr, l, check=check) 

def profile(x: np.ndarray, xle: float, rr: float, x1w=0.15, x1l=0.1, l=0.0170008, tail=0.0):
    '''
    Profile of the delta wing
    
    >>> yy = profile(x, xle, rr, x1w=0.15, x1l=0.1, l=0.0170008)

    0.0 <= xle <= x <= 1.0

    x1w = 0.15   for windward part
    x1l = 0.10   for leeward part
    l   = 0.0170008
    '''
    if xle==1.0:
        return np.zeros_like(x)

    yw = leading_edge(x, xle, rr, x1=x1w, l=l, check=False)
    yl = trailing_edge(x, x1=x1l, l=l, rr=0.0, check=False)

    ii = -1
    for i in range(len(x)-1):
        if yw[i]<=yl[i] and yw[i+1]>yl[i+1]:
            ii = i
            break

    if ii == -1:
        raise Exception('No intersection of LE and TE curves')
        
    yy = np.concatenate((yw[:ii, np.newaxis], yl[ii:, np.newaxis]), axis=0)
    yy = yy[:,0]
    r_ = (x[ii]-x[ii-1])/(x[ii+1]-x[ii-1])
    yy[ii] = (1-r_)*yy[ii-1] + r_*yy[ii+1]
    
    x_ = max(0.9, x[ii])
    for i in range(len(x)):
        if x[i]>=x_:
            r_ = (x[i]-x_)/(1.0-x_)
            yy[i] += r_**2*tail/2.0
    
    return yy

def profile_tip(n: int, i_split: int, xle: float, rr: float, x1w=0.15, x1l=0.1, l=0.0170008, tail=0.0):
    '''
    Profile of the delta wing tip
    
    >>> xx, yy = profile_tip(n, xle, rr, x1w=0.15, x1l=0.1, l=0.0170008, tail=0.0)

    1-x1w-x1l <= xle <= x <= 1.0

    x1w = 0.15   for windward part
    x1l = 0.10   for leeward part
    l   = 0.0170008
    '''
    if xle==1.0:
        return np.ones(n), np.zeros(n)

    #* Locate intersection point
    x_ = dist_clustcos(4*n)
    yw = leading_edge(x_, xle, rr, x1=x1w, l=l, check=False)
    yl = trailing_edge(x_, x1=x1l, l=l, rr=0.0, check=False)
    
    ii = -1
    for i in range(len(x_)-1):
        if yw[i]<=yl[i] and yw[i+1]>yl[i+1]:
            ii = i
            break

    if ii == -1:
        raise Exception('No intersection of LE and TE curves')
    
    pi = intersect_point(np.array([x_[ii],yw[ii]]), np.array([x_[ii+1],yw[ii+1]]),
                        np.array([x_[ii],yl[ii]]), np.array([x_[ii+1],yl[ii+1]]))
    x_split = pi[0]
    y_split = pi[1]
    
    #* Generate actual curve
    xw = dist_clustcos(i_split,     a0=0.010, a1=0.96, beta=2)*(x_split-xle) + xle
    xl = dist_clustcos(n-i_split+1, a0=0.050, a1=0.96)*(1.0-x_split) + x_split
    
    xx = np.concatenate((xw[:, np.newaxis], xl[1:, np.newaxis]), axis=0)
    xx = xx[:,0]
    
    yw = leading_edge(xw, xle, rr, x1=x1w, l=l, check=False)
    yl = trailing_edge(xl, x1=x1l, l=l, rr=0.0, check=False)
    
    for i in range(len(xl)):
        r_ = (xl[i]-x_split)/(1.0-x_split)
        r_ = (1-xle)*r_**2 + xle*r_
        yl[i] += r_*tail/2.0
    
    yy = np.concatenate((yw[:, np.newaxis], yl[1:, np.newaxis]), axis=0)
    yy = yy[:,0]
    
    xx[i_split-1] = x_split
    yy[i_split-1] = y_split

    return xx, yy

def sting_fairing(n=501):
    '''
    The sting is a body of revolution and the sting fairing is designed to 
    emerge from the wing slightly aft of the 60-percent root chord station 
    and to match the constant radius part of the sting slightly ahead of 
    the wing trailing edge. 

    >>> xx, yy = sting_fairing(n=501)
    '''
    x0 =  0.61057
    x1 =  0.36917
    rr =  0.27910261994295
    a  =  0.10040234847327
    b  =  0.33279822819157 
    c  = -0.39554969598736
    d  =  0.13603332984884

    xi = dist_clustcos(n, a1=0.6)
    xx = x0 + xi*x1
    yy = x1*(a*np.sqrt(xi)+b*xi+c*np.power(xi,2)+d*np.power(xi,3))

    return xx, yy

def fore_sting(n=101):
    '''
    The downstream continuation of the sting in the near field of the wing 
    is referred to as the fore-sting. It can be subdivided into the four regions

    >>> xx, yy = fore_sting(n=101)
    '''
    #* region 1: 0.97974 <= x < 1.175
    x1 = np.linspace(0.97974, 1.175, num=n); x1=x1[:-1,np.newaxis]
    y1 = np.ones(n)*0.064119; y1=y1[:-1,np.newaxis]

    #* region 2: 1.175  <= x < 1.253
    x2 = np.linspace(1.175, 1.253, num=n); x2=x2[:-1,np.newaxis]
    y2 = np.linspace(0.064119, 0.06564, num=n); y2=y2[:-1,np.newaxis]

    #* region 3: 1.253  <= x < 1.684
    x3 = np.linspace(1.253, 1.684, num=n); x3=x3[:-1,np.newaxis]
    y3 = np.linspace(0.06564, 0.08258, num=n); y3=y3[:-1,np.newaxis]

    #* region 4: 1.684  <= x < 1.758
    x4 = np.linspace(1.684, 1.758, num=n); x4=x4[:,np.newaxis]
    y4 = np.ones(n)*0.08258; y4=y4[:,np.newaxis]

    xx = np.concatenate((x1, x2, x3, x4), axis=0); xx=xx[:,0]
    yy = np.concatenate((y1, y2, y3, y4), axis=0); yy=yy[:,0]

    return xx, yy

def hind_sting(n=101, rr=0.1):
    '''
    The downstream continuation of the fore-sting.
    Defined by Runze LI.

    >>> xx, yy = hind_sting(n=101, rr=0.1)
    '''
    #* region 4: 1.758 <= x < 1.900
    x0 = 1.758
    x1 = 2.000
    l0 = x1-x0-0.1
    ll = 0.08258
    xx = x0 + (x1-x0)*dist_clustcos(n, a0=0.5)
    yy = general_eqn(x1-xx, 0.0, l0, rr, ll, check=False)

    return xx, yy


if __name__ == '__main__':

    delta = 65
    rr    = 0.05
    tail  = 0.00
    clip_le = 0.0
    clip_tip = 0.00

    span  = 1/np.tan(delta/180.0*np.pi)
    
    nn = 201
    ns = 101

    #* Sting Fairing
    xx, yy = sting_fairing(n=501)
    surf1 = surf_axisymmetric(xx, yy, phi1=180, ns=ns)

    #* Fore Sting
    xx, yy = fore_sting(n=101)
    surf2 = surf_axisymmetric(xx, yy, phi1=180, ns=ns)
    
    #* Hind Sting
    xx, yy = hind_sting(n=301, rr=0.5)
    surf3 = surf_axisymmetric(xx, yy, phi1=180, ns=ns)

    #* Delta Wing
    x_ref = dist_clustcos(nn, a0=0.0079, beta=2)
    
    if True:
        
        surf_x = np.zeros([ns,nn])
        surf_y = np.zeros([ns,nn])
        surf_z = np.zeros([ns,nn])
        tip_x  = np.zeros([ns,nn])
        tip_y  = np.zeros([ns,nn])
        tip_z  = np.zeros([ns,nn])
        
        x0 = 0.75

        z_ref = dist_clustcos(ns, a0=0.1, a1=0.6, beta=2)
        for i in range(ns):
            
            xle = z_ref[i]*x0
            surf_z[i,:] = xle*span
            
            if xle < clip_le:
                r_  = (xle/clip_le)**4
                xle = (1-r_)*0.8*clip_le + r_*clip_le

            surf_x[i,:] = xle + x_ref*(1.0-xle)
            surf_y[i,:] = profile(surf_x[i,:], xle, rr, tail=tail)

        z_ref = dist_clustcos(ns, a0=0.8, a1=0.96)
        for i in range(ns):
            xle = z_ref[i]*(1-x0-clip_tip)+x0
            tip_z[i,:] = xle*span
            tip_x[i,:], tip_y[i,:] = profile_tip(nn, int(nn/2)+1, xle, rr, tail=tail)
        
        if clip_tip>0 or clip_le>0:
            fname='model-clip'
        else:
            fname='model'
        
        output_plot3d(
            [surf_x, surf_x, tip_x, tip_x, surf1[0], surf2[0], surf3[0]], 
            [surf_y,-surf_y, tip_y,-tip_y, surf1[1], surf2[1], surf3[1]], 
            [surf_z, surf_z, tip_z, tip_z, surf1[2], surf2[2], surf3[2]],
            fname+'.grd', scale=1000.0)
        
        plot3d_to_igs(fname=fname)
        
    else:
        
        surf_x = np.zeros([ns,nn])
        surf_y = np.zeros([ns,nn])
        surf_z = np.zeros([ns,nn])
        z_ref  = dist_clustcos(ns, a0=0.5, a1=0.6)
        
        for i in range(ns):

            zle = z_ref[i]*span
            xle = z_ref[i]

            surf_z[i,:] = zle
            surf_x[i,:] = xle + x_ref*(1.0-xle)
            surf_y[i,:] = profile(surf_x[i,:], xle, rr)
        
        output_plot3d(
            [surf_x, surf_x, surf1[0], surf2[0], surf3[0]], 
            [surf_y,-surf_y, surf1[1], surf2[1], surf3[1]], 
            [surf_z, surf_z, surf1[2], surf2[2], surf3[2]], 'model.grd')

    









