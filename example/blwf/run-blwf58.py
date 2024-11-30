
import os
import time

import numpy as np
from cfdpost.section.physical import PhysicalSec

from cst_modeling.basic import rotate
from cst_modeling.foil import cst_foil
from cst_modeling.surface import Surface
from cst_modeling.tools.blwf import BLWF
from scipy.spatial.distance import cdist

def run_BLWF(Coefs: dict, fname='blwf.in'):
    '''
    Generate BLWF input file with given coefficients
    '''

    _,_,_,airfoil_thick,_ = cst_foil(501, Coefs['cst_u'], Coefs['cst_l'])
    
    wing = Surface(n_sec=2, name='Wing', nn=73, ns=11)
    
    #* Root
    wing.secs[0].xLE   = 31.0   #* Fixed
    wing.secs[0].yLE   = -1.5   #* Fixed
    wing.secs[0].zLE   = 0.0    #* Fixed
    wing.secs[0].chord = 12.0   #* Fixed
    wing.secs[0].twist = Coefs['TwistAngle']
    wing.secs[0].tail  = 0.0    #* Fixed
    wing.secs[0].specified_thickness = airfoil_thick*Coefs['RThickRatio']
    
    wing.secs[0].cst_u = Coefs['cst_u'].copy()
    wing.secs[0].cst_l = Coefs['cst_l'].copy()
    
    half_span = Coefs['AspectRatio']*(1+Coefs['TRRatio'])*wing.secs[0].chord/4

    #* Tip
    wing.secs[1].xLE   = wing.secs[0].xLE + half_span*np.tan(Coefs['SweptAngle']/180.0*np.pi)
    wing.secs[1].yLE   = wing.secs[0].yLE + half_span*np.tan(Coefs['Dihedral']/180.0*np.pi)
    wing.secs[1].zLE   = half_span
    wing.secs[1].chord = Coefs['TRRatio']*wing.secs[0].chord
    wing.secs[1].twist = - Coefs['TwistAngle']
    wing.secs[1].tail  = 0.0    #* Fixed
    wing.secs[1].specified_thickness = wing.secs[0].specified_thickness*Coefs['TRThickRatio']

    wing.secs[1].cst_u = Coefs['cst_u'].copy()
    wing.secs[1].cst_l = Coefs['cst_l'].copy()

    wing.geo()

    #* -------------------------------------------------------
    blwf = BLWF(name='Simple Swept Wing')
    blwf.cst_wing(wing)
    blwf.update_input_wing(fname=fname)
    
    lines = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    f = open(fname, 'w')
    f.write(lines[0])
    f.write(lines[1])
    f.write(' %9.4f %9.4f   0.     32000000. %9.4f %9.4f\n'%(
        Coefs['Minf'], Coefs['AoA'], 8.0, 100.))
    for line in lines[3:]:
        f.write(line)
    f.close()

    os.system('blwf58.exe blwf.in')
    
    return wing


def scale_slice(section: np.ndarray):
    '''
    Rotate and scale sliced section to unit chord length airfoil (Z=const)
    '''

    dis = cdist(section[:,:2], section[0:1,:2], metric='euclidean')
    iLE = np.argmax(dis[:,0])
    xLE = section[iLE,0]
    yLE = section[iLE,1]
    xTE = section[0,0]
    yTE = section[0,1]
    
    section[:,0] -= xLE
    section[:,1] -= yLE
    angle = np.arctan((yTE-yLE)/(xTE-xLE))/np.pi*180.0

    xx, yy, _ = rotate(section[:,0], section[:,1], section[:,2], angle=-angle, axis='Z')
    
    yy = yy/xx[0]
    xx = xx/xx[0]
    cp = section[:,3]
    
    if yy[int(iLE*1.5)]<yy[int(iLE*0.5)]:
        xx = xx[::-1]
        yy = yy[::-1]
        cp = cp[::-1]
    
    return xx, yy, cp

def extract_CL(X, Cp):
    
    CL = 0
    for i in range(X.shape[0]-1):
        CL -= 0.5*(Cp[i+1]+Cp[i])*(X[i+1]-X[i])
    
    return CL


if __name__ == "__main__":

    n_cst  = 10
    
    SpanRatios = [0.25, 0.50, 0.75]     # Fixed slicing locations (from root to tip)

    #* Read CST coefficients
    Coefs = {}
    
    t0 = time.perf_counter()
    with open('input.txt', 'r') as f:
        
        lines = f.readlines()
        n_coef = len(lines)-2*n_cst
        
        for i in range(n_coef):
            line = lines[i].split()
            Coefs[line[0]] = float(line[1])

        cst_u = np.zeros(n_cst)
        cst_l = np.zeros(n_cst)

        for i in range(n_cst):
            line = lines[n_coef+i].split()
            cst_u[i] = float(line[-1])

            line = lines[n_coef+i+n_cst].split()
            cst_l[i] = float(line[-1])

    Coefs['cst_u'] = cst_u
    Coefs['cst_l'] = cst_l

    #* Run BLWF
    wing = run_BLWF(Coefs)
    half_span = wing.secs[-1].zLE
    
    #* Extract slices
    Pref = np.array([20, 0, 0])
    dir_norm = np.array([0, 0, 1])
    locations = np.array(SpanRatios)*half_span
    sections, name_var = BLWF.extract_slice(locations, Pref, dir_norm, dir_ref=np.array([-1.,0.01,0.]), arrange_method='join')
    CLs = []
    with open('Mw.dat', 'w') as f:
        
        f.write('Variables= X Y Cp Mw \n')
        
        for k in range(len(sections)):
            
            nn = sections[k].shape[0]
            f.write('zone T="# %d Loc= %.2f" i= %d \n'%(int(Coefs['ID']), SpanRatios[k], nn))
            
            xx, yy, cp = scale_slice(sections[k])
            mw = PhysicalSec.toMw(cp, Coefs['Minf'])
            
            CLs.append(extract_CL(xx, cp))

            for i in range(nn):
                f.write('% 20.10f  %20.10f  %20.10f  %20.10f \n'%(xx[i], yy[i], cp[i], mw[i]))

    with open('output.txt', 'w') as f:
        for i in range(len(CLs)):
            f.write(' CL-%d  %.9f \n'%(i+1, CLs[i]))

    #* Organize files
    # Do this in run.bat

    t3 = time.perf_counter()
    print('Total time= %.2f(s)'%(t3-t0))


