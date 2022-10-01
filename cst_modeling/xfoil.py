'''
Run XFoil executable file and read results

'''

import os
import platform

import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import PchipInterpolator


def run_xfoil(AoAs=None, Cls=None,
    Minf=0.1, Re=1e5, nNodes=161, iterVis=10, MReDependence=1,
    fname_airfoil='airfoil.dat', fname_setting='xfoil-input.txt',
    fname_cp='cp.dat', fname_raw='dump.bin', fname_polar='polar.dat',
    fname_log='_xfoil.log', delete_temp=False
    ):
    '''
    
    ### Inputs:
    ```text
    AoAs:       angle of attack (deg):  float, or list [start, end, division] 
    Cls:        lift coefficient:       float, or list [start, end, division] 
    
    Minf:       free stream Mach number
                > 0: take in account compressibility effects through the Prandtl-Glauert correlation (Re must > 0)
                = 0: no compressibility effects
    Re:         Reynolds number (inviscid when Re<=0)
    
    MReDependence:  Re, Mach dependence, (1,2,3)
                    Type   parameters held constant       varying      fixed
                    ----   ------------------------       -------   -----------
                    1    M          , Re            ..   lift     chord, vel.
                    2    M sqrt(CL) , Re sqrt(CL)   ..   vel.     chord, lift
                    3    M          , Re CL         ..   chord    lift , vel.
    
    nNodes:     Number of panel nodes (upper+lower)
    iterVis:    Viscous-solution iteration limit
    
    ```
    '''
    if os.path.exists(fname_polar):
        os.remove(fname_polar)
    
    if fname_raw is not None:
        if os.path.exists(fname_raw):
            os.remove(fname_raw)
        
    #* Write input file
    with open(fname_setting, 'w') as f:
        
        f.write('LOAD %s \n'%(fname_airfoil))       # Airfoil coordinates
        

        f.write('PPAR \n')                          #*  Start paneling menu
        f.write('N %d \n \n'%(nNodes))              # > Number of panel nodes (upper+lower)
        f.write('\n')                               # > Back to main menu
        
        
        f.write('OPER \n')                          #*  Start operation menu
        if Minf > 0:
            f.write('MACH %f \n'%(Minf))            # > Free stream Mach number
        
        if Re > 0:
            f.write('VISC %f \n'%(Re))              # > Activate viscous mode, set Reynolds number
            f.write('ITER %d \n'%(iterVis))         # > Viscous-solution iteration limit
            f.write('TYPE %d \n'%(MReDependence))   # > Re, Mach dependence
        
        f.write('PACC \n')                          # > Activate polar accumulation
        f.write('%s \n'%(fname_polar))              # > File of polar output 
        if fname_raw is not None:
            f.write('%s \n'%(fname_raw))            # > File of dump output
        else:
            f.write('\n')
        
        if AoAs is not None:                        #* Start calculation
            
            if isinstance(AoAs, float) or isinstance(AoAs, int):
                f.write('ALFA %f \n'%(AoAs))            # > Angle of attack
                if fname_cp is not None:
                    f.write('CPWR %s \n'%(fname_cp))    # > Output x, y, Cp to file
            else:
                f.write('ASeq %f %f %f\n'%(AoAs[0], AoAs[1], AoAs[2])) 
            
        elif Cls is not None:
            
            if isinstance(Cls, float) or isinstance(Cls, int):
                f.write('Cl %f \n'%(Cls))               # > Lift coefficient
                if fname_cp is not None:
                    f.write('CPWR %s \n'%(fname_cp))    # > Output x, y, Cp to file
            else:
                f.write('CSeq %f %f %f\n'%(Cls[0], Cls[1], Cls[2]))

            
        else:
            raise Exception('Must specify AoAs or Cls')
        
        
        f.write('PACC \n')                          # > Deactivate polar accumulation
        f.write('\n')                               # > Back to main menu
        f.write('QUIT \n')
        

    #* Run the XFoil calling command
    if platform.system() in 'Windows':
        os.system('xfoil.exe < %s > %s'%(fname_setting, fname_log))
    else:
        if os.path.exists('xfoil.x'):
            os.system('xfoil.x < %s > %s'%(fname_setting, fname_log))
        else:
            os.system('xfoil < %s > %s'%(fname_setting, fname_log))

    if delete_temp:
        os.remove(fname_setting)
        os.remove(fname_log)

def read_xfoil_dump(fname: str, strip=True) -> dict:
    '''
    Read XFoil dump file (Ue,Dstar,Theta,Cf vs s,x,y)
    
    >>> result = read_xfoil_dump(fname)
    
    '''
    f = FortranFile(fname, 'r')
    f.read_ints()                       # NAME, CODE, VERSION
    MACH, REYN, ACRIT           = tuple(f.read_reals(dtype='f4'))
    MATYP, RETYP                = tuple(f.read_ints())
    IITOT, ILETOT, ITETOT, IIB  = tuple(f.read_ints())

    temp = f.read_reals(dtype='f4')     # (XB(IB), YB(IB), IB=1, IIB)
    XB = np.zeros(IIB)
    YB = np.zeros(IIB)
    for i in range(IIB):
        XB[i] = temp[2*i]
        YB[i] = temp[2*i+1]

    LISES = IITOT != 0                  # T if this is an ISES polar, F if XFOIL polar
    LMACH = (MACH <= 0.0) and LISES     # T if this is a Mach sweep, F if alpha sweep


    ALFA=[]; CL=[]; CD=[]; CDI=[]; CM=[]; XTR=[]; MA=[]; II=[]; ILE=[]; ITE=[]; CPSTAR=[]
    X=[]; CP=[]; THET=[]; DSTR=[]; CF=[]; CTAU=[]

    IA = 0
    while True:
        
        try:
            temp = f.read_reals(dtype='f4')
        except:
            break

        X   .append([[],[]])     # shape [N-AoAs, 2, N-point]
        CP  .append([[],[]])
        THET.append([[],[]])
        DSTR.append([[],[]])
        CF  .append([[],[]])
        CTAU.append([[],[]])

        if LISES:   # ISES dump file read

            ALFA.append( temp[0] )
            CL  .append( temp[1] )
            CD  .append( temp[2] )
            CDI .append( temp[3] )
            CM  .append( temp[4] )
            XTR .append([temp[5],temp[6]])
            
            if LMACH:
                
                MA.append( temp[7] )
                
            else:

                if MATYP == 1:
                    MA.append(MACH)
                elif MATYP == 2:
                    MA.append(MACH/np.sqrt(CL[IA]))
                elif MATYP == 3:
                    MA.append(MACH/CL[IA])
                else:
                    raise Exception

            II .append([IITOT, IITOT ])
            ILE.append([ILETOT,ILETOT])
            ITE.append([ITETOT,ITETOT])

        else:       # XFOIL dump file read
            
            ALFA.append( temp[0] )
            CL  .append( temp[1] )
            CD  .append( temp[2] )
            CDI .append( temp[3] )
            CM  .append( temp[4] )
            XTR .append([temp[5],temp[6]])

            if MATYP == 1:
                MA.append(MACH)
            elif MATYP == 2:
                MA.append(MACH/np.sqrt(CL[-1]))
            elif MATYP == 3:
                MA.append(MACH/CL[-1])
            else:
                raise Exception

            temp = f.read_ints()
            II .append([temp[0], temp[1]])
            ITE.append([temp[2], temp[3]])
            ILE.append([1,1])

        for IS in range(2):
            temp = f.read_reals(dtype='f4')
            ii = 0
            for i in range(II[-1][IS]):
                X   [-1][IS].append(temp[ii]);  ii+=1
                CP  [-1][IS].append(temp[ii]);  ii+=1
                THET[-1][IS].append(temp[ii]);  ii+=1
                DSTR[-1][IS].append(temp[ii]);  ii+=1
                CF  [-1][IS].append(temp[ii]);  ii+=1
                CTAU[-1][IS].append(temp[ii]);  ii+=1

        IA += 1
        if IA >= 100000:
            break
    
    f.close()
    NA = len(ALFA)

    GAM = 1.4
    GM1 = GAM - 1.0
    
    for IA in range(NA):
        CPSTAR.append(-1e3)
        if MA[IA] > 0:
            MACHSQ = MA[IA]**2
            CPSTAR[-1] = ( ((1.0+0.5*GM1*MACHSQ)/(1.0+0.5*GM1))**(GAM/GM1)-1.0 )/(0.5*GAM*MACHSQ)

    #* Delete x, cp, etc., that is not the airfoil
    for IA in range(NA):
        for IS in range(2):
            
            n = len(X[IA][IS])
            if strip:
                for i in range(len(X[IA][IS])):
                    if X[IA][IS][i]>1.0:
                        n = i
                        break

            X   [IA][IS] = np.array(X   [IA][IS][:n])
            CP  [IA][IS] = np.array(CP  [IA][IS][:n])
            CF  [IA][IS] = np.array(CF  [IA][IS][:n])
            THET[IA][IS] = np.array(THET[IA][IS][:n])
            DSTR[IA][IS] = np.array(DSTR[IA][IS][:n])
            CTAU[IA][IS] = np.array(CTAU[IA][IS][:n])

    result = {}
    result['numCase']   = NA            # int
    result['Re']        = REYN*1e6      # float
    result['ACRIT  ']   = ACRIT         # float
    result['numCase']   = NA            # int
    result['AoAs']      = np.array(ALFA)    # array [numCase]
    result['CLs']       = np.array(CL)      # array [numCase]
    result['CDs']       = np.array(CD)      # array [numCase]
    result['CDis']      = np.array(CDI)     # array [numCase]
    result['Cms']       = np.array(CM)      # array [numCase]
    result['xTRs']      = np.array(XTR)     # array [numCase]
    result['Minfs']     = np.array(MA)      # array [numCase]
    result['Cp*']       = np.array(CPSTAR)  # array [numCase]
    result['X']         = X                 # list [numCase, 2, nPoint]
    result['Cp']        = CP                # list [numCase, 2, nPoint]
    result['Cf']        = CF                # list [numCase, 2, nPoint]
    result['Theta']     = THET              # list [numCase, 2, nPoint]
    result['DStr']      = DSTR              # list [numCase, 2, nPoint]
    result['CTau']      = CTAU              # list [numCase, 2, nPoint]

    return result

def read_xfoil_polar(fname: str) -> dict:
    '''
    >>> result = read_xfoil_polar(fname)
    '''
    
    f = open(fname, 'r')
    lines = f.readlines()
        
    line = lines[8].split()
    Minf = float(line[2])
    
    alpha=[]; CL=[]; CD=[]; CDp=[]; CM=[]; Top_Xtr=[]; Bot_Xtr=[]

    nn = 0
    for line in lines:
        nn += 1
        line = line.split()
        if len(line)<=1 or nn<=12:
            continue
        
        alpha.append(float(line[0]))
        CL   .append(float(line[1]))
        CD   .append(float(line[2]))
        CDp  .append(float(line[3]))
        CM   .append(float(line[4]))
        Top_Xtr.append(float(line[5]))
        Bot_Xtr.append(float(line[6]))
    
    f.close()
    
    result = {}
    result['numCase']   = len(alpha)        # int
    result['Minfs']     = np.ones_like(alpha)*Minf
    result['AoAs']      = np.array(alpha)   # array [numCase]
    result['CLs']       = np.array(CL)      # array [numCase]
    result['CDs']       = np.array(CD)      # array [numCase]
    result['CDps']      = np.array(CDp)     # array [numCase]
    result['Cms']       = np.array(CM)      # array [numCase]
    result['xTRs-u']    = np.array(Top_Xtr) # array [numCase]
    result['xTRs-l']    = np.array(Bot_Xtr) # array [numCase]

    return result

def foil_for_XFoil(x, yu, yl, fname='airfoil.dat'):
    '''
    Output airfoil geometry for XFoil input
    '''
    with open(fname, 'w') as f:
        
        n = x.shape[0]
        
        f.write('CST-modeling \n')
        for i in range(n):
            f.write(' %9.5f %9.5f \n'%(x[n-i-1], yu[n-i-1]))
        for i in range(n-1):
            f.write(' %9.5f %9.5f \n'%(x[i+1], yl[i+1]))
        f.write('\n')

def xfoil_reconstruction(xx: np.ndarray, yu: np.ndarray, yl: np.ndarray, dump: dict, ii=0):
    '''
    Reconstruct the X, Cp, Cf distribution from Xfoil dump file
    
    Combine to one piece (starts from lower surface trailing edge)
    
    dump is the return from read_xfoil_dump()
    
    ii is the index of the dump results
    
    >>> AoA, Cl, x_, y_, cp_, cf_ = xfoil_reconstruction(xx, yu, yl, dump, ii)
    '''
    nn = xx.shape[0]
    
    if ii >= dump['numCase']:
        print()
        print('Error: index %d > number of cases %d'%(ii, dump['numCase']))
        print()
        raise Exception
    
    AoA = dump['AoAs'][ii]   #type: float
    Cl  = dump['CLs' ][ii]   #type: float
    Xs  = dump['X'   ][ii]
    Cps = dump['Cp'  ][ii]
    Cfs = dump['Cf'  ][ii]
    
    #* Combine to one piece
    _xx = []; _cp = []; _cf = []
    _nn = Xs[1].shape[0]
    for i in range(_nn):
        _xx.append(Xs [1][_nn-i-1])
        _cp.append(Cps[1][_nn-i-1])
        _cf.append(Cfs[1][_nn-i-1])
    for i in range(Xs[0].shape[0]-1):
        _xx.append(Xs [0][i+1])
        _cp.append(Cps[0][i+1])
        _cf.append(Cfs[0][i+1])
    _xx = np.array(_xx)
    _xx = (_xx-np.min(_xx))/(np.max(_xx)-np.min(_xx))
    
    _xu = []; _pu = []; _fu = []
    _xl = []; _pl = []; _fl = []
    
    #* Split by Leading edge
    isLow = True
    for i in range(_xx.shape[0]):
        
        if _xx[i] <= 0.0:
            isLow = False
            _xl = [_xx[i]] + _xl
            _pl = [_cp[i]] + _pl
            _fl = [_cf[i]] + _fl

        if isLow:
            _xl = [_xx[i]] + _xl
            _pl = [_cp[i]] + _pl
            _fl = [_cf[i]] + _fl
        else:
            _xu.append(_xx[i])
            _pu.append(_cp[i])
            _fu.append(_cf[i])

    fpu = PchipInterpolator(_xu, _pu)
    ffu = PchipInterpolator(_xu, _fu)
    fpl = PchipInterpolator(_xl, _pl)
    ffl = PchipInterpolator(_xl, _fl)

    #* Interpolate
    x_ = np.zeros(2*nn-1)
    y_ = np.zeros(2*nn-1)
    cp_= np.zeros(2*nn-1)
    cf_= np.zeros(2*nn-1)
    
    ii = 0
    for i in range(nn):
        x_ [ii] = xx[nn-i-1]
        y_ [ii] = yl[nn-i-1]
        cp_[ii] = fpl(x_[ii])
        cf_[ii] = ffl(x_[ii])
        ii += 1
        
    for i in range(nn-1):
        x_ [ii] = xx[i+1]
        y_ [ii] = yu[i+1]
        cp_[ii] = fpu(x_[ii])
        cf_[ii] = ffu(x_[ii])
        ii += 1
    
    return AoA, Cl, x_, y_, cp_, cf_


