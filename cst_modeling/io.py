'''
Input and Output functions for the CST modeling package.
'''

import os
import re
from typing import Tuple, List

import numpy as np


#* ===========================================
#* I/O and format transfer
#* ===========================================

def output_curve(x: np.ndarray, y: np.ndarray, fname='curve.dat', ID=0) -> None:
    '''
    Output airfoil data to tecplot ASCII format file.

    Parameters
    -----------
    x, y : ndarray
        coordinates of the curve.
    ID : int
        if `ID`=0, create new file and write header.
        If `ID`>0, append to existed file.
    '''
    nn = x.shape[0]

    if ID == 0:
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  \n ')

    with open(fname, 'a') as f:
        f.write('zone T="%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            f.write('   %20.9f  %20.9f \n'%(x[i], y[i]))
        f.write('\n')

def output_foil(x: np.ndarray, yu: np.ndarray, yl: np.ndarray, fname='airfoil.dat', ID=0) -> None:
    '''
    Output airfoil data to tecplot ASCII format file

    Parameters
    -----------
    x, yu, yl : ndarray
        coordinates of the baseline airfoil.
        
    ID : int
        if `ID`=0, create new file and write header.
        If `ID`>0, append to existed file.
    '''
    nn = x.shape[0]

    if ID == 0:
        # Write header
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  \n ')

    with open(fname, 'a') as f:
        f.write('zone T="Upp-%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            f.write('   %20.9f  %20.9f \n'%(x[i], yu[i]))
            
        f.write('zone T="Low-%d" i= %d \n'%(ID, nn))
        for i in range(nn):
            f.write('   %20.9f  %20.9f \n'%(x[i], yl[i]))

def output_surface(surf: List[np.ndarray], fname: str, ID=0, zone_name=None) -> None:
    '''
    Output the surface to `*.dat` in Tecplot format.
    
    Parameters
    ------------
    surf : List[np.ndarray]
        surface coordinates, i.e., [surf_x, surf_y, surf_z].
    
    fname : str
        name of the output file.
        
    ID : int
        if `ID`=0, create new file and write header.
        If `ID`>0, append to existed file.
        
    zone_name : str
        name of the zone.
    '''
    # surf_x[ns,nt], ns => spanwise

    ns = surf[0].shape[0]
    nt = surf[0].shape[1]

    if ID == 0:
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

    with open(fname, 'a') as f:

        if zone_name is None:
            zone_name = '%d'%(ID)
        
        f.write('zone T=" %s " i= %d j= %d \n'%(zone_name, nt, ns))

        for i in range(ns):
            for j in range(nt):
                f.write('  %.9f   %.9f   %.9f\n'%(surf[0][i,j], surf[1][i,j], surf[2][i,j]))
    
def read_tecplot(fname='tecplot.dat'):
    '''
    Read a tecplot format data file.
    
    Parameters
    ------------
    fname : str
        file name.
    
    Returns
    -----------
    data : list of ndarray
        data of all zones, shape [ni,nj,nk,nv]. 
    name_var : list of str
        name of variables.
    titles : list of str
        title of zones
    
    Examples
    -------------
    >>> data, name_var, titles = read_tecplot(fname='tecplot.dat')
    

    '''
    name_var = []
    data = []
    titles = []
    n_var = 0
    
    with open(fname, 'r') as f:
        
        lines = f.readlines()
        nLine = len(lines)
        iLine = 0
    
        while iLine < nLine:
            
            line = lines[iLine].split()
            if len(line) == 0:
                iLine += 1
                continue
            
            if line[0] in 'Variables=' or line[0] in 'VARIABLES=' :
                
                line = re.split(r'[=",\s]', lines[iLine])
                while '' in line:
                    line.remove('')

                name_var = line[1:]
                n_var = len(name_var)
                iLine += 1
                continue
        
            if line[0] in 'zone' or line[0] in 'ZONE' or line[0] in 'Zone':
                line = re.split(r'[=\s]', lines[iLine])
                while '' in line:
                    line.remove('')
                
                if 'i' in line:
                    ni = int(line[line.index('i')+1])
                elif 'I' in line:
                    ni = int(line[line.index('I')+1])
                else:
                    ni = 1
                    
                if 'j' in line:
                    nj = int(line[line.index('j')+1])
                elif 'J' in line:
                    nj = int(line[line.index('J')+1])
                else:
                    nj = 1
                    
                if 'k' in line:
                    nk = int(line[line.index('k')+1])
                elif 'K' in line:
                    nk = int(line[line.index('K')+1])
                else:
                    nk = 1
                    
                if 'T' in line:
                    # 非贪婪模式：寻找最短的可能匹配 https://www.cnblogs.com/baxianhua/p/8571967.html
                    str_pat = re.compile(r'\"(.*?)\"')
                    name = str_pat.findall(lines[iLine])
                    titles.append(name[0])
                else:
                    titles.append('')
                    
                data_ = np.zeros((ni,nj,nk,n_var))
                iLine += 1
                
                for k in range(nk):
                    for j in range(nj):
                        for i in range(ni):
                            
                            line = ['#']
                            while line[0] == '#':
                                line = lines[iLine].split()
                                iLine += 1
                                
                            for v in range(n_var):
                                data_[i,j,k,v] = float(line[v])
                                
                data.append(data_.copy())
                continue

    return data, name_var, titles


def read_block_plot3d(lines: list, iLine0: int, ni: int, nj: int, nk: int) -> Tuple[np.ndarray, int]:
    '''
    Read block data from lines.
    
    Parameters
    -----------
    lines : list of str
        f.readlines() of the entire plot3d formate file
        
    iLine0 : int
        the first line of this block is lines[iLine0]
        
    ni, nj, nk: int
        size of this block
        
    Returns
    ---------
    xyz : ndarray
        coordinates, shape `[ni,nj,nk,3]`.
        
    iLine0_new : int
        index of line after read.

    Examples
    ----------
    >>> xyz, iLine0_new = read_block_plot3d(lines, iLine0, ni, nj, nk)

    '''
    xyz = np.zeros([ni,nj,nk,3])
    ll  = iLine0
    ii  = 0
    line = []

    for m in range(3):
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):

                    if ii >= len(line)-1:
                        # Need to read the next line
                        line = lines[ll].split()
                        ii = 0
                        ll += 1
                    else:
                        # Read next value
                        ii += 1

                    xyz[i,j,k,m] = float(line[ii])

    iLine0_new = ll

    return xyz, iLine0_new

def output_plot3d(Xs: List[np.ndarray], Ys: List[np.ndarray], Zs: List[np.ndarray], fname='plot3d.xyz', scale=1.0) -> None:
    '''
    Output surface to fname in plot3d format.
    
    Parameters
    -------------
    Xs, Ys, Zs: List[np.ndarray] [n0][ns,nn]
        coordinates of multiple surfaces
        
    fname: str
        the name of the file (`*.xyz`)
        
    scale: float
        scaling factor.
    '''
    # ns: number of span-wise points
    # nn: number of curve points

    n0 = len(Xs)

    with open(fname, 'w') as f:
        f.write('%d \n '%(n0))     # Number of surfaces
        for i_sec in range(n0):
            ns = Xs[i_sec].shape[0]
            nn = Xs[i_sec].shape[1]
            f.write('%d %d 1\n '%(nn, ns))

        for i_sec in range(n0):
            ii = 0
            ns = Xs[i_sec].shape[0]
            nn = Xs[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(Xs[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

            ii = 0
            ns = Ys[i_sec].shape[0]
            nn = Ys[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(Ys[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

            ii = 0
            ns = Zs[i_sec].shape[0]
            nn = Zs[i_sec].shape[1]
            for i in range(ns):
                for j in range(nn):
                    f.write(' %.9f '%(Zs[i_sec][i,j]*scale))
                    ii += 1
                    if ii%3==0 or (i==ns-1 and j==nn-1):
                        f.write(' \n ')

def output_curves_igs(x: np.ndarray, y: np.ndarray, z: np.ndarray, fname='curve.igs', n_degree=3, is_planar=True):
    '''
    Output curves in the Initial Graphics Exchange Specification (IGES) format.
    
    Parameters
    ------------
    x, y, z : ndarray
        coordinates of the curve(s), [:] or [n_curve,:]
        
    fname : str
        file name. 
        
    n_degree : int
        degree of basis functions.
        0=Determined by data; 1=Line; 2=Circular arc; 
        3=Elliptical arc; 4=Parabolic arc; 5=Hyperbolic arc.
        
    is_planar : bool
        whether is planar curve in X-Y plane.
        
    References
    ------------
    https://wiki.eclipse.org/IGES_file_Specification
    
    https://filemonger.com/specs/igs/devdept.com/version6.pdf
    '''
    
    #* Curve dimension
    if len(x.shape) == 1:
        n_curve = 1
        n_point = x.shape[0]
        x = x[None,:]
        y = y[None,:]
        z = z[None,:]
    else:
        n_curve = x.shape[0]
        n_point = x.shape[1]
    
    #* Output IGES format file
    f = open(fname, 'w')

    #* Start section
    f.write('This is igs file generated by LI Runze. All rights reserved.            S      1\n')
    
    #* Global section
    f.write('1H,,1H;,3Higs,7Higs.igs,44HDASSAULT SYSTEMES CATIA V5 R20 - www.3ds.com,G      1\n')
    f.write('27HCATIA Version 5 Release 20 ,32,75,6,75,15,3Higs,1.0,2,2HMM,1000,1.0, G      2\n')
    f.write('15H20180311.223810,0.001,10000.0,5Hyancy,15HDESKTOP-BEPNROH,11,0,15H2018G      3\n')
    f.write('0311.223810,;                                                           G      4\n')

    #* Data entry section
    iType = 126 # Rational B-Spline Curve
    
    iLineStart = 1
    nLineCount = 3 + 3*n_point + n_degree
    
    for ic in range(n_curve):
        
        # iLineStart: is the line number inside the PD section that has the first line of this entity data.
        # nLineCount: specifies the number of lines this entity takes up in the PD section.

        f.write(' %7d %7d %7d %7d %7d %7d %7d %7d'%(iType, iLineStart, n_degree, 0, 0, 0, 0, 0))
        f.write(' %1d %1d %1d %1dD %6d\n'%(0, 0, 0, 0, ic*2+1))
        f.write(' %7d %7d %7d %7d %7d'%(iType, 0, 0, nLineCount, 0))
        f.write('                BSp Curv{:<3d}'.format(ic*2+1) + '    0D %6d\n'%(ic*2+2))
        
        iLineStart += nLineCount
    
    #* Parameter data section
    iLine = 0
    for ic in range(n_curve):
        
        is_closed = False       # is the starting and ending point the same
        is_polynomial = True    # if all weights are equal, otherwise, the curve is rational
        is_periodic = False     # actually has no difference
        
        # Starting
        iLine += 1
        f.write(' %4d, %4d, %4d, %4d, %4d, %4d, %4d,'%(iType, n_point-1, n_degree, is_planar, is_closed, is_polynomial, is_periodic))
        f.write('%30dP %6d\n'%(ic*2+1, iLine))

        # Knot sequence (n_point + n_degree + 1)
        xKnot = knotx(n_point, n_degree+1)
        for ix in range(xKnot.shape[0]):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(xKnot[ix], ic*2+1, iLine))
        ximin = xKnot[0]
        ximax = xKnot[-1]
        
        # Weight sequence (n_point)
        for _ in range(n_point):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(1.0, ic*2+1, iLine))
        
        # Node coordinates (3*n_point)
        for i in range(n_point):
            iLine += 1
            f.write('%19.10e,%19.10e,%19.10e,%12dP %6d\n'%(
                x[ic,i], y[ic,i], z[ic,i], ic*2+1, iLine))
    
        # Ending
        # Start parameter value, End parameter value, Unit normal x, y, z (if planar)
        # (note: '%**dP' must have at least '%9d')
        iLine += 1
        f.write('%10.3f,%10.3f,%12.5f,%12.5f,%12.5f;%11dP %6d\n'%(
            ximin, ximax, 0, 0, 1, ic*2+1, iLine))
    
    
    #* Ending section
    f.write('S %6dG %6dD %6dP %6d %40s %6d\n'%(1, 4, 2*n_curve, iLine, 'T', 1))
    f.close()


def plot3d_to_igs(fname='igs', plot3d_ext='.xyz', igs_ext='.igs'):
    '''
    Converts Plot3d surface grid file [fname.xyz] to IGES file [fname.igs].
    
    Original Fortran version by Prof. Zhang Yufei: zhangyufei@tsinghua.edu.cn.
    
    Ref: https://wiki.eclipse.org/IGES_file_Specification
    '''

    #* Read plot3d format file
    if not os.path.exists(fname+plot3d_ext):
        raise Exception(fname+' does not exist for format transformation')
    
    with open(fname+plot3d_ext, 'r') as f:
        lines = f.readlines()
        line  = lines[0].split()
        num_block = int(line[0])
        nIJK  = np.zeros([num_block, 5], dtype=int)

        for i in range(num_block):
            line  = lines[i+1].split()
            nIJK[i,0] = int(line[0])
            nIJK[i,1] = int(line[1])
            nIJK[i,2] = int(line[2])
            nIJK[i,3] = idataline(nIJK[i,0], nIJK[i,1])

            if nIJK[i,2]!=1:
                raise Exception('Wrong input file: dimension K is not 1')

            if nIJK[i,0]<4 or nIJK[i,1]<4:
                raise Exception('Wrong input file: dimension I or J less than 4')

        nIJK[0,4] = 1
        for i in range(1, num_block):
            nIJK[i,4] = nIJK[i-1,3] + nIJK[i-1,4]

        kLine = num_block+1

    #* Output IGES format file
    f = open(fname+igs_ext, 'w')

    #* Start section
    f.write('This is igs file generated by ZHANG Yufei. All rights reserved.         S      1\n')
    
    #* Global section
    f.write('1H,,1H;,3Higs,7Higs.igs,44HDASSAULT SYSTEMES CATIA V5 R20 - www.3ds.com,G      1\n')
    f.write('27HCATIA Version 5 Release 20 ,32,75,6,75,15,3Higs,1.0,2,2HMM,1000,1.0, G      2\n')
    f.write('15H20180311.223810,0.001,10000.0,5Hyancy,15HDESKTOP-BEPNROH,11,0,15H2018G      3\n')
    f.write('0311.223810,;                                                           G      4\n')

    #* Data entry section
    iType = 128 # Rational B-Spline Surface
    for ib in range(num_block):
        iLineStart = nIJK[ib, 4]
        iLineEnd   = nIJK[ib, 3]

        f.write(' %7d %7d %7d %7d %7d %7d %7d %7d'%(iType, iLineStart, 0, 0, 0, 0, 0, 0))
        f.write(' %1d %1d %1d %1dD %6d\n'%(0, 0, 0, 0, ib*2+1))
        f.write(' %7d %7d %7d %7d %7d'%(iType, 0, 0, iLineEnd, 0))
        f.write('                BSp Surf{:<3d}'.format(ib+1) + '    0D %6d\n'%(ib*2+2))

    #* Parameter data section
    iLine = 0
    for ib in range(num_block):
        ni = nIJK[ib, 0]
        nj = nIJK[ib, 1]
        nk = nIJK[ib, 2]

        # Starting
        iLine += 1
        f.write(' %4d, %4d, %4d, %4d, %4d,'%(iType, ni-1, nj-1, 3, 3))
        f.write(' %4d, %4d, %4d, %4d, %4d, %11dP %6d\n'%(0, 0, 1, 0, 0, ib*2+1, iLine))

        # Node vector
        xKnot = knotx(ni)
        for ix in range(ni+4):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(xKnot[ix], ib*2+1, iLine))
        ximin = xKnot[0]
        ximax = xKnot[-1]

        xKnot = knotx(nj)
        for ix in range(nj+4):
            iLine += 1
            f.write('%19.10e, %51dP %6d\n'%(xKnot[ix], ib*2+1, iLine))
        xjmin = xKnot[0]
        xjmax = xKnot[-1]

        # Node weight
        for j in range(nj):
            for i in range(ni):
                iLine += 1
                f.write('%19.10e, %51dP %6d\n'%(1.0, ib*2+1, iLine))

        # Node coordinates
        xyz, kLine = read_block_plot3d(lines, kLine, ni, nj, nk)
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    iLine += 1
                    f.write('%19.10e,%19.10e,%19.10e,%12dP %6d\n'%(
                        xyz[i,j,k,0], xyz[i,j,k,1], xyz[i,j,k,2],
                        ib*2+1, iLine))

        # Ending
        iLine += 1
        f.write('%14.6e,%14.6e,%14.6e,%14.6e;%12dP %6d\n'%(
            ximin, ximax, xjmin, xjmax, ib*2+1, iLine))

    #* Ending section
    f.write('S %6dG %6dD %6dP %6d %40s %6d\n'%(1, 3, 2*num_block, iLine, 'T', 1))
    f.close()

def idataline(ni: int, nj: int):
    '''
    Function for `plot3d_to_igs`
    '''
    i1 = ni+4
    i2 = nj+4
    i3 = ni*nj
    i4 = ni*nj
    i5 = 1+1

    return i1+i2+i3+i4+i5

def knotx(ni: int, n_offset=4) -> np.ndarray:
    '''
    Function for `plot3d_to_igs`.
    
    Returns [0, 0, 0, 0, ...(ni-n_offset)..., 1.0, 1.0, 1.0, 1.0].
    ''' 

    xKnot = np.zeros(ni+n_offset)

    for i in range(ni-n_offset+1):
        xKnot[i+n_offset] = (i+1.0)/(ni-n_offset+1)

    for i in range(n_offset):
        xKnot[ni+i] = 1.0

    return xKnot


def output_plot3d_for_parts(fname: str, *args, scale=1.0) -> None:
    '''
    Output surfaces of multiple parts in plot3d format.
    
    Parameters
    ------------
    fname : str
        file name.
        
    *args : List[List[np.ndarray]]
        surfaces of different parts, i.e., part1.surfaces, part2.surfaces, ...
        
        part.surfaces = [[surf_x, surf_y, surf_z], [surf_x, surf_y, surf_z], ...]
    
    scale : float
        scaling factor.
    
    Examples
    ------------
    >>> output_plot3d_for_surfaces(fname='plot3d', surf1, surf2, surf3)

    '''
    n_part = len(args)
    
    Xs = []
    Ys = []
    Zs = []
    
    for i_part in range(n_part):
        
        part = args[i_part]
        
        for i_surf in range(len(part)):
            Xs.append(part[i_surf][0])
            Ys.append(part[i_surf][1])
            Zs.append(part[i_surf][2])
        
    output_plot3d(Xs, Ys, Zs, fname, scale=scale)
