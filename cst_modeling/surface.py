'''
Construct surface with sections/open-sections
'''
import os
from typing import List

import numpy as np

from .basic import BasicSurface, rotate
from .section import OpenSection, Section

#* ===========================================
#* CST surfaces
#* ===========================================

class OpenSurface(BasicSurface):
    '''
    Open surface defined by multiple OpenSection objects.
    
    Examples
    ----------
    >>> surf = OpenSurface(n_sec=0, name='Patch', nn=1001, ns=101, projection=True)
    '''
    def __init__(self, n_sec=0, name='Patch', nn=1001, ns=101, projection=True):

        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, projection=projection)

        n_ = max(1, n_sec)
        self.secs = [ OpenSection() for _ in range(n_) ]

    def read_setting(self, fname) -> None:
        '''
        Read in Surface layout and CST parameters from file.

        Parameters
        ----------
        fname : str
            settings file name

        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3}

        found_surf = False
        found_key = 0
        with open(fname, 'r') as f:

            lines = f.readlines()
            iL = 0

            while iL<len(lines):

                line = lines[iL].split()

                if len(line) < 1:
                    iL += 1
                    continue
                
                if not found_surf and len(line) > 1:
                    if '[Surf]' in line[0] and self.name == line[1]:
                        found_surf = True

                elif found_surf and '[Surf]' in line[0]:
                    break

                elif found_surf and found_key == 0:
                    if line[0] in key_dict:
                        found_key = key_dict[line[0]]

                elif found_surf and found_key == 1:
                    for i in range(self.n_sec):
                        iL += 1
                        line = lines[iL].split()
                        self.secs[i].xLE   = float(line[0])
                        self.secs[i].yLE   = float(line[1])
                        self.secs[i].zLE   = float(line[2])
                        self.secs[i].chord = float(line[3])
                        self.secs[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.secs[i].thick_set = float(line[5])

                        if self.l2d:
                            self.secs[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_sec):
                        iL += 2
                        line = lines[iL].split()
                        self.secs[i].cst = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()
                        cst_r = np.zeros(n_cst_refine)

                        i1 = 0

                        for j in range(n_cst_refine):
                            if j>=i_cst_start-1 and i1<len(line1):
                                cst_r[j] = float(line1[i1])
                                i1 += 1

                        self.secs[i].refine = cst_r

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1
        
        print('Read surface [%s] settings'%(self.name))

        self.layout_center()


class Surface(BasicSurface):
    '''
    Surface defined by multiple Section objects, i.e., foils

    Examples
    -----------
    >>> surf = Surface(n_sec=0, name='Wing', nn=1001, ns=101, projection=True)
    '''
    def __init__(self, n_sec=0, name='Wing', nn=1001, ns=101, projection=True):

        super().__init__(n_sec=n_sec, name=name, nn=nn, ns=ns, projection=projection)

        self.secs = [ Section() for _ in range(max(1, n_sec)) ]

    def read_setting(self, fname, tail=0.0) -> None:
        '''
        Read in Surface layout and CST parameters from file

        Parameters
        ----------
        fname : str
            settings file name.
        tail : float or list
            tail thickness (m) of each section.
        '''
        if not os.path.exists(fname):
            raise Exception(fname+' does not exist for surface read setting')
        
        key_dict = {'Layout:': 1, 'CST_coefs:': 2, 'CST_refine:': 3, 'CST_flip:': 4}

        found_surf = False
        found_key = 0
        with open(fname, 'r') as f:

            lines = f.readlines()
            iL = 0

            while iL<len(lines):

                line = lines[iL].split()

                if len(line) < 1:
                    iL += 1
                    continue
                
                if not found_surf and len(line) > 1:
                    if '[Surf]' in line[0] and self.name == line[1]:
                        found_surf = True

                elif found_surf and '[Surf]' in line[0]:
                    break

                elif found_surf and found_key == 0:
                    if line[0] in key_dict:
                        found_key = key_dict[line[0]]

                elif found_surf and found_key == 1:
                    for i in range(self.n_sec):
                        iL += 1
                        line = lines[iL].split()
                        self.secs[i].xLE   = float(line[0])
                        self.secs[i].yLE   = float(line[1])
                        self.secs[i].zLE   = float(line[2])
                        self.secs[i].chord = float(line[3])
                        self.secs[i].twist = float(line[4])

                        if len(line) >= 6:
                            self.secs[i].thick_set = float(line[5])

                        if isinstance(tail, float):
                            self.secs[i].tail  = tail/self.secs[i].chord
                        elif len(tail)==self.n_sec:
                            self.secs[i].tail  = tail[i]/self.secs[i].chord
                        else:
                            raise Exception('tail must be a float or a list with length = section number')
                        
                        if self.secs[i].thick_set <= 0.0:
                            self.secs[i].thick_set = None

                        if self.l2d:
                            self.secs[i].zLE = 0.0

                    found_key = 0

                elif found_surf and found_key == 2:
                    for i in range(self.n_sec):
                        iL += 2
                        line = lines[iL].split()
                        self.secs[i].cst_u = np.array([float(aa) for aa in line])

                        iL += 1
                        line = lines[iL].split()
                        self.secs[i].cst_l = np.array([float(aa) for aa in line])
                    
                    found_key = 0

                elif found_surf and found_key == 3:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])
                    i_cst_start = int(line[1])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()

                        iL += 1
                        line2 = lines[iL].split()

                        cst_ur = np.zeros(n_cst_refine)
                        cst_lr = np.zeros(n_cst_refine)

                        i1 = 0
                        i2 = 0
                        for j in range(n_cst_refine):
                            if j>=i_cst_start-1 and i1<len(line1):
                                cst_ur[j] = float(line1[i1])
                                i1 += 1
                            if j>=i_cst_start-1 and i2<len(line2):
                                cst_lr[j] = float(line2[i2])
                                i2 += 1

                        self.secs[i].refine_fixed_t=True
                        self.secs[i].refine_u = cst_ur
                        self.secs[i].refine_l = cst_lr

                    found_key = 0

                elif found_surf and found_key == 4:
                    iL += 2
                    line = lines[iL].split()
                    n_cst_refine = int(line[0])

                    if n_cst_refine <= 0:
                        iL += self.n_sec*3
                        found_key = 0
                        continue

                    for i in range(self.n_sec):

                        iL += 2
                        line1 = lines[iL].split()

                        iL += 1
                        line2 = lines[iL].split()

                        cst_ur = np.zeros(n_cst_refine)
                        cst_lr = np.zeros(n_cst_refine)

                        i1 = 0
                        i2 = 0
                        for j in range(n_cst_refine):
                            if i1<len(line1):
                                cst_ur[j] = float(line1[i1])
                                i1 += 1
                            if i2<len(line2):
                                cst_lr[j] = float(line2[i2])
                                i2 += 1

                        self.secs[i].cst_flip_u = cst_ur
                        self.secs[i].cst_flip_l = cst_lr

                    found_key = 0

                else:
                    # Lines that are not relevant
                    pass

                iL += 1
        
        print('Read surface [%s] settings'%(self.name))

        # Locate layout center for plot
        self.layout_center()

    def output_tecplot(self, fname=None, one_piece=False, split=False) -> None:
        '''
        Output the surface to `*.dat` in Tecplot format.

        Parameters
        ------------
        fname : str
            name of the output file.
        one_piece : bool
            if True, combine the span-wise sections into one piece.
        split : bool
            if True, split to upper and lower surfaces.
        '''
        if not split:
            super().output_tecplot(fname=fname, one_piece=one_piece)
            return

        if fname is None:
            fname = self.name + '.dat'

        n_sec   = 1 if self.l2d else self.n_sec-1
        n_piece = len(self.surfs)
        
        with open(fname, 'w') as f:
            f.write('Variables= X  Y  Z \n ')

            if not one_piece:

                for isec in range(n_piece):
                    
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    # surf_x[ns,nt], ns => spanwise
                    ns = self.ns
                    nt = int((surf_x.shape[1]+1)/2)

                    f.write('zone T="sec-u %d" i= %d j= %d \n'%(isec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                    f.write('zone T="sec-l %d" i= %d j= %d \n'%(isec, nt, ns))
                    for i in range(ns):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))

            else:
                
                npoint = n_sec*(self.ns-1) + 1
                ns = self.ns
                nt = int((self.surfs[0][0].shape[1]+1)/2)
                
                f.write('zone T="sec-u" i= %d j= %d \n'%(nt, npoint))

                for isec in range(n_piece):
                    
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]

                    if isec>=n_piece-1:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,j+nt-1], surf_y[i,j+nt-1], surf_z[i,j+nt-1]))

                f.write('zone T="sec-l" i= %d j= %d \n'%(nt, npoint))

                for isec in range(n_piece):
                    
                    surf_x = self.surfs[isec][0]
                    surf_y = self.surfs[isec][1]
                    surf_z = self.surfs[isec][2]
                    
                    if isec>=n_piece-1:
                        i_add = 0
                    else:
                        i_add = 1

                    for i in range(ns-i_add):
                        for j in range(nt):
                            f.write('  %.9f   %.9f   %.9f\n'%(surf_x[i,nt-1-j], surf_y[i,nt-1-j], surf_z[i,nt-1-j]))

    def output_plot3d(self, fname=None, split=False) -> None:
        '''
        Output the surface to `*.grd` in plot3d format.

        Parameters
        ------------
        fname : str
            name of the output file.
        split : bool
            if True, split to upper and lower surfaces.
        '''
        if not split:
            super().output_plot3d(fname=fname)
            return

        if fname is None:
            fname = self.name + '.grd'

        n_piece = len(self.surfs)

        # surf_x[ns,nt], ns => spanwise
        ns = self.ns

        with open(fname, 'w') as f:

            f.write('%d \n '%(n_piece*2))   # Number of surfaces
            for isec in range(n_piece):
                nt = int((self.surfs[isec][0].shape[1]+1)/2)
                f.write('%d %d 1\n '%(nt, ns))
                f.write('%d %d 1\n '%(nt, ns))

            for isec in range(n_piece):

                X = self.surfs[isec][0]
                Y = self.surfs[isec][1]
                Z = self.surfs[isec][2]
                nt = int((X.shape[1]+1)/2)

                #* Upper surface
                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(X[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Y[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Z[i,j+nt-1]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                #* Lower surface
                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(X[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Y[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')

                ii = 0
                for i in range(ns):
                    for j in range(nt):
                        f.write(' %.9f '%(Z[i,nt-1-j]))
                        ii += 1
                        if ii%3==0 or (i==ns-1 and j==nt-1):
                            f.write(' \n ')


#* ===========================================
#* Supportive functions
#* ===========================================

def surf_axisymmetric(xx: np.ndarray, yy: np.ndarray, phi0=0.0, phi1=360.0, ns=1001) -> List[np.ndarray]:
    '''
    Axisymmetric surface between curves.
    
    Parameters
    -----------
    xx, yy : ndarray
        generatrix profile
    phi0, phi1 : float
        angle (degree) about X-axis (X-Y plane is 0 degree)
    ns : int
        number of points in the interpolated direction.
    
    Returns
    ---------
    surf : list of ndarray
        coordinates of the surface, `[surf_x, surf_y, surf_z]`,
        `surf_x`'s shape is `[ns, nn]`. 

    Examples
    ------------
    >>> surf = surf_axisymmetric(xx, yy, phi0=0.0, phi1=360.0, ns=1001)

    '''
    nn = xx.shape[0]
    surf_x = np.zeros((ns,nn))
    surf_y = np.zeros((ns,nn))
    surf_z = np.zeros((ns,nn))
    zz = np.zeros(nn)

    for i in range(ns):

        tt    = 1.0*i/(ns-1.0)
        angle = (1-tt)*phi0 + tt*phi1

        surf_x[i,:], surf_y[i,:], surf_z[i,:] = rotate(xx, yy, zz, angle=angle, axis='X')

    return [surf_x, surf_y, surf_z]





