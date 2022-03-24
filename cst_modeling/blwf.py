'''
Interface for BLWF58

Building BLWF input file from wing geometry file and aircraft surface file.

'''
import os
import numpy as np
from .basic import (read_tecplot, intersect_surface_plane, interplot_from_curve,
                    rearrange_points, reconstruct_curve_by_length, BasicSurface)
from .foil import find_circle_3p


def output_curve(curve, fname='curves.dat', append=False):
    
    if append:
        f = open(fname, 'a')
    else:
        f = open(fname, 'w')
        f.write('Variables= X Y Z\n')
    
    if isinstance(curve, np.ndarray):
        
        n = curve.shape[0]
        if n>1:
            f.write('zone i= %d\n'%(n))
            for i in range(n):
                f.write('%20.10f  %20.10f  %20.10f\n'%(curve[i,0],curve[i,1],curve[i,2]))
            
    f.close()


class BLWF():
    '''
    Building input file for BLWF
    
    ### Inputs:
    ```text
    name:   project name
    ITH:    number of horizontal tail, 0, -1, 1, 2
    ITV:    number of vertical tail, 0, -1, 1
    INAC1:  number of inboard nacelle, 0, -1, 1
    INAC2:  number of outboard nacelle, 0, -1, 1
    IGU:    number of upper winglet, 0, -1, 1
    IGL:    number of lower winglet, 0, -1, 1
    ```
    '''
    def __init__(self, name='BLWF Aircraft', ITH=0, ITV=0,
                INAC1=0, INAC2=0, IGU=0, IGL=0):
        
        self.fname_ref = 'blwf-ref.in'
        
        self.name   = name
        self.ITH    = ITH
        self.ITV    = ITV
        self.INAC1  = INAC1
        self.INAC2  = INAC2
        self.IGU    = IGU
        self.IGL    = IGL
        
        #* Mesh plotting parameters for wing-body
        self.xmax_wb = 0.0      #  Coordinates of window boundaries for a plot of wing/body grid view
        self.xmin_wb = 0.0
        self.xmax_t  = 0.0      #  Coordinates of window boundaries for a plot of tail grid view
        self.xmin_t  = 0.0

        #* Wing sections
        self.FNS    = 0         # The number of span station, from root to tip, FNS < 51 
        self.ZLE    = [0.0 for _ in range(self.FNS)]    # Span location
        self.XLE    = [0.0 for _ in range(self.FNS)]    # X of leading edge
        self.YLE    = [0.0 for _ in range(self.FNS)]    # Y of leading edge
        self.CHORD  = [0.0 for _ in range(self.FNS)]    # Local chord
        self.THICK  = [1.0 for _ in range(self.FNS)]    # Y coordinates are multiplied by THICK
        self.EPSIL  = [0.0 for _ in range(self.FNS)]    # The angle through which the section is rotated to introduce twist.
                                                        # This angle is measured in the free stream direction
        self.FSEC_W = [1.0 for _ in range(self.FNS)]    # 1, the coordinates are read from the following lines
        self.YSYM   = [0.0 for _ in range(self.FNS)]    # 0, both upper and lower surfaces are provided, from nose to tail
        self.NU     = [0   for _ in range(self.FNS)]
        self.NL     = [0   for _ in range(self.FNS)]
        self.XSING  = [0.0 for _ in range(self.FNS)]    # Coordinates of the transformation singular line
        self.YSING  = [0.0 for _ in range(self.FNS)]    # XSING = XU[0]+0.5*RLE, YSING=YU[0]=YL[0], RLE = nose profile radius
        self.TRAIL  = [0.0 for _ in range(self.FNS)]    # The included angle at the trailing edge in degrees
        self.SLOPT  = [0.0 for _ in range(self.FNS)]    # The tangent of the slope of the mean camber line at the trailing edge
        self.XU     = [[0.0 for _ in range(self.NU[i])] for i in range(self.FNS)]
        self.YU     = [[0.0 for _ in range(self.NU[i])] for i in range(self.FNS)]
        self.XL     = [[0.0 for _ in range(self.NL[i])] for i in range(self.FNS)]
        self.YL     = [[0.0 for _ in range(self.NL[i])] for i in range(self.FNS)]   # The trailing edge point may be different if the profile has an open tail
        
        self.XLEW   = 0.0       #  X-coordinate of the leading edge of the first input wing section in the body-fitted coordinate system
        self.YLEW   = 0.0       #  X-coordinate of the leading edge of the first input wing section in the body-fitted coordinate system
        
        #* Body sections
        self.NSF    = 0         # number of the body intermediate sections, NSF < 60
        self.XLEF   = 0.0       # coordinate of the leading point of the fuselage 
        self.YLEF   = 0.0       # coordinate of the leading point of the fuselage
        self.XTEF   = 0.0       # coordinate of the trailing point of the fuselage
        self.YTEF   = 0.0       # coordinate of the trailing point of the fuselage
        self.XTEF0  = 0.0       # coordinate of the initiation of the semi-infinite body tail (for computation)
                                # XTEF0 < XTEF, XTEF0 ~ XTEF - 0.05*(XTEF-XLEF)
                                # the tail radius constitutes ~5%-20% of the body radius. 

        self.XF     = [0.0 for _ in range(self.NSF)]    # X location of the section
        self.YF     = [0.0 for _ in range(self.NSF)]    # Y coordinates of the fuselage mid point, = 0.5*(YSF[0]+YSF[-1])
        self.RF     = [0.0 for _ in range(self.NSF)]    # half-width of the fuselage section, = max(abs(ZSF))
        self.FSEC_B = [1.0 for _ in range(self.NSF)]    # 1, the coordinates are read from the following lines
        self.NS     = [0   for _ in range(self.NSF)]    # number of points, NS < 61
        self.YSF    = [[0.0 for _ in range(self.NS[i])] for i in range(self.NSF)]   # from negative to positive
        self.ZSF    = [[0.0 for _ in range(self.NS[i])] for i in range(self.NSF)]   # positive, symmetry is z=0

        #* Vertical tail sections
        self.XLETV  = 0.0       #  X-coordinate of the l.e. of the first input vertical tail section in the body-fitted coordinate system
        self.YLETV  = 0.0
        
        self.NC_VT  = 0         # number of the vertical tail input sections, NC_VT < 11
        self.XLE_VT = [0.0 for _ in range(self.NC_VT)]  # X cooridnate of the leading edge location of each section
        self.YLE_VT = [0.0 for _ in range(self.NC_VT)]  
        self.C_VT   = [0.0 for _ in range(self.NC_VT)]  # local chord length of each section
        self.T_VT   = [1.0 for _ in range(self.NC_VT)]  # thick, Y coordinates are multiplied by THICK
        self.ANT_VT = [0   for _ in range(self.NC_VT)]  # the index of the airfoil at the vertical tail section, ANT <= NA
        
        
        #* AIRFOIL DATA FOR NACELLES AND TAIL
        self.NA     = 0         # the number of airfoils that can be used for the nacelle, tail and winglet, NA<21
        self.A_YSYM = [0   for _ in range(self.NA)]     # 1 denotes a symmetric profile
        self.A_NU   = [0   for _ in range(self.NA)]
        self.A_NL   = [0   for _ in range(self.NA)]
        self.A_XSING= [0.0 for _ in range(self.NA)]     # Coordinates of the transformation singular line
        self.A_YSING= [0.0 for _ in range(self.NA)]     # XSING = XU[0]+0.5*RLE, YSING=YU[0]=YL[0], RLE = nose profile radius
        self.A_TRAIL= [0.0 for _ in range(self.NA)]     # The included angle at the trailing edge in degrees
        self.A_SLOPT= [0.0 for _ in range(self.NA)]     # The tangent of the slope of the mean camber line at the trailing edge
        self.A_XU   = [[0.0 for _ in range(self.A_NU[i])] for i in range(self.NA)]
        self.A_YU   = [[0.0 for _ in range(self.A_NU[i])] for i in range(self.NA)]
        self.A_XL   = [[0.0 for _ in range(self.A_NL[i])] for i in range(self.NA)]
        self.A_YL   = [[0.0 for _ in range(self.A_NL[i])] for i in range(self.NA)]    # The trailing edge point may be different if the profile has an open tail


        return

    @staticmethod
    def read_tecplot(fname='surface-aircraft.dat'):
        '''
        Read the baseline aircraft CFL3D result.
        
        ### Return: 
        ```text
        data:       list of ndarray [ni,nj,nk,nv], data of all zones
        name_var:   list, name of variables
        ```
        '''
        return read_tecplot(fname)
    
    def define_fuselage(self, zone_id: list, n_slice: int, n_point=0, extrapolate2sym=True,
                        fname='surface-aircraft.dat', index_xyz=[0,1,2]):
        '''
        Extract fuselage sections
        
        ### Inputs:
        ```text
        zone_id:    list, index of zones in the tecplot format file, start from 0
        n_slice:    number of sections in the X-axis
        n_point:    number of points to reconstruct the section curve, 0 means no reconstruction
        fname:      file name
        index_xyz:  index of variables in file for XYZ
        ```
        
        ### Return:
        ```text
        fuselage_sections:  list of ndarray [:,3]
        ```
        '''

        #* Read surface data
        data_, _ = read_tecplot(fname)
        if len(zone_id)==0:
            data = data_
        else:
            data = [data_[i] for i in zone_id]
        
        #* Extract range of fuselage
        x_min = 1000.0
        x_max =-1000.0
        y_min = 1000.0
        y_max =-1000.0
        z_min = 0.0
        for data_ in data:
            
            x_min = min(x_min, np.min(data_[:,:,:,index_xyz[0]]))
            x_max = max(x_max, np.max(data_[:,:,:,index_xyz[0]]))
            y_min = min(y_min, np.min(data_[:,:,:,index_xyz[1]]))
            y_max = max(y_max, np.max(data_[:,:,:,index_xyz[1]]))
            z_min = min(z_min, np.min(data_[:,:,:,index_xyz[2]]))
            
        lx = x_max-x_min
        ly = y_max-y_min
        
        #* Intersect sections
        Xs = np.linspace(x_min+0.01*lx, x_max-0.01*lx, n_slice, endpoint=True)
        fuselage_sections = []

        for kk in range(len(Xs)):

            P0 = np.array([Xs[kk], y_min-0.05*ly,   0.0])
            P1 = np.array([Xs[kk], y_min-0.05*ly, z_min])
            P3 = np.array([Xs[kk], y_max+0.05*ly,   0.0])
            
            curves = []
            xi_curves = []
            yt_curves = []
            
            for data_ in data:
                
                surface = np.concatenate((data_[:,:,:,index_xyz[0]:index_xyz[0]+1],
                            data_[:,:,:,index_xyz[1]:index_xyz[1]+1],
                            data_[:,:,:,index_xyz[2]:index_xyz[2]+1]), axis=3)
                surface = surface.squeeze()
                curve, _, xi_curve, yt_curve = intersect_surface_plane(surface, P0, P1, P3, within_bounds=True)
                
                if len(curve) > 0:
                    curves += curve
                    xi_curves += xi_curve.tolist()
                    yt_curves += yt_curve.tolist()

            _, old_index = rearrange_points(np.array(xi_curves), np.array(yt_curves), avg_dir=None)
            curve = np.array([curves[ii] for ii in old_index])
            
            #* Extrapolate to symmetry plane (Z=0)
            if extrapolate2sym:
                if abs(curve[0,2])>1e-3:
                    y_sym = curve[0,1] - (curve[0,1]-curve[1,1])/(curve[0,2]-curve[1,2])*curve[0,2]
                    curve = np.array([[curve[0,0], y_sym, 0.0]]+[curves[ii] for ii in old_index])
                elif abs(curve[-1,2])>1e-3:
                    y_sym = curve[-1,1] - (curve[-1,1]-curve[-2,1])/(curve[-1,2]-curve[-2,2])*curve[-1,2]
                    curve = np.array([curves[ii] for ii in old_index]+[[curve[-1,0], y_sym, 0.0]])
                
            #* Interploate section curves
            if n_point>0:
                curve = reconstruct_curve_by_length(curve, n_point)

            fuselage_sections.append(curve.copy())
            
        #* Update BLWF parameters
        self.xmax_wb = x_max + 0.05*lx
        self.xmin_wb = x_min - 0.05*lx
        
        self.NSF    = n_slice
        self.XLEF   = x_min
        self.XTEF   = x_max + 0.01*lx   # This additional 0.2*lx is to maintain the shape of body bottom
        self.XTEF0  = x_max - 0.02*lx
        
        min_x = 1000.0
        max_x =-1000.0

        for data_ in data:
            
            x_ = data_[:,:,:,index_xyz[0]]
            y_ = data_[:,:,:,index_xyz[1]]
            i_max = np.unravel_index(np.argmax(x_, axis=None), x_.shape)
            i_min = np.unravel_index(np.argmin(x_, axis=None), x_.shape)
            
            if x_[i_min] < min_x:
                min_x = x_[i_min]
                self.YLEF = y_[i_min]
            
            if x_[i_max] > max_x:
                max_x = x_[i_max]
                self.YTEF = y_[i_max]
                
        self.FSEC_B = [1.0 for _ in range(self.NSF)]
        
        for i in range(self.NSF):
            
            curve = fuselage_sections[i]    # ndarray [n,3]
            
            if curve[0,1] > curve[-1,1]:
                curve = np.flip(curve, axis=0)
            
            self.NS.append(curve.shape[0])
            self.XF.append(Xs[i])
            self.YF.append(0.5*(curve[0,1]+curve[-1,1]))
            self.RF.append(np.max(np.abs(curve[:,2])))

            self.YSF.append([])
            self.ZSF.append([])
            for k in range(curve.shape[0]):
                self.YSF[i].append(curve[k,1])
                self.ZSF[i].append(abs(curve[k,2]))

        return fuselage_sections

    def define_vertical_tail(self, zone_id: list, n_slice: int, n_point=0, avg_dir=None, 
                        fname='surface-aircraft.dat', index_xyz=[0,1,2], ratio=0.1):
        '''
        Extract vertical tail sections
        
        ### Inputs:
        ```text
        zone_id:    list, index of zones in the tecplot format file, start from 0
        n_slice:    number of sections in the X-axis
        n_point:    number of points to reconstruct the section curve, 0 means no reconstruction
        avg_dir:    ndarray [2], specified average direction
        fname:      file name
        index_xyz:  index of variables in file for XYZ
        ```
        
        ### Return:
        ```text
        tail_sections:  list of ndarray [:,3]
        ```
        '''
        #* Read surface data
        data_, _ = read_tecplot(fname)
        if len(zone_id)==0:
            data = data_
        else:
            data = [data_[i] for i in zone_id]
        
        #* Extract range of vertical tail
        x_min = 1000.0
        x_max =-1000.0
        y_min = 1000.0
        y_max =-1000.0

        for data_ in data:
            
            x_min = min(x_min, np.min(data_[:,:,:,index_xyz[0]]))
            x_max = max(x_max, np.max(data_[:,:,:,index_xyz[0]]))
            y_min = min(y_min, np.min(data_[:,:,:,index_xyz[1]]))
            y_max = max(y_max, np.max(data_[:,:,:,index_xyz[1]]))

        lx = x_max-x_min
        ly = y_max-y_min
        
        #* Intersect sections
        Ys = np.linspace(y_min+ratio*ly, y_max-0.001*ly, n_slice, endpoint=True)
        tail_sections = []

        for kk in range(len(Ys)):

            P0 = np.array([x_min-0.01*lx, Ys[kk],     0.0])
            P1 = np.array([x_max+0.01*lx, Ys[kk],     0.0])
            P3 = np.array([x_min-0.01*lx, Ys[kk], -0.5*lx])
            
            curves = []
            xi_curves = []
            yt_curves = []
            
            for data_ in data:
                
                surface = np.concatenate((data_[:,:,:,index_xyz[0]:index_xyz[0]+1],
                            data_[:,:,:,index_xyz[1]:index_xyz[1]+1],
                            data_[:,:,:,index_xyz[2]:index_xyz[2]+1]), axis=3)
                surface = surface.squeeze()
                curve, _, xi_curve, yt_curve = intersect_surface_plane(surface, P0, P1, P3, within_bounds=True)
                
                if len(curve) > 0:
                    curves += curve
                    xi_curves += xi_curve.tolist()
                    yt_curves += yt_curve.tolist()

            _, old_index = rearrange_points(np.array(xi_curves), np.array(yt_curves), avg_dir=avg_dir)
            curve = np.array([curves[ii] for ii in old_index])
            
            #* Interploate section curves
            if n_point>0:
                curve = reconstruct_curve_by_length(curve, n_point)

            tail_sections.append(curve.copy())
            
        #* Update BLWF parameters
        if self.ITV != 0:

            self.NC_VT = n_slice
            
            for i in range(n_slice):
                
                curve = tail_sections[i].copy()
                
                self.XLE_VT.append(curve[0,0])
                self.YLE_VT.append(curve[0,1])
                self.C_VT.append(curve[-1,0]-curve[0,0])
                self.T_VT.append(1.0)
                self.ANT_VT.append(i+1)
                
                self.NA += 1
                self.A_YSYM.append(1)
                self.A_NU.append(curve.shape[0])
                self.A_NL.append(curve.shape[0])
                
                curve = curve/self.C_VT[-1]
                curve[:,0] = curve[:,0] - curve[0,0]
                curve[:,2] = -(curve[:,2] - curve[0,2])
                
                rle, te_angle, _ = BLWF.airfoil_info(curve[:,0], curve[:,2], -curve[:,2])
                
                self.A_XSING.append(0.5*rle)
                self.A_YSING.append(0.0)
                self.A_TRAIL.append(te_angle)
                self.A_SLOPT.append(0.0)
                
                self.A_XU.append(curve[:,0].tolist())
                self.A_YU.append(curve[:,2].tolist())
                self.A_XL.append([])
                self.A_YL.append([])
                
                
            self.XLE_VT[0] -= 2*ratio*(self.XLE_VT[-1]-self.XLE_VT[0])
            self.YLE_VT[0] -= 2*ratio*(self.YLE_VT[-1]-self.YLE_VT[0])
            self.C_VT[0]   -= 2*ratio*(self.C_VT[-1]-self.C_VT[0])
            
            self.XLETV = self.XLE_VT[0]
            self.YLETV = self.YLE_VT[0]

        return tail_sections

    def define_wing(self, zone_id: list, n_slice: int, n_point=0, avg_dir=None, ratio_z=0.01,
                        fname='surface-aircraft.dat', index_xyz=[0,1,2]):
        '''
        Extract wing or horizontal tail sections
        
        ### Inputs:
        ```text
        zone_id:    list, index of zones in the tecplot format file, start from 0
        n_slice:    number of sections in the X-axis
        n_point:    number of points to reconstruct the section curve, 0 means no reconstruction
        avg_dir:    ndarray [2], specified average direction
        ratio_z:    ratio to control the range of slices in Z-axis
        fname:      file name
        index_xyz:  index of variables in file for XYZ
        ```
        
        ### Return:
        ```text
        wing_sections:  list of ndarray [:,3]
        ```
        '''
        #* Read surface data
        data_, _ = read_tecplot(fname)
        if len(zone_id)==0:
            data = data_
        else:
            data = [data_[i] for i in zone_id]
        
        #* Extract range of vertical tail
        x_min = 1000.0
        x_max =-1000.0
        y_min = 1000.0
        y_max =-1000.0
        z_min = 0.0
        z_max =-1000.0
        for data_ in data:
            
            x_min = min(x_min, np.min(data_[:,:,:,index_xyz[0]]))
            x_max = max(x_max, np.max(data_[:,:,:,index_xyz[0]]))
            y_min = min(y_min, np.min(data_[:,:,:,index_xyz[1]]))
            y_max = max(y_max, np.max(data_[:,:,:,index_xyz[1]]))
            z_min = min(z_min, np.min(data_[:,:,:,index_xyz[2]]))
            z_max = max(z_max, np.max(data_[:,:,:,index_xyz[2]]))

        lx = x_max-x_min
        ly = y_max-y_min
        lz = z_max-z_min
        
        #* Intersect sections
        Zs = np.linspace(z_min+0.001*lz, z_max-ratio_z*lz, n_slice, endpoint=True)
        wing_sections = []

        for kk in range(len(Zs)):
            
            P0 = np.array([x_min-0.01*lx, y_min-0.01*ly, Zs[kk]])
            P1 = np.array([x_max+0.01*lx, y_min-0.01*ly, Zs[kk]])
            P3 = np.array([x_min-0.01*lx, y_max+0.01*ly, Zs[kk]])
            
            curves = []
            xi_curves = []
            yt_curves = []
            
            for data_ in data:
                
                surface = np.concatenate((data_[:,:,:,index_xyz[0]:index_xyz[0]+1],
                            data_[:,:,:,index_xyz[1]:index_xyz[1]+1],
                            data_[:,:,:,index_xyz[2]:index_xyz[2]+1]), axis=3)
                surface = surface.squeeze()
                curve, _, xi_curve, yt_curve = intersect_surface_plane(surface, P0, P1, P3, within_bounds=True)
                
                if len(curve) > 0:
                    curves += curve
                    xi_curves += xi_curve.tolist()
                    yt_curves += yt_curve.tolist()

            _, old_index = rearrange_points(np.array(xi_curves), np.array(yt_curves), avg_dir=avg_dir)
            curve = np.array([curves[ii] for ii in old_index])
            
            #* Interploate section curves
            if n_point>0:
                curve = reconstruct_curve_by_length(curve, n_point)

            wing_sections.append(curve.copy())

        return wing_sections

    def cst_wing(self, wing: BasicSurface):
        '''
        Define wing parameters with wing
        '''
        self.XLEW = wing.secs[0].xLE
        self.YLEW = wing.secs[0].yLE
        
        self.FNS = wing.n_sec
        
        self.XLE    = [wing.secs[i].xLE - wing.secs[0].xLE for i in range(self.FNS)]
        self.YLE    = [wing.secs[i].yLE - wing.secs[0].yLE for i in range(self.FNS)]
        self.ZLE    = [wing.secs[i].zLE     for i in range(self.FNS)]
        self.CHORD  = [wing.secs[i].chord   for i in range(self.FNS)]
        
        self.EPSIL  = [-wing.secs[i].twist  for i in range(self.FNS)]
        
        self.THICK  = [1.0 for _ in range(self.FNS)]
        self.FSEC_W = [1.0 for _ in range(self.FNS)]
        self.YSYM   = [0.0 for _ in range(self.FNS)]
        self.NU     = [wing.secs[i].xx.shape[0] for i in range(self.FNS)]
        self.NL     = [wing.secs[i].xx.shape[0] for i in range(self.FNS)]
        
        self.XSING  = [0.5*wing.secs[i].RLE for i in range(self.FNS)]
        self.YSING  = [0.0 for _ in range(self.FNS)]
        
        self.TRAIL  = [wing.secs[i].te_angle for i in range(self.FNS)]
        self.SLOPT  = [wing.secs[i].te_slope for i in range(self.FNS)]
        
        for i in range(self.FNS):
            self.XU.append(wing.secs[i].xx.tolist())
            self.YU.append(wing.secs[i].yu.tolist())
            self.XL.append(wing.secs[i].xx.tolist())
            self.YL.append(wing.secs[i].yl.tolist())
    
    def update_input_wing(self, fname='blwf.in'):
        '''
        Update the wing sections in blwf.in
        '''
        lines = []
        with open(fname, 'r') as f0:
            lines = f0.readlines()
        
        for i in range(len(lines)):
            if 'WING/BODY DATA' in lines[i]:
                N1 = i-1
            if 'XLEF' in lines[i] and 'YLEF' in lines[i]:
                N3 = i
        
        f = open(fname, 'w')
        def wt(string):
            f.write(string+'\n')
        
        #* Line   1-215: fixed format, use original input file.
        for i in range(215):
            f.write(lines[i])
        
        #* Line 216-217
        wt('[  XLEW  ][  YLEW  ]')
        wt(' %9.4f %9.4f'%(self.XLEW, self.YLEW))
        
        #* Line 218-N1
        for i in range(217,N1):
            f.write(lines[i])
            
        #* Line: WING DATA
        wt('--------------------------------------------------------------')
        wt('      WING/BODY DATA')
        wt('[ FNS ]')           # The number of span station at which the wing sections 
        wt(' %.1f'%(self.FNS))  # are defined from the root to the wing tip. (FNS<51)
        
        N2 = N1+4
        for i in range(self.FNS):
            wt('[  ZLE   ][   XLE  ][   YLE  ][  CHORD ][  THICK ][  EPSIL ][  FSEC  ]')
            wt(' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'%(
                self.ZLE[i], self.XLE[i], self.YLE[i], self.CHORD[i], 
                self.THICK[i], self.EPSIL[i], self.FSEC_W[i]))
            wt('[  YSYM  ][   NU   ][   NL   ]')
            wt(' %9.4f %9.4f %9.4f'%(
                self.YSYM[i], self.NU[i], self.NL[i]))
            wt('[  XSING ][  YSING ][  TRAIL ][  SLOPT ]')
            wt(' %9.4f %9.4f %9.4f %9.4f'%(
                self.XSING[i], self.YSING[i], self.TRAIL[i], self.SLOPT[i]))
            wt('[   XU   ][   YU   ]')
            for k in range(self.NU[i]):
                wt(' %9.6f %9.6f'%(self.XU[i][k], self.YU[i][k]))
            wt('[   XL   ][   YL   ]')
            for k in range(self.NL[i]):
                wt(' %9.6f %9.6f'%(self.XL[i][k], self.YL[i][k]))
                
            N2 += 8 + self.NU[i] + self.NL[i]
        
        #* Line N3-END
        for i in range(N3, len(lines)):
            f.write(lines[i])
        
        f.close()

    def write_input_file(self, fname='blwf.in'):
        '''
        An example of BLWF input file
        '''
        lines = []
        with open(self.fname_ref, 'r') as f0:
            lines = f0.readlines()
        
        f = open(fname, 'w')
        def wt(string):
            f.write(string+'\n')
        
        #* Line   1-186: fixed format, use reference input file.
        wt('   %s'%(self.name))
        for i in range(185):
            f.write(lines[i+1])

        ii = 185
        #* Line 187-201: MESH PLOTTING PARAMETERS
        while ii<=200:
            ii += 1
            
            if ii==192:     #* Line 193: FOR WING-BODY
                # [  XMIN  ][  XMAX  ][  YMIN  ][  YMAX  ][  ZMIN  ][  ZMAX  ] 
                # YMAX-YMIN =XMAX-XMIN , ZMAX-ZMIN=1.5*(XMAX-XMIN)
                
                xmax = self.xmax_wb
                xmin = self.xmin_wb
                
                ymax = (xmax-xmin)/2.0
                ymin = - ymax
                zmin = 0.0
                zmax = 1.5*(xmax-xmin)
                
                wt(' %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f'%(
                    xmin, xmax, ymin, ymax, zmin, zmax))
                
            elif ii==197:   #* Line 198: FOR TAIL
                # [  XMIN  ][  XMAX  ][  YMIN  ][  YMAX  ][  ZMIN  ][  ZMAX  ]
                # YMAX-YMIN =XMAX-XMIN , ZMAX-ZMIN=1.5*(XMAX-XMIN)
                
                xmax = self.xmax_t
                xmin = self.xmin_t
                
                ymax = (xmax-xmin)/2.0
                ymin = - ymax
                zmin = 0.0
                zmax = 1.5*(xmax-xmin)
                
                wt(' %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f'%(
                    xmin, xmax, ymin, ymax, zmin, zmax))
                
            else:            
                f.write(lines[ii])

        #* Line 202-215: fixed format, use reference input file.
        while ii<=213:
            ii += 1
            f.write(lines[ii])
        
        #* Line 216-217
        wt('[  XLEW  ][  YLEW  ]')
        wt(' %9.4f %9.4f'%(self.XLEW, self.YLEW))

        N1 = 217

        #* Line (N1+1)-N2: HORIZONTAL TAIL MESH PARAMETERS AND HORIZONTAL TAIL POSITION
        wt('-------------------------------------------------------------')
        wt('      HORIZONTAL TAIL MESH PARAMETERS AND HORIZONTAL TAIL POSITION')
        wt('[  ITH   ]')
        if self.ITH == 0:
            # No horizontal tail
            N2 = N1+4
            wt('    0.    ')

        else:
            # ITH = -1, 1, 2
            N2 = N1+12
            wt('    1.    ')
            wt('[ NX_TH ][ NY_TH ][ NZ_TH ][ NT_TH ]')
            wt('    96.      6.       14.      10. ')
            wt('[ PZROOT ][ PZTIP ][ PXLE ][ PXTE ][ PYTE ] ')
            wt('   0.25     0.25      1.0     1.0     1.0 ')
            wt('[ XRB_TH ][ YRB1_TH][ YRB2_TH][ ZRB_TH ] ')
            wt('   0.5      0.5       1.0     1.4 ')
            wt('[ XLETH ][ YLETH ]')
            wt('   18.2     0.0 ')

        #* Line (N2+1)-N3: VERTICAL TAIL MESH PARAMETERS AND VERTICAL TAIL POSITION
        wt('--------------------------------------------------------------')
        wt('      VERTICAL TAIL MESH PARAMETERS AND VERTICAL TAIL POSITION')
        wt('[  ITV   ]')
        if self.ITV == 0:
            # No vertical tail
            N3 = N2+4
            wt('    0.    ')

        else:
            # ITV = -1, 1
            N3 = N2+12
            wt('    1.    ')
            wt('[ NX_TV  ][ NY_TV  ][ NZ_TV  ][ NT_TV ]')
            wt('   96.      8.        14.       10. ')
            wt('[ PZROOT ][ PZTIP  ][ PXLE   ][ PXTE   ][ PYTE   ]')
            wt('   0.4      0.4       1.0       1.0       1.0 ')
            wt('[ XRB_TV ][ ZRB_TV ][ YRB_TV ]')
            wt('   0.5      1.2       1.4 ')
            wt('[ XLETV  ][ YLETV  ]')
            wt(' %9.4f %9.4f'%(self.XLETV, self.YLETV))

        #* Line (N3+1)-N4: FIRST NACELLE MESH PARAMETERS AND NACELLE POSITION
        wt('--------------------------------------------------------------')
        wt('      FIRST NACELLE MESH PARAMETERS AND NACELLE POSITION')
        wt('[  INAC1 ]')
        if self.INAC1 == 0:
            # No nacelle
            N4 = N3+4
            wt('    0.    ')

        else:
            # INAC1 = -1, 1
            N4 = N3+14
            wt('    1.    ')
            wt('[ NYN ][ NZN ]')
            wt('   6.     8.  ')
            wt('[ NXNS ][ NXNW1 ][ NXNW2 ][ NXNA1 ][ NXNA2 ][ DXW1 ]')
            wt('   16.     8.       8.       2.       6.       2.4 ')
            wt('[ PXLEN ][ PYTEN ] [ PXTEN1 ][ PXTEN2 ][ PXNW1 ][ PXNW2 ]')
            wt('   1.       1.        2.        0.3       4.       15.   ')
            wt('[ XLERN ][ RB1 ][ RB2 ][ RB3 ][ RB4 ][ YOB ][ DYBW ][ DYBN ]')
            wt('   1.       1.    0.7    0.8    0.7    0.0    0.06    0.02  ')
            wt('[ XLEN ][ YLEN ][ ZLEN ][ NIUL ]')
            wt('  7.51    -0.4    4.0     -1.0  ')

        #* Line (N4+1)-N5: SECOND NACELLE MESH PARAMETERS AND NACELLE POSITION
        wt('--------------------------------------------------------------')
        wt('      SECOND NACELLE MESH PARAMETERS AND NACELLE POSITION')
        wt('[  INAC2 ]')
        if self.INAC2 == 0:
            # No second nacelle
            N5 = N4+4
            wt('    0.    ')

        else:
            # INAC2 = -1, 1
            N5 = N4+14
            wt('    1.    ')
            wt('[ NYN ][ NZN ]')
            wt('   6.     8.  ')
            wt('[ NXNS ][ NXNW1 ][ NXNW2 ][ NXNA1 ][ NXNA2 ][ DXW1 ]')
            wt('   16.     8.       8.       2.       6.       2.4 ')
            wt('[ PXLEN ][ PYTEN ] [ PXTEN1 ][ PXTEN2 ][ PXNW1 ][ PXNW2 ]')
            wt('   1.       1.        2.        0.3       4.       15.   ')
            wt('[ XLERN ][ RB1 ][ RB2 ][ RB3 ][ RB4 ][ YOB ][ DYBW ][ DYBN ]')
            wt('   1.       1.     1.    0.8    0.65   0.0    0.06    0.02  ')
            wt('[ XLEN ][ YLEN ][ ZLEN ][ NIUL ]')
            wt('  9.232   -0.3    7.0     -1.0  ')

        #* Line (N5+1)-N6: UPPER WINGLET MESH PARAMETERS. WINGLET POSITION
        wt('--------------------------------------------------------------')
        wt('      UPPER WINGLET MESH PARAMETERS. WINGLET POSITION')
        wt('[  IG ]')
        if self.IGU == 0:
            # No upper winglet
            N6 = N5+4
            wt('    0.    ')

        else:
            # IGU = -1, 1
            N6 = N5+12
            wt('    1.    ')
            wt('[ GAMMA ][ FI ]')
            wt('   1.      80. ')
            wt('[ XLEGW ][ NXLEGW ][ PXLEGW ]')
            wt('   0.2       10.      0.52   ')
            wt('[ XTEGW ][ NXTEGW ][ PXTEGW ]')
            wt('   0.9       6.0      0.3    ')
            wt('[ NYJTEG ][ NYJB ][ NZKB ][PZROOTG][ PZTIPG ]')
            wt('   6.0      3.0       3.0    0.1      0.3   ')

        #* Line (N6+1)-N7: LOWER WINGLET MESH PARAMETERS. WINGLET POSITION
        wt('--------------------------------------------------------------')
        wt('      LOWER WINGLET MESH PARAMETERS. WINGLET POSITION')
        wt('[  IG ]')
        if self.IGL == 0:
            # No lower winglet
            N7 = N6+4
            wt('    0.    ')

        else:
            # IGL = -1, 1
            N7 = N6+12
            wt('    1.    ')
            wt('[ GAMMA ][ FI ]')
            wt('   1.      80. ')
            wt('[ XLEGW ][ NXLEGW ][ PXLEGW ]')
            wt('   0.2       10.      0.52   ')
            wt('[ XTEGW ][ NXTEGW ][ PXTEGW ]')
            wt('   0.9       6.0      0.3    ')
            wt('[ NYJTEG ][ NYJB ][ NZKB ][PZROOTG][ PZTIPG ]')
            wt('   6.0      3.0       3.0    0.1      0.3   ')

        #* Line: WING DATA
        wt('--------------------------------------------------------------')
        wt('      WING/BODY DATA')
        wt('[ FNS ]')           # The number of span station at which the wing sections 
        wt(' %.1f'%(self.FNS))  # are defined from the root to the wing tip. (FNS<51)
        
        for i in range(self.FNS):
            wt('[  ZLE   ][   XLE  ][   YLE  ][  CHORD ][  THICK ][  EPSIL ][  FSEC  ]')
            wt(' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'%(
                self.ZLE[i], self.XLE[i], self.YLE[i], self.CHORD[i], 
                self.THICK[i], self.EPSIL[i], self.FSEC_W[i]))
            wt('[  YSYM  ][   NU   ][   NL   ]')
            wt(' %9.4f %9.4f %9.4f'%(
                self.YSYM[i], self.NU[i], self.NL[i]))
            wt('[  XSING ][  YSING ][  TRAIL ][  SLOPT ]')
            wt(' %9.4f %9.4f %9.4f %9.4f'%(
                self.XSING[i], self.YSING[i], self.TRAIL[i], self.SLOPT[i]))
            wt('[   XU   ][   YU   ]')
            for k in range(self.NU[i]):
                wt(' %9.6f %9.6f'%(self.XU[i][k], self.YU[i][k]))
            wt('[   XL   ][   YL   ]')
            for k in range(self.NL[i]):
                wt(' %9.6f %9.6f'%(self.XL[i][k], self.YL[i][k]))

        #* Line: BODY DATA
        # NSF: number of the body intermediate sections ( NSF<60.)
        wt('[  XLEF  ][  YLEF  ][  XTEF  ][  YTEF  ][  XTEF0 ][   NSF  ]')
        wt(' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'%(
            self.XLEF, self.YLEF, self.XTEF, self.YTEF, self.XTEF0, self.NSF))
        
        for i in range(self.NSF):
            wt('[  XF    ][  YF    ][  RF    ][  FSEC  ]')
            wt(' %9.4f %9.4f %9.4f %9.4f'%(
                self.XF[i], self.YF[i], self.RF[i], self.FSEC_B[i]))
            wt('[  NS    ]')
            wt(' %9.4f'%(self.NS[i]))
            wt('[  YSF   ][  ZSF   ]')
            for k in range(self.NS[i]):
                wt(' %9.6f %9.6f'%(self.YSF[i][k], self.ZSF[i][k]))
        
        #* Line: HORIZONTAL TAIL SECTION DATA
        wt('--------------------------------------------------------------')
        wt('      HORIZONTAL TAIL SECTION DATA')
        wt('[  NC ]')   # the number of the horisontal tail input sections (0<=NC<11)
        if self.ITH == 0:
            wt('    0.    ')

        else:
            wt('    2.    ')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ ANT ]')
            wt(' 0.000000  2179.2529 245.80607 253.95923  1.00000   0.0000      1.')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ ANT ]')
            wt(' 420.0000  2537.0166 303.78506 88.892090  1.00000   0.0000      2.')

        #* Line: VERTICAL TAIL SECTION DATA 
        wt('--------------------------------------------------------------')
        wt('      VERTICAL TAIL SECTION DATA')
        wt('[  NC    ]')   # the number of the vertical tail input sections (0<=NC<11)
        if self.ITV == 0:
            wt('    0.    ')

        else:
            wt(' %9.4f'%(self.NC_VT*1.0))
            for i in range(self.NC_VT):
                wt('[  ZLE   ][   XLE  ][   YLE  ][  CHORD ][  THICK ][  EPSIL ][  ANT   ]')
                wt(' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f'%(
                    0.0, self.XLE_VT[i], self.YLE_VT[i], self.C_VT[i], self.T_VT[i], 0.0, self.ANT_VT[i]))

        #* Line: SECTION DATA FOR FIRST NACELLE
        wt('--------------------------------------------------------------')
        wt('      SECTION DATA FOR FIRST NACELLE')
        wt('[  NC ]')
        if self.INAC1 == 0:
            wt('    0.    ')

        else:
            wt('    6.    ')
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('      0.       0.0    1.0517      3.572      1.       0.        1.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('     90.     0.118    1.0517      3.454      1.       0.        2.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    120.     0.118    1.0517      3.4199     1.       0.        4.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    180.     0.138    1.0517      3.384      1.       0.        3.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    240.     0.118    1.0517      3.4199     1.       0.        4.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    270.     0.118    1.0517      3.454      1.       0.        2.    ')  

        #* Line: SECTION DATA FOR SECOND NACELLE 
        wt('--------------------------------------------------------------')
        wt('      SECTION DATA FOR SECOND NACELLE')
        wt('[  NC ]')
        if self.INAC2 == 0:
            wt('    0.    ')

        else:
            wt('    6.    ')
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('      0.       0.0    1.0517      3.572      1.       0.        1.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('     90.     0.118    1.0517      3.454      1.       0.        2.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    120.     0.118    1.0517      3.4199     1.       0.        4.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    180.     0.138    1.0517      3.384      1.       0.        3.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    240.     0.118    1.0517      3.4199     1.       0.        4.    ')              
            wt('[   FN   ][  XNLE  ][  RNLE  ][ CHORDN ][  THICK ][   EN   ][   AN   ]')              
            wt('    270.     0.118    1.0517      3.454      1.       0.        2.    ')  

        #* Line: UPPER WINGLET SECTION DATA
        wt('--------------------------------------------------------------')
        wt('      UPPER WINGLET SECTION DATA')
        wt('[  NC ]')   # the number of the upper winglet input sections (0<=NC<11)
        if self.IGU == 0:
            wt('    0.    ')

        else:
            wt('    2.    ')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ AG ]')
            wt(' 0.000000  2179.2529 245.80607 253.95923  1.00000   0.0000      1.')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ AG ]')
            wt(' 420.0000  2537.0166 303.78506 88.892090  1.00000   0.0000      2.')
        
        #* Line: LOWER WINGLET SECTION DATA 
        wt('--------------------------------------------------------------')
        wt('      LOWER WINGLET SECTION DATA')
        wt('[  NC ]')   # the number of the lower winglet input sections (0<=NC<11)
        if self.IGL == 0:
            wt('    0.    ')

        else:
            wt('    2.    ')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ AG ]')
            wt(' 0.000000  2179.2529 245.80607 253.95923  1.00000   0.0000      1.')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ AG ]')
            wt(' 420.0000  2537.0166 303.78506 88.892090  1.00000   0.0000      2.')
        
        #* Line: AIRFOIL DATA FOR NACELLE, TAIL AND WINGLET SECTIONS.
        wt('--------------------------------------------------------------')
        wt('      AIRFOIL DATA FOR NACELLES AND TAIL ')
        wt('[  NA    ]')
        wt(' %9.4f'%(self.NA))
        if self.NA > 0:
            for i in range(self.NA):
                
                wt('[  YSYM  ][   NU   ][   NL   ]')
                wt(' %9.4f %9.4f %9.4f'%(self.A_YSYM[i], self.A_NU[i], self.A_NL[i]))
                wt('[  XSING ][ YSING  ][  TRAIL ][  SLOPT ]')
                wt(' %9.4f %9.4f %9.4f %9.4f'%(self.A_XSING[i], self.A_YSING[i], self.A_TRAIL[i], self.A_SLOPT[i]))
                wt('[   XU   ][   YU   ]')
                for k in range(self.A_NU[i]):
                    wt(' %9.6f %9.6f'%(self.A_XU[i][k], self.A_YU[i][k]))
                if self.A_YSYM[i]!=1:
                    wt('[   XL   ][   YL   ]')
                    for k in range(self.A_NL[i]):
                        wt(' %9.6f %9.6f'%(self.A_XL[i][k], self.A_YL[i][k]))

        f.close()

    @staticmethod
    def clean_result_files():
        '''
        
        '''
        try:
            os.system('del *.ps')
        except:
            pass
        
        files = ['blwf.sav', 'BLWF56.OUT', 'blwf.msh', 'blwf.lim', 'blwf.pl0', 'blwf.pl4', 'blwf.total',
                 'body_CF.plt', 'main_56.out', 'tp1.plt']
        
        for file in files:
            if os.path.exists(file):
                os.remove(file)
    
    @staticmethod
    def airfoil_info(xx, yu, yl):
        '''
        >>> rle, te_angle, te_slope = airfoil_info(xx, yu, yl)
        '''
        
        #* Calculate leading edge radius
        x_RLE = 0.005
        yu_RLE = interplot_from_curve(x_RLE, xx, yu)
        yl_RLE = interplot_from_curve(x_RLE, xx, yl)
        rle, _ = find_circle_3p([0.0,0.0], [x_RLE,yu_RLE], [x_RLE,yl_RLE])
        
        #* Trailing edge information
        a1 = [xx[-1]-xx[-5], yu[-1]-yu[-5]]
        a2 = [xx[-1]-xx[-5], yl[-1]-yl[-5]]
        te_angle = np.arccos(np.dot(a1,a2)/np.linalg.norm(a1)/np.linalg.norm(a2))/np.pi*180.0
        te_slope = 0.5*((yu[-1]+yl[-1])-(yu[-5]+yl[-5]))/(xx[-1]-xx[-5])
        
        return rle, te_angle, te_slope
        
