'''
Interface for BLWF58

Building BLWF input file from wing geometry file and aircraft surface file.

'''

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
        
        self.name   = name
        self.ITH    = ITH
        self.ITV    = ITV
        self.INAC1  = INAC1
        self.INAC2  = INAC2
        self.IGU    = IGU
        self.IGL    = IGL
        
        self.xmax_wb = 0.0
        self.xmin_wb = 0.0
        self.xmax_t  = 0.0
        self.xmin_t  = 0.0

        #* Wing sections
        self.FNS    = 0
        self.ZLE    = [0.0 for _ in range(self.FNS)]
        self.XLE    = [0.0 for _ in range(self.FNS)]
        self.YLE    = [0.0 for _ in range(self.FNS)]
        self.CHORD  = [0.0 for _ in range(self.FNS)]
        self.THICK  = [0.0 for _ in range(self.FNS)]
        self.EPSIL  = [0.0 for _ in range(self.FNS)]
        self.FSEC_W = [0.0 for _ in range(self.FNS)]
        self.YSYM   = [0.0 for _ in range(self.FNS)]
        self.NU     = [0   for _ in range(self.FNS)]
        self.NL     = [0   for _ in range(self.FNS)]
        self.XSING  = [0.0 for _ in range(self.FNS)]
        self.YSING  = [0.0 for _ in range(self.FNS)]
        self.TRAIL  = [0.0 for _ in range(self.FNS)]
        self.SLOPT  = [0.0 for _ in range(self.FNS)]
        self.XU     = [[0.0 for _ in range(self.NU[i])] for i in range(self.FNS)]
        self.YU     = [[0.0 for _ in range(self.NU[i])] for i in range(self.FNS)]
        self.XL     = [[0.0 for _ in range(self.NL[i])] for i in range(self.FNS)]
        self.YL     = [[0.0 for _ in range(self.NL[i])] for i in range(self.FNS)]
        
        #* Body sections
        self.NSF    = 0
        self.XLEF   = 0.0
        self.YLEF   = 0.0 
        self.XTEF   = 0.0
        self.YTEF   = 0.0
        self.XTEF0  = 0.0
        self.XF     = [0.0 for _ in range(self.NSF)]
        self.YF     = [0.0 for _ in range(self.NSF)]
        self.RF     = [0.0 for _ in range(self.NSF)]
        self.FSEC_B = [0.0 for _ in range(self.NSF)]
        self.NS     = [0   for _ in range(self.NSF)]
        self.YSF    = [[0.0 for _ in range(self.NS[i])] for i in range(self.NSF)]
        self.ZSF    = [[0.0 for _ in range(self.NS[i])] for i in range(self.NSF)]

        return

    def read_surface(self, fname='surface-aircraft.dat'):
        '''
        Read the surface.dat of the baseline aircraft CFL3D result.
        '''


















        


    def write_input_file(self, fname='blwf.in'):
        '''
        An example of BLWF input file
        '''
        lines = []
        with open('blwf-ref.in', 'r') as f0:
            lines = f0.readlines()
        
        f = open(fname, 'w')
        def wt(string):
            f.write(string+'\n')
        
        #* Line   1-186: fixed format, use reference input file.
        wt('   %s'%(self.name))
        for i in range(185):
            wt(lines[i+1])

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
                wt(lines[ii])

        #* Line 202-217: fixed format, use reference input file.
        while ii<=216:
            ii += 1
            wt(lines[ii])

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
            wt('[ NX_TV ][ NY_TV ][ NZ_TV ][ NT_TV ]')
            wt('    96.      8.       14.      10. ')
            wt('[ PZROOT ][ PZTIP ][ PXLE ][ PXTE ][ PYTE ] ')
            wt('   0.4      0.4      1.0     1.0     1.0 ')
            wt('[ XRB_TV ][ ZRB_TV ][ YRB_TV ]  ')
            wt('   0.5      1.2     1.4 ')
            wt('[ XLETV ][ YLETV ]')
            wt('   18.2     0.0 ')

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
        wt(' %.0f'%(self.FNS))  # are defined from the root to the wing tip. (FNS<51)
        
        for i in range(self.FNS):
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ FSEC ]')
            wt(' %.6f %.6f %.6f %.6f %.6f %.6f %.6f'%(
                self.ZLE[i], self.XLE[i], self.YLE[i], self.CHORD[i], 
                self.THICK[i], self.EPSIL[i], self.FSEC_W[i]))
            wt('[  YSYM  ][   NU   ][   NL   ]')
            wt(' %.6f %.6f %.6f'%(
                self.YSYM[i], self.NU[i], self.NL[i]))
            wt('[  XSING ][  YSING ][  TRAIL ][  SLOPT ]')
            wt(' %.6f %.6f %.6f %.6f'%(
                self.XSING[i], self.YSING[i], self.TRAIL[i], self.SLOPT[i]))
            wt('[   XU   ][   YU   ]')
            for k in range(self.NU[i]):
                wt(' %.6f %.6f'%(self.XU[i][k], self.YU[i][k]))
            wt('[   XL   ][   YL   ]')
            for k in range(self.NL[i]):
                wt(' %.6f %.6f'%(self.XL[i][k], self.YL[i][k]))

        #* Line: BODY DATA
        # NSF: number of the body intermediate sections ( NSF<60.)
        wt('[  XLEF  ][  YLEF  ][  XTEF  ][  YTEF  ][  XTEF0 ][   NSF  ]')
        wt(' %.6f %.6f %.6f %.6f %.6f %.6f'%(
            self.XLEF, self.YLEF, self.XTEF, self.YTEF, self.XTEF0, self.NSF))
        
        for i in range(self.NSF):
            wt('[  XF    ][  YF    ][  RF    ][  FSEC  ]')
            wt(' %.6f %.6f %.6f %.6f'%(
                self.XF[i], self.YF[i], self.RF[i], self.FSEC_B[i]))
            wt('[  NS    ]')
            wt(' %.6f'%(self.NS[i]))
            wt('[  YSF   ][  ZSF   ]')
            for k in range(self.NS[i]):
                wt(' %.6f %.6f'%(self.YSF[i][k], self.ZSF[i][k]))
        
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
        wt('[  NC ]')   # the number of the vertical tail input sections (0<=NC<11)
        if self.ITV == 0:
            wt('    0.    ')

        else:
            wt('    2.    ')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ ANT ]')
            wt(' 0.000000  2179.2529 245.80607 253.95923  1.00000   0.0000      1.')
            wt('[ ZLE ][ XLE ][ YLE ][ CHORD ][ THICK ][ EPSIL ][ ANT ]')
            wt(' 0.000000  2537.0166 303.78506 88.892090  1.00000   0.0000      2.')

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
        wt('      AIRFOIL DATA FOR NACELLES AND TAIL (The following example can be deleted) ')
        wt('[  NA ]')
        wt('   1.0 ')
        wt('[  YSYM  ][   NU   ][   NL   ]')
        wt('1.00000  3.00000  3.00000 ')
        wt('[  XSING ][ YSING  ][  TRAIL ][  SLOPT ]')
        wt('0.00293  -0.00333  2.0000   -0.15000')
        wt('[   XU   ][   YU   ]')
        wt('0.00000     0.00000')
        wt('0.50000     0.05000')
        wt('1.00000     0.00000')
        wt('[   XL   ][   YL   ]')
        wt('0.00000     0.00000')
        wt('0.50000    -0.05000')
        wt('1.00000     0.00000')

        f.close()




