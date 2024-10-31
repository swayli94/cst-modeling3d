import numpy as np
import math
import os
import scipy.interpolate as spi
from scipy.integrate import odeint
from cst_modeling.io import plot3d_to_igs
deg = math.pi/180

def Conical_flow(M1:float, shockangle: float, n_out=200, T1=216.65, P1=1200):  
    # Input: Mach number , shock wave angle &  sound velocity
    r = 1.4
    R = 287
    def Taylor_Maccoll(v,t):
        r = 1.4 
        vr  = v[0].copy()
        vw = v[1].copy()
        dvdt = [vw, ((r-1)/2*(1-vr**2-vw**2)*(2*vr+vw/math.tan(t)) - vr*vw**2)/ ((r+1)/2*vw**2+(r-1)/2*vr**2-(r-1)/2)]
        return dvdt
    a1 = math.sqrt(r*R*T1)
    V1 = M1*a1
    V1t = V1*math.cos(shockangle)
    V1n = V1*math.sin(shockangle)
    M1n = M1*math.sin(shockangle)
    M2n = math.sqrt( (1+(r-1)/2*M1n**2)/ (r*M1n**2-(r-1)/2) )
    V2n = (2/(r+1)/M1n**2 + (r-1)/(r+1))*V1n
    V2t = V1t
    V2 = math.sqrt(V2t**2 + V2n**2)
    M2 = V2/V2n*M2n
    P2 = (1+2*r/(r+1)*(M1n**2-1))*P1
    P20 = P2*(1+(r-1)/2*M2**2)**(r/(r-1))
    T0 = (1+(r-1)/2*M1**2)*T1 # total Temperature
    T2 = T0/(1+(r-1)/2*M2**2) 
    Vmax = math.sqrt(V2**2+2*r*R*T2/(r-1)) # Vmax=sqrt(2h0) 
    # V0: nondimensionalized velocity after shock wave, V0=V2/Vmax
    V0 = 1/math.sqrt(2/(r-1)/M2**2+1)
    Vr0 = V2t/V2*V0
    Vw0 = -V2n/V2*V0
    wc = shockangle*0.99
    Vwc = Vw0
    while(Vwc <0 and wc>0):
        wc = wc-0.01*deg
        nt = 200
        thetas = np.linspace(shockangle,wc,nt)
        v = odeint(Taylor_Maccoll, [Vr0,Vw0], thetas)
        vr = (Vmax*v[:,0]).copy()
        vw = (Vmax*v[:,1]).copy()
        nt = len(vr)
        Vwc = vw[nt-1]
    thetas = np.linspace(shockangle,wc,n_out)
    vxs = np.zeros(n_out)
    vys = np.zeros(n_out)
    vs = np.zeros(n_out)
    Ms = np.zeros(n_out)
    Ts = np.zeros(n_out)
    Ps = np.zeros(n_out)
    mus = np.zeros(n_out)  # viscosity coefficient
    v = odeint(Taylor_Maccoll, [Vr0,Vw0], thetas)
    vrs = (Vmax*v[:,0]).copy()
    vws = (Vmax*v[:,1]).copy()
    for i in range(n_out):
        vxs[i] = vrs[i]*math.cos(thetas[i]) - vws[i]*math.sin(thetas[i])
        vys[i] = vrs[i]*math.sin(thetas[i]) + vws[i]*math.cos(thetas[i])
        vs[i] = math.sqrt(vxs[i]**2+vys[i]**2)
        Md = 1.01
        Mu = M2
        Ms[i] = (Mu+Md)/2
        for j in range(20):
            # In isentropic state, V3*sqrt(M3**(-2)+(r-1)/2) = V2*sqrt(M2**(-2)+(r-1)/2)
            if(vs[i]*math.sqrt(Ms[i]**(-2)+(r-1)/2) > V2*math.sqrt(M2**(-2)+(r-1)/2) ):                
                Md = Ms[i]
                Ms[i] = (Mu+Md)/2
            else:
                Mu = Ms[i]
                Ms[i] = (Mu+Md)/2    
        Ts[i] = T0/(1+(r-1)/2*Ms[i]**2)
        mus[i] = 1.458*10**(-6)*Ts[i]**1.5/(Ts[i]+110.4)
        Ps[i] = P20/(1+(r-1)/2*Ms[i]**2)**(r/(r-1))
    return thetas, vxs, vys

# streamline
def streamline(shockangle=30*math.pi/180, Rw=2.0, Ru=0.7, n_stream=1000, vxs=[], vys=[], thetas=[]):
    n_theta = len(thetas)
    x_stream = np.zeros(n_stream)
    y_stream = np.zeros(n_stream)
    x_stream[0] = (Rw-Ru)/math.tan(shockangle)
    y_stream[0] = Rw-Ru
    dx = Ru/math.tan(shockangle)/(n_stream-1)
    for i in range(n_stream-1):
        theta = math.atan(y_stream[i]/x_stream[i]) 
        j = 1
        while(thetas[j]>=theta and j<n_theta-1):  #thetas[] decreases monotonically 
            j = j+1
        vx = (thetas[j-1]-theta)/(thetas[j-1]-thetas[j])*vxs[j] + (theta-thetas[j])/(thetas[j-1]-thetas[j])*vxs[j-1]
        vy = (thetas[j-1]-theta)/(thetas[j-1]-thetas[j])*vys[j] + (theta-thetas[j])/(thetas[j-1]-thetas[j])*vys[j-1]
        dy = vy/vx*dx
        x_stream[i+1] = x_stream[i]+dx
        y_stream[i+1] = y_stream[i]+dy
    for i in range(1,n_stream):
        x_stream[i] = x_stream[i] - x_stream[0]
        y_stream[i] = y_stream[i] - y_stream[0]
    x_stream[0] = 0
    y_stream[0] = 0
    return x_stream,y_stream

'''''
n_arcs:    number of arcs
R_arc[]:   radius of each arc
phi_arc[]: angle span of each arc
R_arc[]:   radius of each arc          
xw0[],yw0[]:           left point position of each arc
phi0[]:                left point angle of each arc 
x_center[],y_center[]: center position of each arc
phiu[]: angle of each osculating plane
Ru[]: distance bitween wave and upline
Rw[];
zw[],yw[]:
zu[],yu[]:
''''' 
if __name__ == "__main__":
#* ============================================
#* Read Input
#* ============================================
    if True:
        with open('./files/input_wave.txt', 'r+') as f1:
            lines = f1.readlines()
            iL = 0
            while iL<len(lines):
                line = lines[iL].split()  
                if len(line) < 1:
                    iL += 1
                    continue              
                elif 'Arc_total' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    n_arcs = math.floor(float(line[0]))
                    iL += 1
                    phi_arc = np.zeros(n_arcs)
                    R_arc = np.zeros(n_arcs)         
                elif 'Phi' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    for i in range(n_arcs):
                        phi_arc[i] = float(line[i])*deg
                    iL += 1
                elif 'R' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    for i in range(n_arcs):
                        R_arc[i] = float(line[i])
                    iL += 1               
                else:
                    iL +=1

        with open('./files/Ma_w0.txt','r') as f:
            lines = f.readlines()
            Ma = float(lines[0].split()[0])
            shockangle = float(lines[0].split()[1])*deg

        
        with open('./files/upline.dat', 'r') as f1:
            lines = f1.readlines()
            nu = len(lines)
            zu = np.zeros(nu)
            yu = np.zeros(nu)
            for i in range(nu):
                line = lines[i].split()
                zu[i] = float(line[0])
                yu[i] = float(line[1])
    
#* ============================================
#* Main
#* ============================================
    if True:
        zw0 = np.zeros(n_arcs+1)
        yw0 = np.zeros(n_arcs+1)
        phi0 = np.zeros(n_arcs+1)
        z_center = np.zeros(n_arcs)
        y_center = np.zeros(n_arcs)    
        for i in range(n_arcs):
            phi0[i+1] = phi0[i] + phi_arc[i]
            zw0[i+1] = zw0[i] + (math.sin(phi0[i+1]) - math.sin(phi0[i]))*R_arc[i]
            yw0[i+1] = yw0[i] + (-math.cos(phi0[i+1]) + math.cos(phi0[i]))*R_arc[i]
            z_center[i] = zw0[i] - math.sin(phi0[i])*R_arc[i]
            y_center[i] = yw0[i] + math.cos(phi0[i])*R_arc[i]          
        for i in range(n_arcs):
            yw0[i] -= yw0[n_arcs]
            y_center[i] -= yw0[n_arcs]
        yw0[n_arcs] = 0
        # Width of upline input is modified to the value of shock wave line.
        for i in range(nu):
            zu[i] = zu[i]*zw0[n_arcs]/zu[nu-1]
            yu[i] = yu[i]*zw0[n_arcs]/zu[nu-1]
        zw = np.zeros(nu)
        yw = np.zeros(nu)
        phiu = np.zeros(nu)
        Ru = np.zeros(nu)
        Rw = np.zeros(nu)
        i = 0
        for j in range(nu-1):
            while(math.atan( (zu[j]-z_center[i])/(y_center[i]-yu[j]) ) -phi0[i] > phi_arc[i]):
                i += 1       
            phiu[j] = math.atan( (zu[j]-z_center[i])/(y_center[i]-yu[j]) )
            zw[j] = zw0[i] + (math.sin(phiu[j]) - math.sin(phi0[i]))*R_arc[i]
            yw[j] = yw0[i] + (-math.cos(phiu[j]) + math.cos(phi0[i]))*R_arc[i]
            Rw[j] = R_arc[i]
            Ru[j] = math.sqrt((zw[j]-zu[j])**2 + (yw[j]-yu[j])**2) 
        zw[nu-1] = zw0[n_arcs]
        Rw[nu-1] = R_arc[n_arcs-1]
        
        # Lower surface
        thetas, vxs, vys= Conical_flow(Ma, shockangle, n_out=1000)
        n_stream = 1000
        Lower_surface = np.zeros((nu, n_stream, 3))
        for i in range(nu-1):
            [x_stream,y_stream] = streamline(shockangle, Rw[i], Ru[i], n_stream, vxs, vys, thetas)
            x0 = -x_stream[n_stream-1]
            z0 = zu[i]
            y0 = yu[i]
            for j in range(n_stream):
                Lower_surface[i,j,0] = x0 + x_stream[j]
                Lower_surface[i,j,1] = y0 - y_stream[j]*math.cos(phiu[i])
                Lower_surface[i,j,2] = z0 + y_stream[j]*math.sin(phiu[i])
        for j in range(n_stream):
            Lower_surface[nu-1,j,2] = zu[nu-1]

        # Upper Surface
        Upper_surface = np.zeros((nu, n_stream, 3))
        for i in range(nu):
            for j in range(n_stream):
                Upper_surface[i,j,0] = Lower_surface[i,j,0]
                Upper_surface[i,j,1] = yu[i]         
                Upper_surface[i,j,2] = zu[i]
        
        # trailing
        n_trailing = 9
        Trailing_surface = np.zeros((nu, n_trailing, 3))
        zd = np.zeros(nu)    #Lower surface trailing curve
        yd = np.zeros(nu)   
        for i in range(nu):
            zd[i] = Lower_surface[i,n_stream-1,2]
            yd[i] = Lower_surface[i,n_stream-1,1]
        for i in range(nu):
            for j in range(n_trailing):
                Trailing_surface[i,j,0] = 0
                Trailing_surface[i,j,2] = (zu[i]-zd[i])/(n_trailing-1)*j + zd[i]
                Trailing_surface[i,j,1] = (yu[i]-yd[i])/(n_trailing-1)*j + yd[i]

        '''''
        connect
        '''''
        def arc_line(theta: float, length:float, n0: int):
            x = np.zeros(n0)
            y = np.zeros(n0)
            dtheta = theta/(n0-1)
            for i in range(n0):
                x[i] = math.sin(theta) - math.sin(theta-i*dtheta)
                y[i] = math.cos(theta-i*dtheta) - math.cos(theta) 
            x = np.multiply(x, length/x[n0-1]).copy()
            y = np.multiply(y, length/x[n0-1]).copy()
            return(x,y)

        # nomalize
        factor = 1/zu[nu-1]
        Lower_surface = Lower_surface*factor
        Upper_surface = Upper_surface*factor
        Trailing_surface = Trailing_surface*factor
        zu = zu*factor
        yu = yu*factor
        zw = zw*factor
        yw = yw*factor
        zd = zd*factor
        yd = yd*factor

        # Wave surface
        Wave_surface = np.zeros((nu,n_stream,3))
        for i in range(nu):
            Wave_surface[i,:,0] = np.linspace(Lower_surface[i,0,0], 0, n_stream)
            Wave_surface[i,:,1] = np.linspace(Lower_surface[i,0,1], yw[i], n_stream)
            Wave_surface[i,:,2] = np.linspace(Lower_surface[i,0,2], zw[i], n_stream)

        # Volume & Planform & Section
        planform = 0
        volume = 0
        ipo3 = spi.splrep(zu,yu,k=3)
        for i in range(nu-1):
            for j in range(n_stream-1):
                vector1 = np.array([Lower_surface[i+1,j+1,0]-Lower_surface[i,j,0],Lower_surface[i+1,j+1,2]-Lower_surface[i,j,2]])
                vector2 = np.array([Lower_surface[i,j+1,0]-Lower_surface[i+1,j,0],Lower_surface[i,j+1,2]-Lower_surface[i+1,j,2]])     
                thickness = np.sum(spi.splev(Lower_surface[i:i+2,j:j+2,2],ipo3))
                thickness -= np.sum(Lower_surface[i:i+2,j:j+2,1])
                thickness = thickness/4
                planform += abs(np.cross(vector1,vector2))/2
                volume += abs(np.cross(vector1,vector2))/2*thickness
        ipo3 = spi.splrep(zu,yu,k=3)
        yu2 = spi.splev(zd,ipo3)
        section = 0.0
        for i in range(nu-1):
            section += (yu2[i]-yd[i]+yu2[i+1]-yd[i+1])*(zd[i+1]-zd[i])/2

#* ============================================
#* Output
#* ============================================
    if True:
        with open('./dump/Planform_volume_section.txt','w+') as f:
            f.write('Planform  '+ str('%.4f' % planform) +'\n')
            f.write('Volume    '+ str('%.4f' % volume) +'\n')
            f.write('Section   '+ str('%.4f' % section) +'\n')
        
        # trailing
        with open('trailing_tec.dat', 'w+') as f2:
            f2.write('Variables= x y\n')
            f2.write('Zone T="Up"\n')
            for i in range(nu):
                f2.write( str(zu[i]) +'  '+ str(yu[i]) +'\n' )   
            f2.write('Zone T="Shock wave"\n')
            for i in range(nu):
                f2.write( str(zw[i]) +'  '+ str(yw[i]) +'\n' )
            f2.write('Zone T="Down"\n')
            for i in range(nu):
                f2.write( str(zd[i]) +'  '+ str(yd[i]) +'\n' )

        # tecplot
        with open('./dump/waverider_tec.dat','w+') as f4:
            f4.write('Variables= x y z\n')
            f4.write('Zone T="Lower surface" i='+ str(n_stream) +'  j='+ str(nu) +'\n')
            for i in range(nu):
                for j in range(n_stream):
                    f4.write(str(Lower_surface[i,j,0]) +'  '+ str(Lower_surface[i,j,1]) +'  '+ str(Lower_surface[i,j,2]) +'\n')
            f4.write('Zone T="Upper surface" i='+ str(n_stream) +'  j='+ str(nu) +'\n')
            for i in range(nu):
                for j in range(n_stream):
                    f4.write(str(Upper_surface[i,j,0]) +'  '+ str(Upper_surface[i,j,1]) +'  '+ str(Upper_surface[i,j,2]) +'\n')
            f4.write('Zone T="Trailing" i=' + str(n_trailing) + '  j='+ str(nu) +'\n')
            for i in range(nu):
                for j in range(n_trailing):
                    f4.write(str(Trailing_surface[i,j,0]) +'  '+ str(Trailing_surface[i,j,1]) +'  '+ str(Trailing_surface[i,j,2]) +'\n')     
        # grd
        scaley=910.00
        translationx=6400.00
        jump = 10
        n_x = math.floor(n_stream/jump)   
        index=0
        for i in range(nu):
            if Lower_surface[i][-1][2]*scaley >= 700 and index == 0:
                index=i
        ratio=(Lower_surface[index][-1][2]*scaley-700)/(Lower_surface[index][-1][2]-Lower_surface[index-1][-1][2])/scaley
        translationz=-(Lower_surface[index][-1][1]*(1-ratio)+Lower_surface[index-1][-1][1]*ratio)*scaley
        
        with open('./dump/lower_surface.xyz','w+') as f1:
            f1.write('1\n')
            f1.write(str(n_x) +'  '+ str(nu) + '  1\n')
            for i in [0,2,1]:
                for j in range(nu):
                    for k in range(n_x-1):
                        if i == 0 :
                            f1.write(str(Lower_surface[j][k*jump][i]*scaley+translationx))
                        elif i == 1:
                            f1.write(str(Lower_surface[j][k*jump][i]*scaley+translationz))
                        else:
                            f1.write(str(Lower_surface[j][k*jump][i]*scaley))
                        if (i*nu*n_x + j*n_x + k+1)%3 ==0:
                            f1.write('\n')
                        else:
                            f1.write('  ')
                    if i == 0 :
                        f1.write(str(Lower_surface[j][n_stream-1][i]*scaley+translationx))
                    elif i == 1:
                        f1.write(str(Lower_surface[j][n_stream-1][i]*scaley+translationz))
                    else:
                        f1.write(str(Lower_surface[j][n_stream-1][i]*scaley))
                    if (i*nu*n_x + (j+1)*n_x)%3 ==0:
                        f1.write('\n')
                    else:
                        f1.write('  ')
                                                
        with open('./dump/upper_surface.xyz','w+') as f1:
            f1.write('1\n')
            f1.write(str(n_x) +'  '+ str(nu) + '  1\n')
            for i in [0,2,1]:
                for j in range(nu):
                    for k in range(n_x-1):
                        if i == 0 :
                            f1.write(str(Upper_surface[j][k*jump][i]*scaley+translationx))
                        elif i == 1:
                            f1.write(str(Upper_surface[j][k*jump][i]*scaley+translationz))
                        else:
                            f1.write(str(Upper_surface[j][k*jump][i]*scaley))
                        if (i*nu*n_x + j*n_x + k+1)%3 ==0:
                            f1.write('\n')
                        else:
                            f1.write('  ')
                    if i == 0 :
                        f1.write(str(Upper_surface[j][n_stream-1][i]*scaley+translationx))
                    elif i == 1:
                        f1.write(str(Upper_surface[j][n_stream-1][i]*scaley+translationz))
                    else:
                        f1.write(str(Upper_surface[j][n_stream-1][i]*scaley))
                    if (i*nu*n_x + (j+1)*n_x)%3 ==0:
                        f1.write('\n')
                    else:
                        f1.write('  ')

        with open('./dump/TE.txt','w') as f1:
            for j in range(0,nu):
                f1.write(str(Lower_surface[j][-1][2]*scaley))
                f1.write('  ')
                f1.write(str(Lower_surface[j][-1][1]*scaley+translationz))
                f1.write('\n')

        with open('./dump/LE.txt','w') as f1:
            for j in range(0,nu):
                f1.write(str(Lower_surface[j][0][0]*scaley+translationx))
                f1.write('  ')
                f1.write(str(Lower_surface[j][0][2]*scaley))
                f1.write('  ')
                f1.write(str(Lower_surface[j][0][1]*scaley+translationz))
                f1.write('\n')
        with open('./dump/SLOPE.txt','w') as f1:
            slope=(Lower_surface[0][-1][1]-Lower_surface[0][-2][1])/(Lower_surface[0][-1][0]-Lower_surface[0][-2][0])
            f1.write(str(slope))

        # igs
        plot3d_to_igs(fname='./dump/upper_surface')
        plot3d_to_igs(fname='./dump/lower_surface')