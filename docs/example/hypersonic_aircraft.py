from cmath import pi
from multiprocessing.connection import answer_challenge
from re import A
from statistics import harmonic_mean
from string import hexdigits
from cst_modeling.surface import OpenSurface
from cst_modeling.io import plot3d_to_igs
from scipy.interpolate import CubicSpline
import math
import os
import numpy as np

class Transition(OpenSurface):
    def trans(self, leader : list, tailer: list, q : float):
        # self(0:side , 1:bottom)
        #* Specify the leading edge and tailing edge
        ns = self.surfs[0][0].shape[0]
        nn = self.surfs[0][0].shape[1]
        for j in range(1, ns-1):
            self.surfs[0][0][j,0] = leader[j][0]
            self.surfs[0][1][j,0] = leader[j][1]
            self.surfs[0][2][j,0] = leader[j][2]

            self.surfs[0][0][j,-1] = tailer[j][0]
            self.surfs[0][1][j,-1] = tailer[j][1]
            self.surfs[0][2][j,-1] = tailer[j][2]
        # make sure that the boundary coordinate is monotonic
        # check in sequential order first
        #start the smooth transient
        # interpolation between the edge
        surfy=[[]for i in range(nn)]
        surfz=[[]for i in range(nn)]
        scaley=np.zeros(nn)
        scalez=np.zeros(nn)
        for i_surf in range(1):
            for i in range(1,nn-1):
                tt = 1.0*i/(nn-1)
                surfy[i] = tt*self.surfs[i_surf][1][:,-1]+(1-tt)*self.surfs[i_surf][1][:,0]
                surfz[i] = tt*self.surfs[i_surf][2][:,-1]+(1-tt)*self.surfs[i_surf][2][:,0]
                scaley[i] = (self.surfs[i_surf][1][-1,i]-self.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                scalez[i] = (self.surfs[i_surf][2][-1,i]-self.surfs[i_surf][2][0,i])/(surfz[i][-1]-surfz[i][0])
                self.surfs[i_surf][1][1:ns-1,i] = self.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                self.surfs[i_surf][2][1:ns-1,i] = self.surfs[i_surf][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
        # modify the transition to avoid intersect
        ratioy=[[] for i in range(nn)]
        ratioy_=ratioy.copy()
        modify=np.ones(ns)
        weight=[]
        position=int(np.floor(nn/21.0))
        height=q
        for i in range(nn):
            if i<=position:
                weight.append(height*(math.sin((math.pi/2/position)*i))**1.0)
            else:
                weight.append(height*(math.cos((math.pi/2/(nn-1-position))*(i-position)))**1.0)
        for i in range(ns):
            modify[i]=0.5*(1-(1-(i/51))**20)+0.5*((i/51)**(1/20))
        for k in range(nn):
            ratioy[k]=(self.surfs[0][2][:,k]-self.surfs[0][2][0,k])/(self.surfs[0][2][-1,k]-self.surfs[0][2][0,k])
        for k in range(nn):
            ratioy_[k]=ratioy[k]*(1-weight[k])+modify*weight[k]
        for k in range(1,nn-1):
            self.surfs[0][2][:,k]=ratioy_[k]*(self.surfs[0][2][-1,k]-self.surfs[0][2][0,k])+self.surfs[0][2][0,k]
def volume_caculation():
    ns=forebody.ns
    nn=forebody.nn
    volume_forebody=0.0
    volume_forebody2=0.0
    for i in range(nn):
        flag=0
        if forebody.surfs[0][1][ns-1,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] < forebody.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(index_i,ns):
                    dv=(forebody.surfs[0][2][-1,i]-forebody.surfs[0][2][j,i]) * (forebody.surfs[0][1][j,i]-forebody.surfs[0][1][j-1,i]) * (forebody.surfs[0][0][j,i+1]-forebody.surfs[0][0][j,i])
                    volume_forebody+=dv*1e-9
                for j in range(1,index_i):
                    dv=(forebody.surfs[0][2][-1,i]-forebody.surfs[0][2][j,i])*(forebody.surfs[0][1][j,i]-forebody.surfs[0][1][j-1,i])*(forebody.surfs[0][0][j,i+1]-forebody.surfs[0][0][j,i])
                    volume_forebody2+=dv*1e-9
        else:
            if i <= nn-2:
                for j in range(1,ns):
                    dv=(forebody.surfs[0][2][-1,i]-forebody.surfs[0][2][j,i])*(forebody.surfs[0][1][j,i]-forebody.surfs[0][1][j-1,i])*(forebody.surfs[0][0][j,i+1]-forebody.surfs[0][0][j,i])
                    volume_forebody2+=dv*1e-9
    for i in range(nn):
        flag=0
        if forebody.surfs[2][1][0,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] > forebody.surfs[2][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(1,index_i+1):
                    dv=(forebody.surfs[2][2][j,i]-forebody.surfs[0][2][-1,i]) * (forebody.surfs[2][1][j-1,i]-forebody.surfs[2][1][j,i]) * (forebody.surfs[2][0][j,i+1]-forebody.surfs[2][0][j,i])
                    volume_forebody+=dv*1e-9
                for j in range(index_i+1,ns):
                    dv=(forebody.surfs[2][2][j,i]-forebody.surfs[0][2][-1,i])*(forebody.surfs[2][1][j-1,i]-forebody.surfs[2][1][j,i])*(forebody.surfs[2][0][j,i+1]-forebody.surfs[2][0][j,i])
                    volume_forebody2+=dv*1e-9
        else:
            if i <= nn-2:
                for j in range(1,ns):
                    dv=(forebody.surfs[2][2][j,i]-forebody.surfs[0][2][-1,i])*(forebody.surfs[2][1][j-1,i]-forebody.surfs[2][1][j,i])*(forebody.surfs[2][0][j,i+1]-forebody.surfs[2][0][j,i])
                    volume_forebody2+=dv*1e-9
    ns=body1.ns
    nn=body1.nn
    volume_body1=0.0
    volume_body12=0.0
    for i in range(nn):
        flag=0
        if body1.surfs[0][1][ns-1,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] < body1.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(index_i,ns):
                    dv=(body1.surfs[0][2][-1,i]-body1.surfs[0][2][j,i]) * (body1.surfs[0][1][j,i]-body1.surfs[0][1][j-1,i]) * (body1.surfs[0][0][j,i+1]-body1.surfs[0][0][j,i])
                    volume_body1+=dv*1e-9
                for j in range(1,index_i):
                    dv=body1.surfs[0][2][-1,i] * (body1.surfs[0][1][j,i]-body1.surfs[0][1][j-1,i]) * (body1.surfs[0][0][j,i+1]-body1.surfs[0][0][j,i])
                    volume_body12+=dv*1e-9
    for i in range(nn):
        flag=0
        if body1.surfs[2][1][0,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] > body1.surfs[2][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(1,index_i+1):
                    dv=(body1.surfs[2][2][j,i]-body1.surfs[0][2][-1,i]) * (body1.surfs[2][1][j-1,i]-body1.surfs[2][1][j,i]) * (body1.surfs[2][0][j,i+1]-body1.surfs[2][0][j,i])
                    volume_body1+=dv*1e-9
                for j in range(index_i+1,ns):
                    dv=(body1.surfs[2][2][j,i]-body1.surfs[0][2][-1,i])*(body1.surfs[2][1][j-1,i]-body1.surfs[2][1][j,i])*(body1.surfs[2][0][j,i+1]-body1.surfs[2][0][j,i])
                    volume_body12+=dv*1e-9   
    ns=lowersurf.ns
    nn=lowersurf.nn
    volume_lowersurf=0.0
    volume_lowersurf2=0.0
    for i in range(nn):
        flag=0
        if lowersurf.surfs[0][1][ns-1,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] < lowersurf.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(index_i,ns):
                    dv=(wing_leading_z-lowersurf.surfs[0][2][j,i]) * (lowersurf.surfs[0][1][j,i]-lowersurf.surfs[0][1][j-1,i]) * (lowersurf.surfs[0][0][j,i+1]-lowersurf.surfs[0][0][j,i])
                    volume_lowersurf+=dv*1e-9
                for j in range(1,index_i):
                    dv=wing_leading_z * (lowersurf.surfs[0][1][j,i]-lowersurf.surfs[0][1][j-1,i]) * (lowersurf.surfs[0][0][j,i+1]-lowersurf.surfs[0][0][j,i])
                    volume_lowersurf2+=dv*1e-9

    ns=uppersurf.ns
    nn=uppersurf.nn
    volume_uppersurf=0.0
    volume_uppersurf2=0.0
    for i in range(nn):
        flag=0
        if uppersurf.surfs[0][1][0,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] > uppersurf.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(1,index_i+1):
                    dv=(uppersurf.surfs[0][2][j,i]-wing_leading_z) * (uppersurf.surfs[0][1][j-1,i]-uppersurf.surfs[0][1][j,i]) * (uppersurf.surfs[0][0][j,i+1]-uppersurf.surfs[0][0][j,i])
                    volume_uppersurf+=dv*1e-9
                for j in range(index_i+1,ns):
                    dv=(uppersurf.surfs[0][2][j,i]-wing_leading_z) * (uppersurf.surfs[0][1][j-1,i]-uppersurf.surfs[0][1][j,i]) * (uppersurf.surfs[0][0][j,i+1]-uppersurf.surfs[0][0][j,i])
                    volume_uppersurf2+=dv*1e-9

    ns=ending_lower.ns
    nn=ending_lower.nn
    volume_ending_lower=0.0
    volume_ending_lower2=0.0
    for i in range(nn):
        flag=0
        if ending_lower.surfs[0][1][ns-1,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] < ending_lower.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(index_i,ns):
                    dv=(wing_leading_z-ending_lower.surfs[0][2][j,i]) * (ending_lower.surfs[0][1][j,i]-ending_lower.surfs[0][1][j-1,i]) * (ending_lower.surfs[0][0][j,i+1]-ending_lower.surfs[0][0][j,i])
                    volume_ending_lower+=dv*1e-9    
                for j in range(1,index_i):
                    zz=max(0 , ending_lower.surfs[0][2][0,i])
                    dv=(wing_leading_z-zz) * (ending_lower.surfs[0][1][j,i]-ending_lower.surfs[0][1][j-1,i]) * (ending_lower.surfs[0][0][j,i+1]-ending_lower.surfs[0][0][j,i])
                    volume_ending_lower2+=dv*1e-9    
    volume_ending_lower2+=(wing_leading_z-ending_lower.surfs[0][2][index0,-1])**2*0.7307*inlet_edge[1]*1e-9
    ns=ending_upper.ns
    nn=ending_upper.nn
    volume_ending_upper=0.0
    volume_ending_upper2=0.0
    for i in range(nn):
        flag=0
        if ending_upper.surfs[0][1][0,i] > inlet_edge[1]:
            for j in range(1,ns):
                if inlet_edge[1] > ending_upper.surfs[0][1][j,i] and flag == 0:
                    index_i=j
                    flag=1
            if i <= nn-2:
                for j in range(1,index_i):
                    dv=(ending_upper.surfs[0][2][j,i]-wing_leading_z)*(ending_upper.surfs[0][1][j-1,i]-ending_upper.surfs[0][1][j,i])*(ending_upper.surfs[0][0][j,i+1]-ending_upper.surfs[0][0][j,i])
                    volume_ending_upper+=dv*1e-9
                if ending_upper.surfs[0][0][-1,i] < 23289:
                    for j in range(index_i+1,ns):
                        dv=(ending_upper.surfs[0][2][j,i]-wing_leading_z)*(ending_upper.surfs[0][1][j-1,i]-ending_upper.surfs[0][1][j,i])*(ending_upper.surfs[0][0][j,i+1]-ending_upper.surfs[0][0][j,i])
                        volume_ending_upper2+=dv*1e-9
                else:
                    for j in range(index_i+1,ns):
                        zz=400+0.6843*(ending_upper.surfs[0][0][-1,i]-23289)
                        dv=(ending_upper.surfs[0][2][j,i]-zz) * (ending_upper.surfs[0][1][j-1,i]-ending_upper.surfs[0][1][j,i]) * (ending_upper.surfs[0][0][j,i+1]-ending_upper.surfs[0][0][j,i])
                        volume_ending_upper2+=dv*1e-9
    volume_tank=volume_uppersurf     + volume_lowersurf + volume_body1 + volume_forebody  +  volume_ending_upper  +  volume_ending_lower + 0.0054 #0.0054 is a fixed correction
    volume_fuselage=volume_uppersurf2+ volume_lowersurf2+ volume_body12+ volume_forebody2 +  volume_ending_upper2 +  volume_ending_lower2
    volume_total=volume_tank+volume_fuselage
    return volume_total
def slope_upper_wing(x: float):
    # Make sure smooth transition to the upper surface of the wing
    y = 6.423*x**6 - 17.729*x**5 + 16.257*x**4 - 4.8813*x**3 + 0.2013*x**2 - 0.0308*x - 0.2506
    return y
def slope_lower_wing(x: float):
    # Make sure smooth transition to the lower surface of the wing
    y = 39.662*x**6 - 132.44*x**5 + 169.57*x**4 - 103.17*x**3 + 29.834*x**2 - 3.9306*x + 0.4914
    return y
def landing_gear_limit(x: float):
    # Make sure the swelling of the tank is enough to hold the landing gear
    y = -390.28*x**4 + 835.09*x**3 - 989.08*x**2 + 741.62*x + 1090.1
    z = -2324.9*x**4 + 4758.9*x**3 - 2025.7*x**2 - 359.56*x - 274.81
    return y, z
def top_crv1(x: float):
    # Define the top_curve(on the symetry plane) of forebody&body1 part
    scalez=vertex2_z-vertex1_z
    scalex=wing_leading_x-vertex1_x
    b1=scalex/scalez*math.tan(slope_top_forebody/180*pi)  # decided the slope of the top_line(leading edge)
    x_norm=(x-vertex1_x)/scalex
    z = (1-(1-x_norm)**b1)*scalez+vertex1_z+2.72
    return z
def top_crv2(x: float):
    # Define the top_curve(on the symetry plane) of uppersurf&ending_upper part
    x_norm=(x-wing_leading_x)/(tail_A[0][0]-wing_leading_x)
    cs = CubicSpline([0, 1], [1, 0], bc_type=((1, 0), (2, 0)))
    z_norm = cs(x_norm)
    z = z_norm*(vertex2_z-tail_A[0][2])+tail_A[0][2]
    return z
def shape_lower_wing(x: float):
    # Define the shape of the lower wing
    interp_x=np.linspace(0, 1, 51)
    interp_y=np.array([-2.662, -26.1626391, -45.93780588, -62.91401961, -77.81788335, -91.20650792, -103.4982956, -114.9993435, -125.9245618, -136.414878,\
        -146.5510201, -156.3652628, -165.8526109, -174.9807689, -182.7181361, -187.3205887, -190.4011362, -192.6163518, -194.2658047, -195.5142406,\
        -196.4649427, -197.185544, -197.7253063, -198.1191227, -198.3937713, -198.5674154, -198.6493841, -198.5791451, -198.1083454, -197.2071342,\
        -195.8679012, -194.079879, -191.8294678, -189.100568, -185.8749248, -182.132487, -177.8517631, -173.0101171, -167.5841866, -161.5501975,\
        -154.884257, -147.5626831, -139.5623883, -130.8612359, -121.4384627, -111.2749667, -100.3536475, -88.6597576, -76.181273, -62.90924934, -48.8382297])
    cs = CubicSpline(interp_x, interp_y)
    z = cs(x)
    return z

    return z
def build():
    if True:
    #* ============================================
    #* forebody
    #* ============================================
        if True:
            ns=forebody.ns
            nn=forebody.nn
            # produce and specify the shape of secs[0], which is the bottomline on symetry plane
            if True:
                forebody.secs[0].xLE=vertex1[0]
                forebody.secs[0].yLE=vertex1[1]
                forebody.secs[0].zLE=vertex1[2]-2.72  
                forebody.secs[0].xTE=inlet_edge[0]
                forebody.secs[0].yTE=0
                forebody.secs[0].zTE=0
                forebody.secs[0].section(nn=nn)
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    forebody.secs[0].y[j]=0
                    forebody.secs[0].x[j]=(1-tt)*forebody.secs[0].xLE+tt*forebody.secs[0].xTE
                    forebody.secs[0].z[j]=(1-tt)*forebody.secs[0].zLE+tt*forebody.secs[0].zTE
            # produce and specify the shape of secs[1], which is the leading edge1
            if True:
                forebody.secs[1].section(nn=nn)
                forebody.secs[1].x[0]=vertex1[0]
                forebody.secs[1].x[-1]=inlet_edge[0]
                forebody.secs[1].y[0]=vertex1[1]
                forebody.secs[1].z[0]=vertex1[2]-2.72
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    forebody.secs[1].x[j]=(1-tt)*forebody.secs[1].x[0]+tt*forebody.secs[1].x[-1]
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].y[:])
                forebody.secs[1].y[1:] = curv(forebody.secs[1].x[1:])
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].z[:])
                forebody.secs[1].z[1:] = curv(forebody.secs[1].x[1:])-2.72
            # produce and specify the shape of secs[2], which is the leading edge2
            if True:
                forebody.secs[2].section(nn=nn)
                forebody.secs[2].x[0]=vertex1[0]
                forebody.secs[2].x[-1]=inlet_edge[0]
                forebody.secs[2].y[0]=vertex1[1]
                forebody.secs[2].z[0]=vertex1[2]+2.72
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    forebody.secs[2].x[j]=(1-tt)*forebody.secs[2].x[0]+tt*forebody.secs[2].x[-1]
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].y[:])
                forebody.secs[2].y[1:] = curv(forebody.secs[2].x[1:])
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].z[:])
                forebody.secs[2].z[1:] = curv(forebody.secs[2].x[1:])+2.72
            # produce and specify the shape of secs[2], which is the topline on symetry plane
            if True:
                forebody.secs[3].section(nn=nn)
                forebody.secs[3].x[0]=vertex1[0]
                forebody.secs[3].x[-1]=inlet_edge[0]
                forebody.secs[3].z[0]=vertex1[2]+2.72
                for j in range(1, nn):
                    tt = 1.0*j/(nn-1)
                    forebody.secs[3].x[j]=(1-tt)*forebody.secs[3].x[0]+tt*forebody.secs[3].x[-1]
                    forebody.secs[3].y[j]=0.0
                    forebody.secs[3].z[j] = top_crv1(forebody.secs[3].x[j])
            # produce the geometry & obtain the number of secs between secs[1] and secs[2]
            forebody.geo(update_sec=False)
            # define the shape of the baseline ( upper part )
            if True:
                z0=forebody.secs[3].z[-1]
                y1=forebody.secs[2].y[-1]
                z1=forebody.secs[2].z[-1]
                scale=z0-z1
                shapez=np.zeros(ns)
                for j in range(1,ns-1):
                    y=forebody.surfs[2][1][j,-1]/y1
                    shapez[j]=(1-y**a_upper)**b_upper
                    forebody.surfs[2][2][j,-1]=shapez[j]*scale+z1
            # define the shape of the passivated leading edge
            if True:
                for i in range(0,nn):
                    for j in range(1,ns-1):
                        tt = 1.0*j/(ns-1)
                        forebody.surfs[1][0][j,i] -= -a_ellipse*(tt-0.5)**2+0.25*a_ellipse
            # define the shape of the baseline ( lower part )
            if True:
                index=int(np.floor(inlet_edge[1]/forebody.secs[1].y[-1]*ns))
                for j in range(index):
                    forebody.surfs[0][2][j,-1] = forebody.secs[0].z[-1]
                forebody.surfs[0][0][index,-1] = forebody.secs[0].x[-1]
                forebody.surfs[0][1][index,-1] = inlet_edge[1]
                forebody.surfs[0][2][index,-1] = forebody.secs[0].z[-1]
                z0=forebody.secs[0].z[-1]
                y1=forebody.secs[1].y[-1]
                z1=forebody.secs[1].z[-1]
                scalez=z1-z0
                scaley=y1-inlet_edge[1]
                k_norm = lower_k0/scalez*scaley
                cs = CubicSpline([0, lower_w0_norm,1], [0, lower_h0_norm,1], bc_type=((1, 0), (1, k_norm)))
                shapez=np.zeros(ns)
                for j in range(index+1,ns-1):
                    y=(forebody.surfs[0][1][j,-1]-inlet_edge[1])/scaley
                    shapez[j]=cs(y)
                    forebody.surfs[0][2][j,-1]=shapez[j]*scalez+z0
            # interpolation between the edge
            surfy=[[]for i in range(nn)]
            surfz=[[]for i in range(nn)]
            scaley=np.zeros(nn)
            scalez=np.zeros(nn)
            for i_surf in [0,2]:
                for i in range(1,nn-1):
                    tt = 1.0*i/(nn-1)
                    surfy[i] = tt*forebody.surfs[i_surf][1][:,-1]+(1-tt)*forebody.surfs[i_surf][1][:,0]
                    surfz[i] = tt*forebody.surfs[i_surf][2][:,-1]+(1-tt)*forebody.surfs[i_surf][2][:,0]
                    if surfy[i][-1]-surfy[i][0] ==0:
                        scaley[i] = 0
                    else:
                        scaley[i] = (forebody.surfs[i_surf][1][-1,i]-forebody.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                    if surfz[i][-1]-surfz[i][0] ==0:
                        scalez[i]=0
                    else:
                        scalez[i] = (forebody.surfs[i_surf][2][-1,i]-forebody.surfs[i_surf][2][0,i])/(surfz[i][-1]-surfz[i][0])
                    forebody.surfs[i_surf][1][1:ns-1,i] = forebody.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                    forebody.surfs[i_surf][2][1:ns-1,i] = forebody.surfs[i_surf][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
    #* ============================================
    #* body1
    #* ============================================
        if True:
            ns=body1.ns
            nn=body1.nn
            # produce and specify the shape of secs[0], which is the bottomline on symetry plane
            if True:
                body1.secs[0].section(nn=nn)
                body1.secs[0].x[0]=forebody.secs[0].x[-1]
                body1.secs[0].y[0]=forebody.secs[0].y[-1]
                body1.secs[0].z[0]=forebody.secs[0].z[-1]
                body1.secs[0].x[-1]=vertex2[0]
                for j in range(1,nn-1):
                    tt = 1.0*j/(nn-1)
                    body1.secs[0].x[j]=(1-tt)*body1.secs[0].x[0]+tt*body1.secs[0].x[-1]
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].y[:])
                body1.secs[0].y[1:nn] = curv(body1.secs[0].x[1:nn])
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].z[:])
                body1.secs[0].z[1:nn] = curv(body1.secs[0].x[1:nn])
            # produce and specify the shape of secs[1], which is the leading edge
            if True:
                body1.secs[1].section(nn=nn)
                body1.secs[1].x[0]=forebody.secs[1].x[-1]
                body1.secs[1].y[0]=forebody.secs[1].y[-1]
                body1.secs[1].z[0]=forebody.secs[1].z[-1]
                body1.secs[1].x[-1]=wing_leading[0]
                body1.secs[1].y[-1]=wing_leading[1]
                body1.secs[1].z[-1]=wing_leading[2]-2.72
                for j in range(1,nn-1):
                    tt = 1.0*j/(nn-1)
                    body1.secs[1].x[j]=(1-tt)*body1.secs[1].x[0]+tt*body1.secs[1].x[-1]
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].y[:])
                body1.secs[1].y[1:nn-1] = curv(body1.secs[1].x[1:nn-1])
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].z[:])
                body1.secs[1].z[1:nn-1] = curv(body1.secs[1].x[1:nn-1])-2.72
            # produce and specify the shape of secs[1], which is the leading edge
            if True:
                body1.secs[2].section(nn=nn)
                body1.secs[2].x[0]=forebody.secs[2].x[-1]
                body1.secs[2].y[0]=forebody.secs[2].y[-1]
                body1.secs[2].z[0]=forebody.secs[2].z[-1]
                body1.secs[2].x[-1]=wing_leading[0]
                body1.secs[2].y[-1]=wing_leading[1]
                body1.secs[2].z[-1]=wing_leading[2]+2.72
                for j in range(1,nn-1):
                    tt = 1.0*j/(nn-1)
                    body1.secs[2].x[j]=(1-tt)*body1.secs[2].x[0]+tt*body1.secs[2].x[-1]
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].y[:])
                body1.secs[2].y[1:nn-1] = curv(body1.secs[2].x[1:nn-1])
                curv = CubicSpline(support.secs[1].x[:], support.secs[1].z[:])
                body1.secs[2].z[1:nn-1] = curv(body1.secs[2].x[1:nn-1])+2.72
            # produce and specify the shape of secs[2], which is the topline on symetry plane
            if True:
                body1.secs[3].section(nn=nn)
                body1.secs[3].x[0]=forebody.secs[2].x[-1]
                body1.secs[3].z[0]=forebody.secs[2].z[-1]
                body1.secs[3].x[-1]=wing_leading_x
                body1.secs[3].z[-1]=vertex2_z
                for j in range(1,nn):
                    tt = 1.0*j/(nn-1)
                    body1.secs[3].x[j]=(1-tt)*body1.secs[3].x[0]+tt*body1.secs[3].x[-1]
                    body1.secs[3].y[j]=0.0
                    body1.secs[3].z[j] = top_crv1(body1.secs[3].x[j])
            # produce the geometry & obtain the number of secs between secs[1] and secs[2]
            body1.geo(update_sec=False)
            # define the shape in y-z plane ( lower part )
            if True:
                # front side
                body1.surfs[0][0][:,0]=forebody.surfs[0][0][:,-1]
                body1.surfs[0][1][:,0]=forebody.surfs[0][1][:,-1]
                body1.surfs[0][2][:,0]=forebody.surfs[0][2][:,-1]
                # other y-z curve
                lower_w1_norm=(lower_w1-inlet_edge[1])/(body1.secs[1].y[-1]-inlet_edge[1])
                lower_h1_norm=(lower_h1-body1.secs[0].z[-1])/(body1.secs[1].z[-1]-body1.secs[0].z[-1])
                for i in range(1,nn):
                    tt = 1.0*i/(nn-1)
                    index=int(np.floor(inlet_edge[1]/body1.secs[1].y[i]*ns))
                    for j in range(index):
                        body1.surfs[0][2][j,i] = body1.secs[0].z[i]
                    body1.surfs[0][0][index,i] = body1.secs[0].x[i]
                    body1.surfs[0][1][index,i] = inlet_edge[1]
                    body1.surfs[0][2][index,i] = body1.secs[0].z[i]
                    z0=body1.secs[0].z[i]
                    y1=body1.secs[1].y[i]
                    z1=body1.secs[1].z[i]
                    scalez=z1-z0
                    scaley=y1-inlet_edge[1]
                    k_norm=(tt*slope_lower_wing(0)+(1-tt)*lower_k0)/scalez*scaley
                    w_norm=tt*lower_w1_norm+(1-tt)*lower_w0_norm
                    h_norm=tt*lower_h1_norm+(1-tt)*lower_h0_norm
                    cs = CubicSpline([0, w_norm, 1], [0, max(h_norm,0), 1], bc_type=((1, 0), (1, k_norm)))
                    shapez=np.zeros(ns)
                    for j in range(index+1,ns-1):
                        y=(body1.surfs[0][1][j,i]-inlet_edge[1])/scaley
                        shapez[j]=cs(y)
                        body1.surfs[0][2][j,i]=shapez[j]*scalez+z0
            # define the shape of the passivated leading edge
            if True:
                # front side
                body1.surfs[1][0][:,0]=forebody.surfs[1][0][:,-1]
                body1.surfs[1][1][:,0]=forebody.surfs[1][1][:,-1]
                body1.surfs[1][2][:,0]=forebody.surfs[1][2][:,-1]
                # back side
                for i in range(1,nn):
                    for j in range(1,ns-1):
                        tt = 1.0*j/(ns-1)
                        body1.surfs[1][0][j,i] -= -a_ellipse*(tt-0.5)**2+0.25*a_ellipse
            # define the shape of the baseline ( upper part )
            if True:
                # front side
                body1.surfs[2][0][:,0]=forebody.surfs[2][0][:,-1]
                body1.surfs[2][1][:,0]=forebody.surfs[2][1][:,-1]
                body1.surfs[2][2][:,0]=forebody.surfs[2][2][:,-1]
                # back side
                z0=body1.secs[3].z[-1]
                y1=body1.secs[2].y[-1]
                z1=body1.secs[2].z[-1]
                scalez=z0-z1
                scaley=y1
                k_norm=slope_upper_wing(0)/scalez*scaley
                w_norm=upper_w1/scaley
                h_norm=(upper_h1-z1)/scalez
                wm_norm=upper_w1m/scaley
                hm_norm=(upper_h1m-z1)/scalez
                cs = CubicSpline([0, w_norm, wm_norm, 1], [1, h_norm, hm_norm, 0], bc_type=((1, 0), (1, k_norm)))
                shapez=np.zeros(ns)
                for j in range(1,ns-1):
                    y=body1.surfs[2][1][j,-1]/scaley
                    shapez[j]=cs(y)
                    body1.surfs[2][2][j,i]=shapez[j]*scalez+z1
            # interpolation between the edge, only the side and upper part, 
            # lower part & side part has already been determined, only need to build the upper part
            surfy=[[]for i in range(nn)]
            surfz=[[]for i in range(nn)]
            scaley=np.zeros(nn)
            scalez=np.zeros(nn)
            for i in range(1,nn-1):
                tt = 1.0*i/(nn-1)
                surfy[i] = tt*body1.surfs[2][1][:,-1]+(1-tt)*body1.surfs[2][1][:,0]
                surfz[i] = tt*body1.surfs[2][2][:,-1]+(1-tt)*body1.surfs[2][2][:,0]
                if surfy[i][-1]-surfy[i][0] ==0:
                    scaley[i] = 0
                else:
                    scaley[i] = (body1.surfs[2][1][-1,i]-body1.surfs[2][1][0,i])/(surfy[i][-1]-surfy[i][0])
                if surfz[i][-1]-surfz[i][0] ==0:
                    scalez[i] = 0
                else:
                    scalez[i] = (body1.surfs[2][2][-1,i]-body1.surfs[2][2][0,i])/(surfz[i][-1]-surfz[i][0])
                body1.surfs[2][1][1:ns-1,i] = body1.surfs[2][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                body1.surfs[2][2][1:ns-1,i] = body1.surfs[2][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
    #* ============================================
    #* uppersurf
    #* ============================================
        if True:
            ns=uppersurf.ns
            nn=uppersurf.nn
            # produce and specify the shape of secs[0], which is the side section which connect to the wing
            if True:
                uppersurf.secs[0].xLE=body1.secs[2].x[-1]
                uppersurf.secs[0].yLE=body1.secs[2].y[-1]
                uppersurf.secs[0].zLE=body1.secs[2].z[-1]          
                uppersurf.secs[0].xTE=wing_tailing_up[0]
                uppersurf.secs[0].yTE=wing_tailing_up[1]
                uppersurf.secs[0].zTE=wing_tailing_up[2]
                uppersurf.secs[0].chord = math.sqrt((uppersurf.secs[0].xTE-uppersurf.secs[0].xLE)**2+\
                                                    (uppersurf.secs[0].zTE-uppersurf.secs[0].zLE)**2)
                uppersurf.secs[0].section(nn=nn)
                for i in range(nn):
                    t = 1.0*i/(nn-1)
                    uppersurf.secs[0].x[i] = (1-t)*uppersurf.secs[0].xLE+t*uppersurf.secs[0].xTE
                    uppersurf.secs[0].y[i] = wing_leading_y
                    uppersurf.secs[0].z[i] = 474.58*t**6 - 1170.2*t**5 + 850.53*t**4 - 438.01*t**3 - 162.82*t**2 + 397.5*t + 2.7691 + wing_leading_z
            # produce and specify the shape of secs[1], which is the upper line on the symetry plane
            if True:
                uppersurf.secs[1].xLE=body1.secs[3].x[-1]
                uppersurf.secs[1].yLE=body1.secs[3].y[-1]
                uppersurf.secs[1].zLE=body1.secs[3].z[-1]          
                uppersurf.secs[1].xTE=wing_tailing_up[0]
                uppersurf.secs[1].yTE=0.0
                uppersurf.secs[1].zTE=top_crv2(wing_tailing_up[0])
                uppersurf.secs[1].section(nn=nn)
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    uppersurf.secs[1].y[j]=0
                    uppersurf.secs[1].x[j]=(1-tt)*uppersurf.secs[1].xLE+tt*uppersurf.secs[1].xTE
                    uppersurf.secs[1].z[j]=top_crv2(uppersurf.secs[1].x[j])
            # produce the geometry 
            uppersurf.geo(update_sec=False)
            # define the shape of the baseline ( upper part )
            if True:
                # front side
                uppersurf.surfs[0][0][:,0]=body1.surfs[2][0][:,-1]
                uppersurf.surfs[0][1][:,0]=body1.surfs[2][1][:,-1]
                uppersurf.surfs[0][2][:,0]=body1.surfs[2][2][:,-1]        
            # interpolation between the edge
            surfy=[[]for i in range(nn)]
            surfz=[[]for i in range(nn)]
            scaley=np.zeros(nn)
            scalez=np.zeros(nn)
            for i in range(1,nn):
                tt = 1.0*i/(nn-1)
                z0=uppersurf.surfs[0][2][-1,i]
                y1=uppersurf.surfs[0][1][0,i]
                z1=uppersurf.surfs[0][2][0,i]
                scalez=z0-z1
                scaley=y1
                k_norm=slope_upper_wing(tt)/scalez*scaley
                w_norm=(tt*upper_w2+(1-tt)*upper_w1)/scaley
                h_norm=(tt*upper_h2+(1-tt)*upper_h1-z1)/scalez
                wm_norm=(tt*upper_w2m+(1-tt)*upper_w1m)/scaley
                hm_norm=(tt*upper_h2m+(1-tt)*upper_h1m-z1)/scalez
                cs = CubicSpline([0, w_norm, wm_norm, 1], [1, h_norm, hm_norm, 0], bc_type=((1, 0), (1, k_norm)))
                shapez=np.zeros(ns)
                for j in range(1,ns-1):
                    y=uppersurf.surfs[0][1][j,-1]/scaley
                    shapez[j]=cs(y)
                    uppersurf.surfs[0][2][j,i]=shapez[j]*scalez+z1
    #* ============================================
    #* lowersurf
    #* ============================================
        if True:
            ns=lowersurf.ns
            nn=lowersurf.nn
            # produce and specify the shape of secs[0], which is the lower line on the symetry plane 
            if True:
                lowersurf.secs[0].section(nn=nn)
                lowersurf.secs[0].x[0]=body1.secs[0].x[-1]
                lowersurf.secs[0].y[0]=body1.secs[0].y[-1]
                lowersurf.secs[0].z[0]=body1.secs[0].z[-1]
                lowersurf.secs[0].x[-1]=wing_tailing_low[0]
                for j in range(1,nn-1):
                    tt = 1.0*j/(nn-1)
                    lowersurf.secs[0].x[j]=(1-tt)*lowersurf.secs[0].x[0]+tt*lowersurf.secs[0].x[-1]
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].y[:])
                lowersurf.secs[0].y[1:nn] = curv(lowersurf.secs[0].x[1:nn])
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].z[:])
                lowersurf.secs[0].z[1:nn] = curv(lowersurf.secs[0].x[1:nn])
            # produce and specify the shape of secs[1], which is the side section which connect to the wing
            if True:
                lowersurf.secs[1].xLE=body1.secs[1].x[-1]
                lowersurf.secs[1].yLE=body1.secs[1].y[-1]
                lowersurf.secs[1].zLE=body1.secs[1].z[-1]          
                lowersurf.secs[1].xTE=wing_tailing_low[0]
                lowersurf.secs[1].yTE=wing_tailing_low[1]
                lowersurf.secs[1].zTE=wing_tailing_low[2]
                lowersurf.secs[1].section(nn=nn)
                for i in range(nn):
                    t = 1.0*i/(nn-1)
                    lowersurf.secs[1].x[i] = (1-t)*lowersurf.secs[1].xLE+t*lowersurf.secs[1].xTE
                    lowersurf.secs[1].y[i] = wing_leading_y
                    lowersurf.secs[1].z[i] = shape_lower_wing(t)+ wing_leading_z
            
            # produce the geometry & obtain the number of secs between secs[1] and secs[2]
            lowersurf.geo(update_sec=False)
            # define the shape of the baseline 
            if True:
                # front side
                lowersurf.surfs[0][0][:,0]=body1.surfs[0][0][:,-1]
                lowersurf.surfs[0][1][:,0]=body1.surfs[0][1][:,-1]
                lowersurf.surfs[0][2][:,0]=body1.surfs[0][2][:,-1]
                # back side 
                for i in range(1,nn):
                    tt = 1.0*i/(nn-1)
                    index=int(np.floor(inlet_edge[1]/lowersurf.secs[1].y[i]*ns))
                    for j in range(index):
                        lowersurf.surfs[0][2][j,i] = lowersurf.secs[0].z[i]
                    lowersurf.surfs[0][0][index,i] = lowersurf.secs[0].x[i]
                    lowersurf.surfs[0][1][index,i] = inlet_edge[1]
                    lowersurf.surfs[0][2][index,i] = lowersurf.secs[0].z[i]
                    z0=lowersurf.secs[0].z[i]
                    y1=lowersurf.secs[1].y[i]
                    z1=lowersurf.secs[1].z[i]
                    scalez=z1-z0
                    scaley=y1-inlet_edge[1]
                    k_norm=slope_lower_wing(tt)/scalez*scaley
                    y_limit, z_limit= landing_gear_limit(tt)
                    w_norm=(y_limit-inlet_edge[1])/scaley
                    h_norm=(z_limit-z0)/scalez
                    cs = CubicSpline([0, w_norm, 1], [0, h_norm, 1], bc_type=((1, 0), (1, k_norm)))
                    shapez=np.zeros(ns)
                    for j in range(index+1,ns-1):
                        y=(lowersurf.surfs[0][1][j,i]-inlet_edge[1])/scaley
                        shapez[j]=cs(y)
                        lowersurf.surfs[0][2][j,i]=shapez[j]*scalez+z0
    #* ============================================
    #* ending_upper
    #* ============================================
        if True:
            ns=ending_upper.ns
            nn=ending_upper.nn
            refine_factor=int((ns-1)/(uppersurf.ns-1))
            # produce and specify the shape of secs[0], which is the side section
            if True:
                ending_upper.secs[0].xLE=uppersurf.secs[0].xTE
                ending_upper.secs[0].yLE=uppersurf.secs[0].yTE
                ending_upper.secs[0].zLE=uppersurf.secs[0].zTE        
                ending_upper.secs[0].xTE=tail_A[5][0]
                ending_upper.secs[0].yTE=tail_A[5][1]
                ending_upper.secs[0].zTE=tail_A[5][2]
                ending_upper.secs[0].section(nn=nn)
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    ending_upper.secs[0].x[j]=(1-tt)*ending_upper.secs[0].xLE+tt*ending_upper.secs[0].xTE
                    ending_upper.secs[0].y[j]=(1-tt)*ending_upper.secs[0].yLE+tt*ending_upper.secs[0].yTE
                    ending_upper.secs[0].z[j]=(1-tt)*ending_upper.secs[0].zLE+tt*ending_upper.secs[0].zTE
            # produce and specify the shape of secs[1], which is the upper line on the symetry plane
            if True:
                ending_upper.secs[1].xLE=uppersurf.secs[1].xTE
                ending_upper.secs[1].yLE=uppersurf.secs[1].yTE
                ending_upper.secs[1].zLE=uppersurf.secs[1].zTE          
                ending_upper.secs[1].xTE=tail_A[0][0]
                ending_upper.secs[1].yTE=tail_A[0][1]
                ending_upper.secs[1].zTE=tail_A[0][2]
                ending_upper.secs[1].section(nn=nn)
                for j in range(nn):
                    tt=1.0*j/(1.0*(nn-1))
                    ending_upper.secs[1].x[j]=(1-tt)*ending_upper.secs[1].xLE+tt*ending_upper.secs[1].xTE
                    ending_upper.secs[1].y[j]=0
                    ending_upper.secs[1].z[j]=top_crv2(ending_upper.secs[1].x[j])
            # produce the geometry & obtain the number of secs between secs[1] and secs[2]
            ending_upper.geo(update_sec=False)
            # define the shape of the baseline ( upper part )
            
            if True:
                # front side
                for k in range(0,3):
                    for i in range(int((ns-1)/refine_factor)):
                        for j in range(refine_factor):
                            tt=1.0*j/(1.0*refine_factor)
                            ending_upper.surfs[0][k][refine_factor*i+j,0]=(1-tt)*uppersurf.surfs[0][k][i,-1]+tt*uppersurf.surfs[0][k][i+1,-1]
                # back side
                delta_n=32
                index0=int(np.floor((ending_upper.secs[0].yTE-tail_A[4][1])/ending_upper.secs[0].yTE*ns))-delta_n
                index1=int(np.floor((ending_upper.secs[0].yTE-tail_A[2][1])/ending_upper.secs[0].yTE*ns))-delta_n
                index2=int(np.floor((ending_upper.secs[0].yTE-tail_A[1][1])/ending_upper.secs[0].yTE*ns))

                for j in range(0,index0):
                    tt=1.0*j/index0
                    ending_upper.surfs[0][0][j,-1]=tt*tail_A[4][0]+(1-tt)*tail_A[5][0]
                    ending_upper.surfs[0][1][j,-1]=tt*tail_A[4][1]+(1-tt)*tail_A[5][1]
                    ending_upper.surfs[0][2][j,-1]=tt*tail_A[4][2]+(1-tt)*tail_A[5][2]
                
                for j in range(index0,index1):
                    tt=1.0*(j-index0)/(index1-index0)
                    ending_upper.surfs[0][0][j,-1]=(1-tt)*tail_A[4][0]+tt*tail_A[3][0]
                    ending_upper.surfs[0][1][j,-1]=(1-tt)*tail_A[4][1]+tt*tail_A[3][1]
                    ending_upper.surfs[0][2][j,-1]=(1-tt)*tail_A[4][2]+tt*tail_A[3][2]
                for j in range(index1,index1+delta_n):
                    tt=1.0*(j-index1)/delta_n
                    ending_upper.surfs[0][0][j,-1]=(1-tt)*tail_A[3][0]+tt*tail_A[2][0]
                    ending_upper.surfs[0][1][j,-1]=(1-tt)*tail_A[3][1]+tt*tail_A[2][1]
                    ending_upper.surfs[0][2][j,-1]=(1-tt)*tail_A[3][2]+tt*tail_A[2][2]
                for j in range(index1+delta_n,index2):
                    tt=1.0*(j-index1-delta_n)/(index2-index1-delta_n)
                    ending_upper.surfs[0][0][j,-1]= -52.149*tt**4 + 73.261*tt**3 - 35.161*tt**2 + 60.192*tt + 23965
                    ending_upper.surfs[0][1][j,-1]=-19.277*tt**4 + 5.4014*tt**3 - 8.2513*tt**2 - 1.178*tt + 703.02
                    ending_upper.surfs[0][2][j,-1]= -38.948*tt**4 + 52.052*tt**3 - 20.553*tt**2 + 41.702*tt + 861.38
                for j in range(index2,ns):
                    tt=1.0*(j-index2)/(ns-index2-1)
                    ending_upper.surfs[0][0][j,-1]=(1-tt)*tail_A[1][0]+tt*tail_A[0][0]
                    ending_upper.surfs[0][1][j,-1]=(1-tt)*tail_A[1][1]+tt*tail_A[0][1]
                    ending_upper.surfs[0][2][j,-1]=(1-tt)*tail_A[1][2]+tt*tail_A[0][2]

            # interpolation between the edge
            surfy=[[]for i in range(nn)]
            surfz=[[]for i in range(nn)]
            scaley=np.zeros(nn)
            scalez=np.zeros(nn)
            for i_surf in range(1):
                for i in range(1,nn-1):
                    tt = 1.0*i/(nn-1)
                    surfy[i] = tt*ending_upper.surfs[i_surf][1][:,-1]+(1-tt)*ending_upper.surfs[i_surf][1][:,0]
                    surfz[i] = tt*ending_upper.surfs[i_surf][2][:,-1]+(1-tt)*ending_upper.surfs[i_surf][2][:,0]
                    scaley[i] = (ending_upper.surfs[i_surf][1][-1,i]-ending_upper.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                    scalez[i] = (ending_upper.surfs[i_surf][2][-1,i]-ending_upper.surfs[i_surf][2][0,i])/(surfz[i][-1]-surfz[i][0])
                    ending_upper.surfs[i_surf][1][1:ns-1,i] = ending_upper.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                    ending_upper.surfs[i_surf][2][1:ns-1,i] = ending_upper.surfs[i_surf][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
                    ending_upper.surfs[i_surf][0][1:ns-1,i] = (1-tt)*ending_upper.surfs[i_surf][0][1:ns-1,0]+tt*ending_upper.surfs[i_surf][0][1:ns-1,-1]
    #* ============================================
    #* ending_lower
    #* ============================================
        if True:
            ns=ending_lower.ns
            nn=ending_lower.nn
            refine_factor=int((ns-1)/(lowersurf.ns-1))
            # produce and specify the shape of secs[0], which is the lower line on the symetry plane
            if True:
                ending_lower.secs[0].section(nn=nn)
                ending_lower.secs[0].x[0]=lowersurf.secs[0].x[-1]
                ending_lower.secs[0].y[0]=lowersurf.secs[0].y[-1]
                ending_lower.secs[0].z[0]=lowersurf.secs[0].z[-1]
                ending_lower.secs[0].x[-1]=tail_B[0][0]
                for j in range(1,nn-1):
                    tt = 1.0*j/(nn-1)
                    ending_lower.secs[0].x[j]=(1-tt)*ending_lower.secs[0].x[0]+tt*ending_lower.secs[0].x[-1]
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].y[:])
                ending_lower.secs[0].y[1:nn] = curv(ending_lower.secs[0].x[1:nn])
                curv = CubicSpline(support.secs[2].x[:], support.secs[2].z[:])
                ending_lower.secs[0].z[1:nn] = curv(ending_lower.secs[0].x[1:nn])
            # produce and specify the shape of secs[1], which is the side section
            if True:
                ending_lower.secs[1].xLE=lowersurf.secs[1].xTE
                ending_lower.secs[1].yLE=lowersurf.secs[1].yTE
                ending_lower.secs[1].zLE=lowersurf.secs[1].zTE          
                ending_lower.secs[1].xTE=tail_B[3][0]
                ending_lower.secs[1].yTE=tail_B[3][1]
                ending_lower.secs[1].zTE=tail_B[3][2]
                
                ending_lower.secs[1].chord = math.sqrt((ending_lower.secs[1].xTE-ending_lower.secs[1].xLE)**2+\
                                                    (ending_lower.secs[1].yTE-ending_lower.secs[1].yLE)**2)
                ending_lower.secs[1].section(nn=nn)
                angle1_xz_rad=math.atan((ending_lower.secs[1].zTE-ending_lower.secs[1].zLE)/(ending_lower.secs[1].xTE-ending_lower.secs[1].xLE))
                dz=ending_lower.secs[1].yy*math.cos(angle1_xz_rad)*ending_lower.secs[1].chord
                dx=ending_lower.secs[1].yy*math.sin(angle1_xz_rad)*ending_lower.secs[1].chord
                for j in range(nn):
                    tt = 1.0*j/(nn-1)
                    ending_lower.secs[1].y[j]=(1-tt)*ending_lower.secs[1].yLE+tt*ending_lower.secs[1].yTE
                    ending_lower.secs[1].x[j]=(1-tt)*ending_lower.secs[1].xLE+tt*ending_lower.secs[1].xTE
                    ending_lower.secs[1].z[j]=(1-tt)*ending_lower.secs[1].zLE+tt*ending_lower.secs[1].zTE
                    ending_lower.secs[1].x[j]=ending_lower.secs[1].x[j]-dx[j]
                    ending_lower.secs[1].z[j]=ending_lower.secs[1].z[j]+dz[j]
            # produce the geometry & obtain the number of secs between secs[1] and secs[2]
            ending_lower.geo(update_sec=False)
            # define the shape of the baseline ( lower part )
            if True:
                # front side
                for k in range(0,3):
                    for i in range(int((ns-1)/refine_factor)):
                        for j in range(refine_factor):
                            tt=1.0*j/(1.0*refine_factor)
                            ending_lower.surfs[0][k][refine_factor*i+j,0]=(1-tt)*lowersurf.surfs[0][k][i,-1]+tt*lowersurf.surfs[0][k][i+1,-1]
                for k in range(0,3):
                    ending_lower.surfs[0][k][-1,0]=lowersurf.surfs[0][k][-1,-1]
                # back side
                # A parametic section of baseline which can control the shape of wing-body connection
                # k: deside the angle(below) of the leading_edge
                index0=int(np.floor(tail_B[1][1]/ending_lower.secs[1].yTE*ns))
                index1=int(np.floor(tail_B[2][1]/ending_lower.secs[1].yTE*ns))
                for j in range(0,index0+1):
                    tt=1.0*j/index0
                    ending_lower.surfs[0][0][j,-1]=tail_B[0][0]
                    ending_lower.surfs[0][1][j,-1]=(1-tt)*tail_B[0][1]+tt*tail_B[1][1]
                    ending_lower.surfs[0][2][j,-1]=tail_B[0][2]
                for j in range(index0+1,index1+1):
                    tt=1.0*(j-index0)/(index1-index0)
                    ending_lower.surfs[0][0][j,-1]=tt*tail_B[2][0]+(1-tt)*tail_B[1][0]
                    ending_lower.surfs[0][1][j,-1]=tt*tail_B[2][1]+(1-tt)*tail_B[1][1]
                    ending_lower.surfs[0][2][j,-1]=tt*tail_B[2][2]+(1-tt)*tail_B[1][2]
                for j in range(index1+1,ns):
                    tt=1.0*(j-index1)/(ns-index1-1)
                    ending_lower.surfs[0][0][j,-1]=tt*tail_B[3][0]+(1-tt)*tail_B[2][0]
                    ending_lower.surfs[0][1][j,-1]=tt*tail_B[3][1]+(1-tt)*tail_B[2][1]
                    ending_lower.surfs[0][2][j,-1]=tt*tail_B[3][2]+(1-tt)*tail_B[2][2]
            # interpolation between the edge
            surfy=[[]for i in range(nn)]
            surfz=[[]for i in range(nn)]
            scaley=np.zeros(nn)
            scalez=np.zeros(nn)
            for i_surf in range(1):
                for i in range(1,nn-1):
                    tt = 1.0*i/(nn-1)
                    surfy[i] = tt*ending_lower.surfs[i_surf][1][:,-1]+(1-tt)*ending_lower.surfs[i_surf][1][:,0]
                    surfz[i] = tt*ending_lower.surfs[i_surf][2][:,-1]+(1-tt)*ending_lower.surfs[i_surf][2][:,0]
                    scaley[i] = (ending_lower.surfs[i_surf][1][-1,i]-ending_lower.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                    scalez[i] = (ending_lower.surfs[i_surf][2][-1,i]-ending_lower.surfs[i_surf][2][0,i])/(surfz[i][-1]-surfz[i][0])
                    ending_lower.surfs[i_surf][1][1:ns-1,i] = ending_lower.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                    ending_lower.surfs[i_surf][2][1:ns-1,i] = ending_lower.surfs[i_surf][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
                    ending_lower.surfs[i_surf][0][1:ns-1,i] = (1-tt)*ending_lower.surfs[i_surf][0][1:ns-1,0]+tt*ending_lower.surfs[i_surf][0][1:ns-1,-1]
    #* ============================================
    #* return
    #* ============================================  
        if True:
            return forebody, body1, uppersurf, lowersurf, ending_upper, ending_lower
def output():
        print('************Start to Output************')
        # output the 'bottom' geometry 
        bottom.output_tecplot(fname='bottom.dat', one_piece=False)
        bottom.output_plot3d(fname='bottom.xyz')
        plot3d_to_igs(fname='bottom')
        # Output the 'connect' geometry
        connect.output_tecplot(fname='connect.dat', one_piece=False)
        connect.output_plot3d(fname='connect.xyz')
        plot3d_to_igs(fname='connect')
        # output the 'support' geometry
        support.output_tecplot(fname='support.dat')
        support.output_plot3d(fname='support.xyz')
        plot3d_to_igs(fname='support')
        # output the 'forebody' geometry
        forebody.output_tecplot(fname='forebody.dat')
        forebody.output_plot3d(fname='forebody.xyz')
        plot3d_to_igs(fname='forebody')
        # output the 'body1' geometry
        body1.output_tecplot(fname='body1.dat')
        body1.output_plot3d(fname='body1.xyz')
        plot3d_to_igs(fname='body1')
        # output the 'uppersurf' geometry
        uppersurf.output_tecplot(fname='uppersurf.dat')
        uppersurf.output_plot3d(fname='uppersurf.xyz')
        plot3d_to_igs(fname='uppersurf')
        # output the 'lowersurf' geometry
        lowersurf.output_tecplot(fname='lowersurf.dat')
        lowersurf.output_plot3d(fname='lowersurf.xyz')
        plot3d_to_igs(fname='lowersurf')
        # output the 'ending_upper' geometry
        ending_upper.output_tecplot(fname='ending_upper.dat')
        ending_upper.output_plot3d(fname='ending_upper.xyz')
        plot3d_to_igs(fname='ending_upper')
        # output the 'ending_lower' geometry
        ending_lower.output_tecplot(fname='ending_lower.dat')
        ending_lower.output_plot3d(fname='ending_lower.xyz')
        plot3d_to_igs(fname='ending_lower')
        # output the 'side' geometry
        side.output_tecplot(fname='side.dat')
        side.output_plot3d(fname='side.xyz')
        plot3d_to_igs(fname='side')
        print('Construct Completed')
if __name__ == "__main__":
    general_ns=52
    n1=0.8
    n2=1.0
    # parameters in y-z plane
    b_upper=2   #  z=(1-y**a_upper)**b_upper

    upper_w1=680  # fixed
    upper_h1=800  # constraint by volume
    upper_w1m=1000  # fixed
    upper_h1m=425  # parameter, now fixed, 425 is a moderate number
    upper_w2=680  # fixed
    upper_h2=800  # constraint by volume
    upper_w2m=1000  # fixed
    upper_h2m=upper_h1m # = upper_h1m

    lower_k0=0.0
    lower_w0_norm=0.3  # fixed now
    lower_h0_norm=0.28  # waverider-side-edge parameter
    lower_w1=1094.5  # fixed because of landing-gear constraint
    lower_h1=-277.972 # still fixed because of landing-gear constraint now, but maybe modify in the future
    # x-y plane
    interp_x=6400 # fixed now
    # x-z plane
    slope_leading_tank = 1.7 #should bigger than 1.52, which is the slope angle of the original shape
    slope_ending_tank = 9.2
#* ============================================
#* Parameters
#* ============================================
    if True:
        with open('./files/inputHA.txt', 'r+') as f1:
            lines = f1.readlines()
            iL = 0
            while iL<len(lines):
                line = lines[iL].split()
                if len(line) < 1:
                    iL += 1
                    continue          
                elif 'forebody_length' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    forebody_length = float(line[0])
                    vertex1_x = 6400-forebody_length
                elif 'forebody_thickness' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    forebody_thickness = float(line[0])
                    a_upper = forebody_thickness
                elif 'forebody_width' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    forebody_width = float(line[0])
                    interp_y = forebody_width
                elif 'nose_height' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    nose_height = float(line[0])
                    vertex1_z = nose_height
                elif 'nose_upper_angle' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    nose_upper_angle = float(line[0])
                    slope_top_forebody = nose_upper_angle
                elif 'nose_side_angle' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    nose_side_angle = float(line[0])
                    slope_side_forebody = float(line[0])
                elif 'radius_LE' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    radius_LE = float(line[0])
                    a_ellipse = 2.72*2.72*2/radius_LE# passivation radius of the leading edge
                elif 'ideal_volume' in line[0]:
                    iL += 1
                    line = lines[iL].split()
                    ideal_volume = float(line[0])
                else:
                    iL += 1
#* ============================================
#* bottom(sec[0] is on the symetry, sec[1] is the side section)
#* ============================================
    if True:
        bottom = OpenSurface(n_sec=2, name='bottom', nn=101, ns=general_ns)
        bottom.read_setting('./files/fuselage.txt')
        ns=bottom.ns
        nn=bottom.nn
        # produce and specify the shape of secs[0], which is on the symetry plane
        if True:
            bottom.secs[0].xLE=11397.243
            bottom.secs[0].yLE=0.0
            bottom.secs[0].zLE=-1280.817
            bottom.secs[0].xTE=20652.15
            bottom.secs[0].yTE=0.0
            bottom.secs[0].zTE=-1400.233
            scalex=bottom.secs[0].xTE-bottom.secs[0].xLE
            bottom.secs[0].section(nn=nn)
            #cut to three section
            j_1=int(np.floor(nn/3))
            j_2=int(np.floor(nn*2/3))
            # parameter of the first section
            scalez=280
            scalex=bottom.secs[0].x[j_1]-bottom.secs[0].x[0]
            h=0.8 # should less than 1
            e1=8*h-6
            e2=11-16*h
            e3=8*h-4
            # parameter of the third section
            scalex2=bottom.secs[0].x[nn-1]-bottom.secs[0].x[j_2+1]
            scalez2=scalez+(bottom.secs[0].zTE-bottom.secs[0].zLE)
            k2=math.tan(10/180*math.pi)*scalex2/scalez2    # decide the slope of the tailing_edge, parameter
            for j in range(nn):
                tt = 1.0*j/(nn-1)
                if j<= j_1:
                    t = 1.0*j/(1.0*j_1)
                    bottom.secs[0].x[j]=(1-tt)*bottom.secs[0].xLE+tt*bottom.secs[0].xTE
                    bottom.secs[0].y[j]=(1-tt)*bottom.secs[0].yLE+tt*bottom.secs[0].yTE
                    bottom.secs[0].z[j]=bottom.secs[0].z[0]-(e1*t**3+e2*t**2+e3*t)*scalez
                elif j<=j_2:
                    bottom.secs[0].x[j]=(1-tt)*bottom.secs[0].xLE+tt*bottom.secs[0].xTE
                    bottom.secs[0].y[j]=(1-tt)*bottom.secs[0].yLE+tt*bottom.secs[0].yTE
                    bottom.secs[0].z[j]=bottom.secs[0].z[j-1]
                else:
                    t = 1.0*(j-j_2)/(1.0*(nn-j_2-1))
                    bottom.secs[0].x[j]=(1-tt)*bottom.secs[0].xLE+tt*bottom.secs[0].xTE
                    bottom.secs[0].y[j]=(1-tt)*bottom.secs[0].yLE+tt*bottom.secs[0].yTE
                    bottom.secs[0].z[j]=bottom.secs[0].z[j_2]+(t**k2)*scalez2
                #print(bottom.secs[0].x[j],  bottom.secs[0].z[j])
            #print('\n\n\n')
        # produce and specify the shape of secs[1], which is the side edge
        if True:       
            bottom.secs[1].xLE=11397.243
            bottom.secs[1].yLE=670.426
            bottom.secs[1].zLE=-1280.817
            bottom.secs[1].xTE=20652.15
            bottom.secs[1].yTE=680.209
            bottom.secs[1].zTE=-1400.233
            bottom.secs[1].section(nn=nn)
            for j in range(nn):
                tt = 1.0*j/(nn-1)
                bottom.secs[1].x[j]=(1-tt)*bottom.secs[1].xLE+tt*bottom.secs[1].xTE
                bottom.secs[1].y[j]=(1-tt)*bottom.secs[1].yLE+tt*bottom.secs[1].yTE
                bottom.secs[1].z[j]=bottom.secs[0].z[j]
        # produce the geometry & obtain the number of secs between secs[1] and secs[2]
        bottom.geo(update_sec=False)
        # define the shape of the baseline 
        if True:
            # front side
            bottom.surfs[0][0][:,0]=np.linspace( bottom.secs[0].xLE  ,  bottom.secs[1].xLE, 52 )
            bottom.surfs[0][1][:,0]=np.linspace( bottom.secs[0].yLE  ,  bottom.secs[1].yLE, 52 )
            bottom.surfs[0][2][:,0]=np.linspace( bottom.secs[0].zLE  ,  bottom.secs[1].zLE, 52 )
            # back side
            bottom.surfs[0][0][:,-1]=np.linspace( bottom.secs[0].xTE  ,  bottom.secs[1].xTE, 52 )
            bottom.surfs[0][1][:,-1]=np.linspace( bottom.secs[0].yTE  ,  bottom.secs[1].yTE, 52 )
            bottom.surfs[0][2][:,-1]=np.linspace( bottom.secs[0].zTE  ,  bottom.secs[1].zTE, 52 )
        # interpolation between the edge
        surfy=[[]for i in range(nn)]
        surfz=[[]for i in range(nn)]
        scaley=np.zeros(nn)
        scalez=np.zeros(nn)
        for i_surf in range(1):
            for i in range(1,nn-1):
                tt = 1.0*i/(nn-1)
                surfy[i] = tt*bottom.surfs[i_surf][1][:,-1]+(1-tt)*bottom.surfs[i_surf][1][:,0]
                scaley[i] = (bottom.surfs[i_surf][1][-1,i]-bottom.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                bottom.surfs[i_surf][1][1:ns-1,i] = bottom.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                bottom.surfs[i_surf][2][1:ns-1,i] = bottom.surfs[i_surf][2][0,i]
#* ============================================
#* connect
#* ============================================
    if True:
        connect = Transition(n_sec=2, name='connect', nn=101, ns=general_ns)
        ns=connect.ns
        nn=connect.nn
        #* modify the connect.sec[0].x,y,z(which should match the side surface) & connect.sec[1].x,y,z(which should match the bottom surface)
        connect.secs[0].section(nn=nn)
        connect.secs[1].section(nn=nn)
        connect.secs[0].x=np.linspace( 11408.557  ,  20678.591, nn )
        connect.secs[0].y=np.linspace( 705.969   ,  702.993, nn  )
        connect.secs[0].z=np.linspace( -1245.914  ,  -1378.991, nn )
        connect.secs[1].x=bottom.secs[1].x
        connect.secs[1].y=bottom.secs[1].y
        connect.secs[1].z=bottom.secs[1].z
        connect.geo(update_sec=False)

        #* specify the leading edge and tailing edge
        leader=[]
        for i in range(ns):
            t=1.0*i/(1.0*(ns-1))
            leaderx= -2.7992*t**4 + 23.001*t**3 - 28.67*t**2 - 2.859*t + 11409
            leadery= 4.7554*t**4 + 5.8241*t**3 - 47.219*t**2 + 1.1018*t + 705.97
            leaderz= -6.7261*t**4 + 25.934*t**3 + 2.9733*t**2 - 57.088*t - 1245.9 
            leader.append([leaderx,leadery,leaderz])
        tailer=[]
        for i in range(ns):
            t=1.0*i/(1.0*(ns-1))
            tailerx = 41.184*t**4 - 81.289*t**3 + 72.942*t**2 - 59.156*t + 20679
            tailery = 25.216*t**4 - 53.55*t**3 + 10.142*t**2 - 4.6272*t + 703.01
            tailerz = 32.78*t**4 - 57.835*t**3 + 44.792*t**2 - 40.868*t - 1379.1 
            tailer.append([tailerx,tailery,tailerz])

        # build the transition to connect two surface
        # q is the parameter to modify the transition to avoid intersect with inlet
        # q supposed to be 0.9, if there is still intersect, please increase q, but do not exceed 1.0
        connect.trans(leader, tailer,0.9)
    interp_z=min(vertex1_z*0.4463 , (interp_y-708.937)*0.5656)
    wing_leading_y=1700
    wing_leading_x=11957.5
    wing_leading_z=18.379
    vertex2_z=1200
    vertex1=[vertex1_x, 0.0, vertex1_z]  # the leading_point of the whole aircraft /*x:Parameter  y:Fixed  z:Parameter*/
    inlet_edge=[6400 , 705.937 , 0] # /*x:Fixed  y:Fixed  z:Fixed*/
    inlet_edge[1]+=3.0
    wing_leading=[wing_leading_x, wing_leading_y, wing_leading_z] # /*x:Parameter  y:Fixed  z:Parameter*/
    vertex2=[wing_leading_x, 0.0 ,vertex2_z]  # the tailing_point of the top_line of the body1 /*x:Fixed  y:Fixed  z:Parameter*/
    wing_tailing_up=[21312.162, wing_leading_y, wing_leading_z-45.838]
    wing_tailing_low=[21312.162, wing_leading_y, wing_leading_z-48.838]
    #------------------
    tail_A4_y=825.447
    tail_A4_z=wing_tailing_up[2]-3
    tail_A3_z=tail_A4_z+129.209
    tail_B1_z=tail_A4_z-173.56
    tail_A3_x=22701.6+1.467*tail_A3_z
    tail_A4_x=tail_A3_x+389.286
    tail_A5_x=tail_A3_x-277.673
    #------------------
    tail_A=[[24011.355, 0, 895.581], [24011.355, 679.791, 895.581], [23965.349, 702.993, 861.438],\
        [tail_A3_x,702.993,tail_A3_z], [tail_A4_x, tail_A4_y, tail_A4_z], [tail_A5_x, wing_leading_y, tail_A4_z]]
    tail_B=[[22701.6+1.467*tail_B1_z, 0, tail_B1_z], [22701.6+1.467*tail_B1_z, 702.993, tail_B1_z], \
        [tail_A4_x, tail_A4_y, tail_A[4][2]-3], [tail_A5_x, wing_leading_y, tail_A[4][2]-3]]
#* ============================================
#* support
#* ============================================
    if True:
        support = OpenSurface(n_sec=3, name='support', nn=601, ns=general_ns)
        support.read_setting('./files/fuselage.txt')
        ns=support.ns
        nn=support.nn
        # define the top_line (support)
        if True:
            support.secs[0].xLE=vertex1[0]
            support.secs[0].yLE=vertex1[1]
            support.secs[0].zLE=vertex1[2]+2.72
            support.secs[0].xTE=vertex2[0]
            support.secs[0].yTE=vertex2[1]
            support.secs[0].zTE=vertex2[2]
            support.secs[0].section(nn=nn)
            scalez=support.secs[0].zTE-support.secs[0].zLE
            scalex=support.secs[0].xTE-support.secs[0].xLE
            b1=scalex/scalez*math.tan(slope_top_forebody/180*pi)  # decided the slope of the top_line(leading edge)
            for j in range(nn):
                tt=1.0*j/(1.0*(nn-1))
                support.secs[0].x[j]=(1-tt)*support.secs[0].xLE+tt*support.secs[0].xTE
                support.secs[0].y[j]=0
                support.secs[0].z[j]=(1-(1-tt)**b1)*scalez+support.secs[0].zLE
        # define the side_line (support)
        if True:
            support.secs[1].xLE=vertex1[0]
            support.secs[1].yLE=vertex1[1]
            support.secs[1].zLE=vertex1[2]
            support.secs[1].xTE=wing_leading[0]
            support.secs[1].yTE=wing_leading[1]
            support.secs[1].zTE=wing_leading[2]
            scalex=support.secs[1].xTE-support.secs[1].xLE
            scaley=support.secs[1].yTE-support.secs[1].yLE
            scalez=support.secs[1].zTE-support.secs[1].zLE
            k1y=math.tan(slope_side_forebody/180*math.pi)*scalex/scaley    # decide the slope of the leading_edge in xy-plane, parameter
            k2y=0.2232*scalex/scaley  # decide the slope of the tailing_edge, which connect to the wing in xy-plane, fixed
            k1z=0.0 # decide the slope of the leading_edge in yz-plane, parameter,
            k2z=-0.0151*scalex/scalez # decide the slope of the leading_edge in yz-plane, fixed
            support.secs[1].section(nn=nn)
            interp_norm_x=(interp_x-support.secs[1].xLE)/scalex
            interp_norm_y=interp_y/scaley
            interp_norm_z=(interp_z-support.secs[1].zLE)/scalez
            csy = CubicSpline([0, interp_norm_x, 1], [0, interp_norm_y, 1], bc_type=((1, k1y), (1, k2y)))
            csz = CubicSpline([0, interp_norm_x, 1], [0, interp_norm_z, 1], bc_type=((2, k1z), (1, k2z)))
            for j in range(nn):
                tt = 1.0*j/(nn-1)
                support.secs[1].x[j]=(1-tt)*support.secs[1].xLE+tt*support.secs[1].xTE
                support.secs[1].z[j]=(1-tt)*support.secs[1].zLE+tt*support.secs[1].zTE
                support.secs[1].y[j]=csy(tt)*scaley
                support.secs[1].z[j]=csz(tt)*scalez+support.secs[1].zLE
        # define the lower_line (support)
        if True:
            support.secs[2].xLE=inlet_edge[0]
            support.secs[2].yLE=0.0
            support.secs[2].zLE=inlet_edge[2]
            support.secs[2].xTE=22701.6+1.467*tail_B1_z
            support.secs[2].yTE=0.0
            support.secs[2].zTE=tail_B1_z
            scalex=support.secs[2].xTE-support.secs[2].xLE
            scalez=support.secs[2].zTE-support.secs[2].zLE
            k1=min(math.tan(slope_leading_tank/180*math.pi) , 0.2448)*scalex/scalez    # decide the slope of the leading_edge, parameter
            k2=math.tan(slope_ending_tank/180*math.pi)*scalex/scalez    # decide the slope of the leading_edge, parameter
            interp1=[10886.169452089,-249.546254846]
            interp2=[16969.96687, -657.1929723]
            interp1_norm=[(interp1[0]-support.secs[2].xLE)/scalex, (interp1[1]-support.secs[2].zLE)/scalez]
            interp2_norm=[(interp2[0]-support.secs[2].xLE)/scalex, (interp2[1]-support.secs[2].zLE)/scalez]
            cs = CubicSpline([0, interp1_norm[0], interp2_norm[0], 1], [0, interp1_norm[1], interp2_norm[1],  1], bc_type=((1, k1), (1, k2)))
            support.secs[2].section(nn=nn)
            for j in range(nn):
                tt = 1.0*j/(nn-1)
                support.secs[2].x[j]=(1-tt)*support.secs[2].xLE+tt*support.secs[2].xTE
                support.secs[2].y[j]=0.0
                support.secs[2].z[j]=support.secs[2].zLE+cs(tt)*scalez

        # interpolation between the edge
        support.geo(update_sec=False)     
#* ============================================
#* build the main fuselage, by adjusting 'a_lower' to match 'ideal_volume', the bigger, the fuller tank is, so that increasing the volume total
#* ============================================
    if True:
        forebody = OpenSurface(n_sec=4, name='forebody', nn=301, ns=general_ns)
        forebody.read_setting('./files/fuselage.txt')
        body1 = OpenSurface(n_sec=4, name='body1', nn=301, ns=general_ns)
        body1.read_setting('./files/fuselage.txt')
        uppersurf = OpenSurface(n_sec=2, name='uppersurf', nn=201, ns=general_ns)
        uppersurf.read_setting('./files/fuselage.txt')
        lowersurf = OpenSurface(n_sec=2, name='lowersurf', nn=201, ns=general_ns)
        lowersurf.read_setting('./files/fuselage.txt')
        ending_upper = OpenSurface(n_sec=2, name='ending_upper', nn=101, ns=(general_ns-1)*16+1)
        ending_upper.read_setting('./files/fuselage.txt')
        ending_lower = OpenSurface(n_sec=2, name='ending_lower', nn=101, ns=(general_ns-1)*16+1)
        ending_lower.read_setting('./files/fuselage.txt')
        forebody, body1, uppersurf, lowersurf, ending_upper, ending_lower = build()
        index0=int(np.floor(tail_B[1][1]/ending_lower.secs[1].yTE*ns))
        volume_total=volume_caculation()
        delta=abs(volume_total-ideal_volume)
        print('************Volume Modify************')
        print(volume_total)
        error=0.001
        if delta >= error:
            step=50
            ratio=vertex2_z/upper_h2
            while delta >= error:
                delta_temp=delta
                volume_temp=volume_total
                if volume_total-ideal_volume > 0.0:
                    upper_h1-=step
                    upper_h2-=step
                    upper_h1m-=step/3
                    upper_h2m-=step/3
                    vertex2_z-=step*ratio
                else:
                    upper_h1+=step
                    upper_h2+=step
                    upper_h1m+=step/3
                    upper_h2m+=step/3
                    vertex2_z+=step*ratio
                forebody, body1, uppersurf, lowersurf, ending_upper, ending_lower = build()
                index0=int(np.floor(tail_B[1][1]/ending_lower.secs[1].yTE*ns))
                volume_total=volume_caculation()
                delta=abs(volume_total-volume_temp)
                if delta == delta_temp:
                    step=step/2
                print(volume_total)
#* ============================================
#* side
#* ============================================
    if True:
        side = OpenSurface(n_sec=2, name='side', nn=201, ns=general_ns)
        side.read_setting('./files/fuselage.txt')
        ns=side.ns
        nn=side.nn
        # produce and specify the shape of secs[0], which is connect to the tank
        if True:
            side.secs[0].section(nn=nn)
            side.secs[0].x[0]=6400
            side.secs[0].y[0]=705.937
            side.secs[0].z[0]=0
            side.secs[0].x[-1]=22701.6+1.467*tail_B1_z
            side.secs[0].y[-1]=702.993
            side.secs[0].z[-1]=tail_B1_z
            for j in range(1,nn-1):
                tt = 1.0*j/(nn-1)
                side.secs[0].x[j] = (1-tt)*side.secs[0].x[0]+tt*side.secs[0].x[-1]
                side.secs[0].y[j] = (1-tt)*side.secs[0].y[0]+tt*side.secs[0].y[-1]
            curv = CubicSpline(support.secs[2].x[:], support.secs[2].z[:])
            side.secs[0].z[1:nn] = curv(side.secs[0].x[1:nn])
        # produce and specify the shape of secs[1], which is the side section which is the lower curve
        if True:
            side.secs[1].section(nn=nn)
            side.secs[1].x[0]=11408.557
            side.secs[1].y[0]=705.969
            side.secs[1].z[0]=-1245.914          
            side.secs[1].x[-1]=20678.591
            side.secs[1].y[-1]=702.993
            side.secs[1].z[-1]=-1378.991
            for j in range(1,nn-1):
                tt = 1.0*j/(nn-1)
                side.secs[1].x[j] = (1-tt)*side.secs[1].x[0]+tt*side.secs[1].x[-1]
                side.secs[1].y[j] = (1-tt)*side.secs[1].y[0]+tt*side.secs[1].y[-1]
                side.secs[1].z[j] = (1-tt)*side.secs[1].z[0]+tt*side.secs[1].z[-1]
        # produce the geometry & obtain the number of secs between secs[1] and secs[2]
        side.geo(update_sec=False)
        # define the shape of the baseline 
        if True:
            # front side
            index=int(np.floor(-1221.009/side.secs[1].z[-1]*ns))
            for j in range(1,index+1):
                tt=1.0*j/index
                side.surfs[0][0][j,0]=(1-tt)*6400+tt*11388.739
                side.surfs[0][1][j,0]=(1-tt)*705.937+tt*705.969
                side.surfs[0][2][j,0]=(1-tt)*0.0+tt*(-1221.009)
            for j in range(index+1,ns-1):
                tt = 1.0*(j-index)/(ns-index-1)
                side.surfs[0][0][j,0] =  -6.308*tt**3 - 8.0534*tt**2 + 34.178*tt + 11389
                side.surfs[0][1][j,0] = 705.969
                side.surfs[0][2][j,0] = 8.0502*tt**3 - 26.024*tt**2 - 6.9313*tt - 1221
           # back side
            for j in range(1,ns-1):
                tt = 1.0*j/(ns-1)
                side.surfs[0][0][j,-1] = (1-tt)*side.surfs[0][0][0,-1]+tt*side.surfs[0][0][-1,-1]
                side.surfs[0][1][j,-1] = (1-tt)*side.surfs[0][1][0,-1]+tt*side.surfs[0][1][-1,-1]
                side.surfs[0][2][j,-1] = (1-tt)*side.surfs[0][2][0,-1]+tt*side.surfs[0][2][-1,-1]

        # interpolation between the edge
        surfx=[[]for i in range(nn)]
        surfy=[[]for i in range(nn)]
        surfz=[[]for i in range(nn)]
        scalex=np.zeros(nn)
        scaley=np.zeros(nn)
        scalez=np.zeros(nn)
        for i_surf in range(1):
            for i in range(1,nn-1):
                tt = 1.0*i/(nn-1)
                surfx[i] = tt*side.surfs[i_surf][0][:,-1]+(1-tt)*side.surfs[i_surf][0][:,0]
                surfy[i] = tt*side.surfs[i_surf][1][:,-1]+(1-tt)*side.surfs[i_surf][1][:,0]
                surfz[i] = tt*side.surfs[i_surf][2][:,-1]+(1-tt)*side.surfs[i_surf][2][:,0]
                scalex[i] = (side.surfs[i_surf][0][-1,i]-side.surfs[i_surf][0][0,i])/(surfx[i][-1]-surfx[i][0])
                scaley[i] = (side.surfs[i_surf][1][-1,i]-side.surfs[i_surf][1][0,i])/(surfy[i][-1]-surfy[i][0])
                scalez[i] = (side.surfs[i_surf][2][-1,i]-side.surfs[i_surf][2][0,i])/(surfz[i][-1]-surfz[i][0])
                side.surfs[i_surf][0][1:ns-1,i] = side.surfs[i_surf][0][0,i]+(surfx[i][1:ns-1]-surfx[i][0])*scalex[i]
                side.surfs[i_surf][1][1:ns-1,i] = side.surfs[i_surf][1][0,i]+(surfy[i][1:ns-1]-surfy[i][0])*scaley[i]
                side.surfs[i_surf][2][1:ns-1,i] = side.surfs[i_surf][2][0,i]+(surfz[i][1:ns-1]-surfz[i][0])*scalez[i]
#* ============================================
#* display the parameter
#* ============================================
    if True:
        print('************Volume Constrain************')
        print('Height of the Fuselage = ',vertex2_z)
        print('model_volume = ', volume_caculation())
        print('ideal_volume = ', ideal_volume)
        print('Height of the Fuselage = ',vertex2_z)
        print('************Design Variable************')
        print('forebody_length = ',forebody_length)
        print('forebody_thickness = ',forebody_thickness)
        print('forebody_width = ',forebody_width)
        print('nose_height = ',nose_height)
        print('nose_upper_angle = ',nose_upper_angle)
        print('nose_side_angle = ',nose_side_angle)
        print('radius_LE = ',radius_LE)
        print('ideal_volume = ',ideal_volume)
#* ============================================
#* output the result
#* ============================================
    if True:
        output()
