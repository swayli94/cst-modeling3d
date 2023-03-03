'''
Jet engine by Long Yizhu (2023.01)
'''

import numpy as np
import pandas as pd
import math

from cst_modeling.basic import BasicSection, BasicSurface
from cst_modeling.section import cst_curve 


class Section(BasicSection):

    def __init__(self, n_curve = 4):
        self.x = np.zeros(1)
        self.y = np.zeros(1)
        self.z = np.zeros(1)
        self.n_ = max(1,n_curve)
        self.curve = [curve() for _ in range(self.n_)]

    def set_Section(self, sec_x, sec_y, sec_z):
        self.x = sec_x
        self.y = sec_y
        self.z = sec_z

    def join_curve(self):
        for i in range(self.n_):
            self.x = np.append(self.x, self.curve[i].curve_x)
            self.y = np.append(self.y, self.curve[i].curve_y)
            self.z = np.append(self.z, self.curve[i].curve_z)

        self.x = np.delete(self.x,0)
        self.y = np.delete(self.y,0)
        self.z = np.delete(self.z,0)

        self.x,self.y,self.z = order_curve3d(self.x,self.y,self.z)

        self.x[0] = self.x[-1]
        self.y[0] = self.y[-1]
        self.z[0] = self.z[-1]

    def circle(self, p1=np.array([1,0]), p2=np.array([0,1]), p3=np.array([0,-1]), nn=500):
        A = p1[0]*(p2[1]-p3[1]) - p1[1]*(p2[0]-p3[0]) + p2[0]*p3[1] - p3[0]*p2[1]
        if np.abs(A) <= 1E-20:
            raise Exception('Finding circle: 3 points in one line')
        
        p1s = p1[0]**2 + p1[1]**2
        p2s = p2[0]**2 + p2[1]**2
        p3s = p3[0]**2 + p3[1]**2

        B = p1s*(p3[1]-p2[1]) + p2s*(p1[1]-p3[1]) + p3s*(p2[1]-p1[1])
        C = p1s*(p2[0]-p3[0]) + p2s*(p3[0]-p1[0]) + p3s*(p1[0]-p2[0])
        D = p1s*(p3[0]*p2[1]-p2[0]*p3[1]) + p2s*(p1[0]*p3[1]-p3[0]*p1[1]) + p3s*(p2[0]*p1[1]-p1[0]*p2[1])

        cen_x = -B/2/A
        cen_y = -C/2/A

        R = np.sqrt(B**2+C**2-4*A*D)/2/np.abs(A)

        for i in range(nn):
            theta = i*2*np.pi/(nn-1)
            self.x = np.append(self.x, cen_x + R*math.cos(theta))
            self.y = np.append(self.y, cen_y + R*math.sin(theta))
        self.x = np.delete(self.x,0)
        self.y = np.delete(self.y,0)


    def order_curve3d(x,y,z,axis = 'z'):
        if axis == 'x':
            x1 = y
            x2 = z
            x3 = x
        if axis == 'y':
            x1 = x
            x2 = z
            x3 = y
        if axis == 'z':  
            x1 = x
            x2 = y
            x3 = z

        cen_x = np.mean(x1)
        cen_y = np.mean(x2)
        
        x1s = []
        x2s = []
        x3s = []
        for i in range(len(x1)):
            dx = x1[i] - cen_x
            dy = x2[i] - cen_y
            ang = np.arctan2(dy,dx)
            if ang < 0:
                ang += np.pi*2
            x1s.append([x1[i],ang])
            x2s.append([x2[i],ang])
            x3s.append([x3[i],ang])

        x1s = sorted(x1s, key = lambda x:x[1])
        x2s = sorted(x2s, key = lambda x:x[1])
        x3s = sorted(x3s, key = lambda x:x[1])

        x1 = np.array([x[0] for x in x1s])
        x2 = np.array([x[0] for x in x2s])
        x3 = np.array([x[0] for x in x3s])

        if axis == 'x':
            return x3,x1,x2
        if axis == 'y':
            return x1,x3,x2
        if axis == 'z':  
            return x1,x2,x3

class curve():

    def __init__(self, r0=np.array([0,0,0]), r1=np.array([1,0,0]), normal = np.array([0,0,1]), nn=125, coef = np.array([0,0,0,0]), xn1=0.75, xn2=1.0):
        #确定曲线增加的方向
        dr = r1 - r0
        D = np.cross(normal,dr)
        D = D / np.linalg.norm(D)

        #确定基本直线
        L0 = np.linalg.norm(dr)
        base_x = np.zeros(nn)
        base_y = np.zeros(nn)
        base_z = np.zeros(nn)
        for i in range(nn):
            tt = i / (nn-1)
            base_x[i] = (1-tt)*r0[0] + tt*r1[0]
            base_y[i] = (1-tt)*r0[1] + tt*r1[1]
            base_z[i] = (1-tt)*r0[2] + tt*r1[2]

        #生成cst曲线
        cst_x,cst_y = cst_curve(nn,coef,None,xn1,xn2)
        cst_y = cst_y * L0

        self.curve_x = base_x + D[0]*cst_y
        self.curve_y = base_y + D[1]*cst_y
        self.curve_z = base_z + D[2]*cst_y

def read_sec_file(filename,scale = 1):
    
    df = pd.read_excel(filename,header=0)
    data = df.values
    z = data[1:,1]
    y = data[1:,2]
    x = data[1:,3]

    x,y,z = order_curve3d(x,y,z)

    x = np.append(x,x[0])
    y = np.append(y,y[0])
    z = np.append(z,z[0])

    return x*scale,y*scale,z*scale

def order_curve3d(x,y,z,axis = 'z'):
    
    if axis == 'x':
        x1 = y
        x2 = z
        x3 = x
    if axis == 'y':
        x1 = x
        x2 = z
        x3 = y
    if axis == 'z':  
        x1 = x
        x2 = y
        x3 = z

    cen_x = np.mean(x1)
    cen_y = np.mean(x2)
    
    x1s = []
    x2s = []
    x3s = []
    for i in range(len(x1)):
        dx = x1[i] - cen_x
        dy = x2[i] - cen_y
        ang = np.arctan2(dy,dx)
        if ang < 0:
            ang += np.pi*2
        x1s.append([x1[i],ang])
        x2s.append([x2[i],ang])
        x3s.append([x3[i],ang])

    x1s = sorted(x1s, key = lambda x:x[1])
    x2s = sorted(x2s, key = lambda x:x[1])
    x3s = sorted(x3s, key = lambda x:x[1])

    x1 = np.array([x[0] for x in x1s])
    x2 = np.array([x[0] for x in x2s])
    x3 = np.array([x[0] for x in x3s])

    if axis == 'x':
        return x3,x1,x2
    if axis == 'y':
        return x1,x3,x2
    if axis == 'z':  
        return x1,x2,x3


if __name__== "__main__":

    #create outward
    n_sec = 5

    outward = BasicSurface(n_sec=5, name='out', nn=201, ns=51, projection = False)

    outward.secs[0] = Section(n_curve = 4)
    outward.secs[0].curve[0] = curve(r0=np.array([0,     -1.505,15    ]), r1=np.array([0,     -3.030,15    ]), nn=125)
    outward.secs[0].curve[1] = curve(r0=np.array([0,     -3.030,15    ]), r1=np.array([-1.395,-3.025,18.532]), nn=125)
    outward.secs[0].curve[2] = curve(r0=np.array([-1.395,-3.025,18.532]), r1=np.array([-1.395,-1.505,18.532]), nn=125)
    outward.secs[0].curve[3] = curve(r0=np.array([-1.395,-1.505,18.532]), r1=np.array([0,     -1.505,15    ]), nn=125)
    outward.secs[0].join_curve()

    outward.secs[-1] = Section(n_curve = 4)
    outward.secs[-1].curve[0] = curve(r0=np.array([0,     -1.765,29.570]), r1=np.array([0,     -2.765,29.570]), nn=125)
    outward.secs[-1].curve[1] = curve(r0=np.array([0,     -2.765,29.570]), r1=np.array([-1.445,-2.765,27.068]), nn=125)
    outward.secs[-1].curve[2] = curve(r0=np.array([-1.445,-2.765,27.068]), r1=np.array([-1.445,-1.765,27.068]), nn=125)
    outward.secs[-1].curve[3] = curve(r0=np.array([-1.445,-1.765,27.068]), r1=np.array([0,     -1.765,29.570]), nn=125)
    outward.secs[-1].join_curve()

    outward.secs[1] = Section(n_curve = 4)
    outward.secs[1].curve[0] = curve(r0=np.array([0,-1.50,20.307]), r1=np.array([0,-3.030,20.307]), normal=np.array([0,0,1]), nn=125, xn1=0.5, xn2=0.5, 
                                     coef=1.0*np.array([0.22706027, 1.1, 0.19973993, 0.9007924,  0.44730042, 0.44728627, 0.90080244, 0.19973533, 1.1,  0.22706104]))
    outward.secs[1].curve[1] = curve(r0=np.array([0,     -3.030,20.307]), r1=np.array([-1.676,-3.025,20.307]), nn=125)
    outward.secs[1].curve[2] = curve(r0=np.array([-1.676,-3.025,20.307]), r1=np.array([-1.676,-1.505,20.307]), nn=125)
    outward.secs[1].curve[3] = curve(r0=np.array([-1.676,-1.505,20.307]), r1=np.array([0,     -1.50, 20.307]), nn=125)
    outward.secs[1].join_curve()

    outward.secs[2] = Section(n_curve = 4)
    outward.secs[2].curve[0] = curve(r0=np.array([0,-1.50,23.307]), r1=np.array([0,-3.030,23.307]), normal=np.array([0,0,1]), nn=125, xn1=0.5, xn2=0.5,
                                     coef=1.2*np.array([0.22706027, 1.1, 0.19973993, 0.9007924,  0.44730042, 0.44728627, 0.90080244, 0.19973533, 1.1,  0.22706104]))
    outward.secs[2].curve[1] = curve(r0=np.array([0,     -3.030,23.307]), r1=np.array([-1.573,-2.910,23.307]), nn=125)
    outward.secs[2].curve[2] = curve(r0=np.array([-1.573,-2.910,23.307]), r1=np.array([-1.573,-1.620,23.307]), nn=125)
    outward.secs[2].curve[3] = curve(r0=np.array([-1.573,-1.620,23.307]), r1=np.array([0,     -1.50, 23.307]), nn=125)
    outward.secs[2].join_curve()

    outward.secs[3] = Section(n_curve = 4)
    outward.secs[3].curve[0] = curve(r0=np.array([0,-1.50,26.197]), r1=np.array([0,-3.030,26.197]), normal=np.array([0,0,1]), nn=125, xn1=0.5, xn2=0.5, 
                                     coef=0.8*np.array([0.22706027, 1.1, 0.19973993, 0.9007924,  0.44730042, 0.44728627, 0.90080244, 0.19973533, 1.1,  0.22706104]))
    outward.secs[3].curve[1] = curve(r0=np.array([0,     -3.030,26.197]), r1=np.array([-1.474,-2.789,26.197]), nn=125)
    outward.secs[3].curve[2] = curve(r0=np.array([-1.474,-2.789,26.197]), r1=np.array([-1.474,-1.732,26.197]), nn=125)
    outward.secs[3].curve[3] = curve(r0=np.array([-1.474,-1.732,26.197]), r1=np.array([0,     -1.50, 26.197]), nn=125)
    outward.secs[3].join_curve()

    outward.geo(update_sec = False)

    outward.smooth(0,n_sec-1,smooth0=True,ratio_end=-1)

    outward.translate(dX=0, dY=4.53)

    outward.flip(axis = '+Y')
    outward.flip(axis = '+X')

    outward.output_tecplot(fname = './dump/engine_outward.dat')


    #create inward
    n_sec = 9

    inward = BasicSurface(n_sec=9, name='out', nn=201, ns=51, projection = False)

    inward.secs[0] = Section(n_curve = 4)
    inward.secs[0].curve[0] = curve(r0=np.array([0,     -1.505+4.53,15    ]), r1=np.array([0,     -3.030+4.53,15    ]), nn=125)
    inward.secs[0].curve[1] = curve(r0=np.array([0,     -3.030+4.53,15    ]), r1=np.array([-1.395,-3.025+4.53,18.532]), nn=125)
    inward.secs[0].curve[2] = curve(r0=np.array([-1.395,-3.025+4.53,18.532]), r1=np.array([-1.395,-1.505+4.53,18.532]), nn=125)
    inward.secs[0].curve[3] = curve(r0=np.array([-1.395,-1.505+4.53,18.532]), r1=np.array([0,     -1.505+4.53,15    ]), nn=125)
    inward.secs[0].join_curve()

    inward.secs[1] = Section(n_curve = 4)
    inward.secs[1].curve[0] = curve(r0=np.array([-0.314,3.025,17.006]), r1=np.array([-0.314,1.505,17.006]), nn=125)
    inward.secs[1].curve[1] = curve(r0=np.array([-0.314,1.505,17.006]), r1=np.array([-1.549,1.505,19.502]), nn=125)
    inward.secs[1].curve[2] = curve(r0=np.array([-1.549,1.505,19.502]), r1=np.array([-1.549,3.025,19.502]), nn=125)
    inward.secs[1].curve[3] = curve(r0=np.array([-1.549,3.025,19.502]), r1=np.array([-0.314,3.025,17.006]), nn=125)
    inward.secs[1].join_curve()

    inward.secs[2] = Section(n_curve = 4)
    inward.secs[2].curve[0] = curve(r0=np.array([-1.174,3.025,19.368]), r1=np.array([-1.174,1.505,19.368]), nn=125)
    inward.secs[2].curve[1] = curve(r0=np.array([-1.174,1.505,19.368]), r1=np.array([-1.580,1.505,19.904]), nn=125)
    inward.secs[2].curve[2] = curve(r0=np.array([-1.580,1.505,19.904]), r1=np.array([-1.580,3.025,19.904]), nn=125)
    inward.secs[2].curve[3] = curve(r0=np.array([-1.580,3.025,19.904]), r1=np.array([-1.174,3.025,19.368]), nn=125)
    inward.secs[2].join_curve()

    inward.secs[3] = Section(n_curve = 4)
    inward.secs[3].curve[0] = curve(r0=np.array([-1.340,3.025,20.307]), r1=np.array([-1.340,1.505,20.307]), nn=125)
    inward.secs[3].curve[1] = curve(r0=np.array([-1.340,1.505,20.307]), r1=np.array([-1.591,1.505,20.307]), nn=125)
    inward.secs[3].curve[2] = curve(r0=np.array([-1.591,1.505,20.307]), r1=np.array([-1.591,3.025,20.307]), nn=125)
    inward.secs[3].curve[3] = curve(r0=np.array([-1.591,3.025,20.307]), r1=np.array([-1.340,3.025,20.307]), nn=125)
    inward.secs[3].join_curve()

    inward.secs[4] = Section(n_curve = 4)
    inward.secs[4].circle(np.array([-1.025,2.507]),np.array([-1.025,2.023]),np.array([-1.305,2.023]), nn=500)
    inward.secs[4].z = 23.000 * np.ones_like(inward.secs[4].x)
    #将曲线排序
    inward.secs[4].x,inward.secs[4].y,inward.secs[4].z = order_curve3d(inward.secs[4].x,inward.secs[4].y,inward.secs[4].z)
    #让曲线封闭
    inward.secs[4].x[0] = inward.secs[4].x[-1]
    inward.secs[4].y[0] = inward.secs[4].y[-1]
    inward.secs[4].z[0] = inward.secs[4].z[-1]

    inward.secs[5] = Section(n_curve = 4)
    inward.secs[5].circle(np.array([-1.025,2.507]),np.array([-1.025,2.023]),np.array([-1.305,2.023]), nn=500)
    inward.secs[5].z = 24.907 * np.ones_like(inward.secs[5].x)
    #将曲线排序
    inward.secs[5].x,inward.secs[5].y,inward.secs[5].z = order_curve3d(inward.secs[5].x,inward.secs[5].y,inward.secs[5].z)
    #让曲线封闭
    inward.secs[5].x[0] = inward.secs[5].x[-1]
    inward.secs[5].y[0] = inward.secs[5].y[-1]
    inward.secs[5].z[0] = inward.secs[5].z[-1]

    #sec_x, sec_y, sec_z = read_sec_file('inward-sec4.xlsx',0.001)
    inward.secs[6] = Section(n_curve = 4)
    inward.secs[6].curve[0] = curve(r0=np.array([-1.045,2.765,25.670]), r1=np.array([-1.045,1.765,25.670]), nn=125)
    inward.secs[6].curve[1] = curve(r0=np.array([-1.045,1.765,25.670]), r1=np.array([-1.323,1.765,25.670]), nn=125)
    inward.secs[6].curve[2] = curve(r0=np.array([-1.323,1.765,25.670]), r1=np.array([-1.323,2.765,25.670]), nn=125)
    inward.secs[6].curve[3] = curve(r0=np.array([-1.323,2.765,25.670]), r1=np.array([-1.045,2.765,25.670]), nn=125)
    inward.secs[6].join_curve()

    #inlet_x,inlet_y,inlet_z = read_sec_file('inward-outlet1.xlsx',0.001)
    inward.secs[7] = Section(n_curve = 4)
    inward.secs[7].curve[0] = curve(r0=np.array([-0.361,2.765,27.984]), r1=np.array([-0.361,1.765,27.984]), nn=125)
    inward.secs[7].curve[1] = curve(r0=np.array([-0.361,1.765,27.984]), r1=np.array([-1.372,1.765,26.234]), nn=125)
    inward.secs[7].curve[2] = curve(r0=np.array([-1.372,1.765,26.234]), r1=np.array([-1.372,2.765,26.234]), nn=125)
    inward.secs[7].curve[3] = curve(r0=np.array([-1.372,2.765,26.234]), r1=np.array([-0.361,2.765,27.984]), nn=125)
    inward.secs[7].join_curve()

    #sec_x, sec_y, sec_z = read_sec_file('inward-outlet2.xlsx',0.001)
    inward.secs[8] = Section(n_curve = 4)
    inward.secs[8].curve[0] = curve(r0=np.array([0,     -1.765+4.53,29.570]), r1=np.array([0,     -2.765+4.53,29.570]), nn=125)
    inward.secs[8].curve[1] = curve(r0=np.array([0,     -2.765+4.53,29.570]), r1=np.array([-1.445,-2.765+4.53,27.068]), nn=125)
    inward.secs[8].curve[2] = curve(r0=np.array([-1.445,-2.765+4.53,27.068]), r1=np.array([-1.445,-1.765+4.53,27.068]), nn=125)
    inward.secs[8].curve[3] = curve(r0=np.array([-1.445,-1.765+4.53,27.068]), r1=np.array([0,     -1.765+4.53,29.570]), nn=125)
    inward.secs[8].join_curve()

    inward.geo(update_sec = False)

    inward.smooth(2,4,smooth0=True,ratio_end=-1)
    inward.smooth(5,6,smooth0=True,ratio_end=-1)
    inward.smooth(6,8,smooth1=True,ratio_end=-1)

    inward.flip(axis = '+Y')
    inward.flip(axis = '+X')

    inward.output_tecplot(fname = './dump/engine_inward.dat')
