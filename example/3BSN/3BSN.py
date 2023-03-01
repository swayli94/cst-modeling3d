
import numpy as np

from cst_modeling.basic import output_curve, rotation_3d
from cst_modeling.surface import BasicSurface

from scipy.interpolate import CubicSpline

import os


if __name__ == "__main__":


    # the length of the 4 sections (3 rotation section and one convergent section)
    lsec = (0.6, 0.8, 0.5, 0.6)
    # number of cross-sections in each section
    nzsec = (2, 2, 2, 2)
    # the radius of the cross-sections bet. each sections
    rsec = (0.6, 0.7, 0.7, 0.6, 0.5)

    max_theta = 95                      # maximum deflection angle of the nozzle

    sec_theta = max_theta / 4.
    gamma = sec_theta / 180 * np.pi

    for irot in range(75):

        rotate_angle = 90 - 90 / 75 * irot
        # axises for section 2 and section 3
        axis2 =  np.array([[0.0, 0.0, lsec[0]], [0.0, -np.sin(gamma), lsec[0]+np.cos(gamma)]])
        axis3 =  np.array([[0.0, 0.0, lsec[0]+lsec[1]], [0.0, np.sin(gamma), lsec[0]+lsec[1]+np.cos(gamma)]])

        #* ============================================
        #* Section 1
        #* ============================================

        sec1 = BasicSurface(n_sec=nzsec[0], name='Sec1', nn=101, ns=31, project=False)
        aa = np.linspace(0.0, 2*np.pi, num=sec1.nn)
        
        for i in range(nzsec[0]):
            sec1.secs[i].set_params(xLE=0.0, yLE=0.0, zLE=lsec[0] / (nzsec[0] - 1) * i)
            sec1.secs[i].xx = np.cos(aa) * (rsec[0] + (rsec[1] - rsec[0]) / (nzsec[0] - 1) * i)
            sec1.secs[i].yy = np.sin(aa) * (rsec[0] + (rsec[1] - rsec[0]) / (nzsec[0] - 1) * i)

        sec1.geo_secs()
        
        for i in range(nzsec[0]):
            sec1.secs[i].rotate(angle=sec_theta / (nzsec[0] - 1) * i, origin=[0.0, 0.0, lsec[0] / (nzsec[0] - 1) * i], axis='X')

        sec1.geo(update_sec=False)
        sec1.rotate(origin=np.array([0,0,0]), axis=np.array([0,0,1]), angle=rotate_angle)
        sec1.flip(axis='+Y')

        sec1.output_tecplot(fname='plot/sec1.dat')

        #* ============================================
        #* Section 2
        #* ============================================

        sec2 = BasicSurface(n_sec=nzsec[1], name='Sec2', nn=101, ns=31, project=False)
        aa = np.linspace(0.0, 2*np.pi, num=sec2.nn)

        for i in range(nzsec[1]):
            sec2.secs[i].set_params(xLE=0.0, yLE=0.0, zLE=lsec[0] + lsec[1] / (nzsec[1] - 1) * i)
            sec2.secs[i].xx = np.cos(aa) * (rsec[1] + (rsec[2] - rsec[1]) / (nzsec[1] - 1) * i)
            sec2.secs[i].yy = np.sin(aa) * (rsec[1] + (rsec[2] - rsec[1]) / (nzsec[1] - 1) * i)
        
        sec2.geo_secs()

        for i in range(nzsec[1]):
            sec2.secs[i].rotate(angle=sec_theta - 2 *  sec_theta/ (nzsec[1] - 1) * i, 
                                origin=[0.0, 0.0, lsec[0] + lsec[1] / (nzsec[1] - 1) * i], axis='X')
            
        sec2.geo(update_sec=False)

        sec2.rotate(origin=np.array([0,0,0]), axis=np.array([0,0,1]), angle=rotate_angle)
        axis2 = rotation_3d(axis2, origin=np.array([0,0,0]), axis=np.array([0,0,1]), angle=rotate_angle)
        sec2.rotate(origin=axis2[0], axis=axis2[1] - axis2[0], angle=-2 * rotate_angle)
        sec2.flip(axis='+Y')

        sec2.output_tecplot(fname='plot/sec2.dat')

        #* ============================================
        #* Section 3 and Covergent nozzle
        #* ============================================

        sec3 = BasicSurface(n_sec=sum(nzsec[2:])-len(nzsec)+3, name='Sec3', nn=101, ns=21, project=False)
        aa = np.linspace(0.0, 2*np.pi, num=sec3.nn)

        for isec in range(2, 4):
            for i in range(nzsec[isec]):
                sec_no = sum(nzsec[2:isec]) + i - (isec - 2)
                sec3.secs[sec_no].set_params(xLE=0.0, yLE=0.0, zLE=sum(lsec[0:isec]) + lsec[isec] / (nzsec[isec] - 1) * i)
                sec3.secs[sec_no].xx = np.cos(aa) * (rsec[isec] + (rsec[isec+1] - rsec[isec]) / (nzsec[isec] - 1) * i)
                sec3.secs[sec_no].yy = np.sin(aa) * (rsec[isec] + (rsec[isec+1] - rsec[isec]) / (nzsec[isec] - 1) * i)

        sec3.geo_secs()

        for i in range(nzsec[2]):
            sec3.secs[i].rotate(angle=-sec_theta + sec_theta / (nzsec[1] - 1) * i, 
                                origin=[0.0, 0.0, sum(lsec[0:2]) + lsec[2] / (nzsec[2] - 1) * i], axis='X')

        sec3.geo(update_sec=False)
        # sec3.smooth(1, 4, smooth0=False, smooth1=False)

        sec3.rotate(origin=np.array([0,0,0]), axis=np.array([0,0,1]), angle=rotate_angle)
        axis3 = rotation_3d(axis3, origin=np.array([0,0,0]), axis=np.array([0,0,1]), angle=rotate_angle)
        sec3.rotate(origin=axis2[0], axis=axis2[1] - axis2[0], angle=-2 * rotate_angle)
        axis3 = rotation_3d(axis3, origin=axis2[0], axis=axis2[1] - axis2[0], angle=-2 * rotate_angle)
        sec3.rotate(origin=axis3[0], axis=axis3[1] - axis3[0], angle=2*rotate_angle)
        sec3.flip(axis='+Y')

        sec3.output_tecplot(fname='plot/sec3.dat', one_piece=True)

        #* ============================================
        #* Output for animation
        #* ============================================

        with open('plot/plot.mcr', 'r') as fin, open('plot/nplot.mcr', 'w') as fout:
            lines = fin.readlines()
            lines[10] = '$!ExportSetup ExportFName = \'%d.png\'\n' % int(irot+75)
            print(lines[10])
            fout.writelines(lines)
        
        print(os.getcwd())
        os.system('tec360 -b -p nplot.mcr')

