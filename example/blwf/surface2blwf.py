
from matplotlib.pyplot import axis
import numpy as np
from cst_modeling.basic import rotate
from cst_modeling.foil import unify_foil, cst_foil, cst_foil_fit
from cst_modeling.tools.blwf import BLWF, output_curve, read_tecplot, extract_slice
from cst_modeling.surface import Surface
import os
from scipy.spatial.distance import cdist


def extract_wing_section_info(section: np.array, n_cst=10):
    '''
    Extract wing section information (XY plane)
    
    LE XYZ, chord length, twist(deg), rel-thick, CST
    
    section [n,3]
    '''
    xTE = 0.5*(section[0,0]+section[-1,0])
    yTE = 0.5*(section[0,1]+section[-1,1])
    dis = np.sum((section[:,:2] - np.array([xTE,yTE]))**2, axis=1)
    
    iLE = np.argmax(dis)
    xLE = section[iLE,0]
    yLE = section[iLE,1]
    zLE = section[iLE,2]
        
    if section[int(1.2*iLE),1] > section[int(0.8*iLE),1]:
        xu = section[iLE:,0]
        yu = section[iLE:,1]
        xl = section[:iLE+1,0]; xl=xl[::-1]
        yl = section[:iLE+1,1]; yl=yl[::-1]
    else:
        xl = section[iLE:,0]
        yl = section[iLE:,1]
        xu = section[:iLE+1,0]; xu=xu[::-1]
        yu = section[:iLE+1,1]; yu=yu[::-1]

    xu, yu, xl, yl, twist, chord, tail = unify_foil(xu, yu, xl, yl)
    
    cst_u, cst_l = cst_foil_fit(xu, yu, xl, yl, n_cst=n_cst)
    
    _, _, _, thick, rLE = cst_foil(501, cst_u, cst_l)
    
    sec_info = {'xLE':xLE, 'yLE':yLE, 'zLE':zLE, 'chord':chord, 'twist':twist,
                'rel-thick':thick, 'rel-tail':tail, 'rel-rLE':rLE, 'cst_u':cst_u, 'cst_l':cst_l}

    return sec_info

def scale_slice(section: np.ndarray):
    '''
    Rotate and scale sliced section to unit chord length airfoil (Z=const)
    '''

    dis = cdist(section[:,:2], section[0:1,:2], metric='euclidean')
    iLE = np.argmax(dis[:,0])
    xLE = section[iLE,0]
    yLE = section[iLE,1]
    xTE = 0.5*(section[0,0]+section[-1,0])
    yTE = 0.5*(section[0,1]+section[-1,1])
    
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


if __name__ == '__main__':
    
    np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=200)
    
    
    zone_wing = [46, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59,
                 60, 61, 62, 63, 64, 65, 67, 68, 70, 71, 72, 73, 
                 74, 75, 76, 77, 79, 85, 87, 88,
                 91, 92, 95, 96, 97, 98]
    
    zone_wTE  = [81, 82, 83, 84, 89, 90, 93, 94]
    
    zone_tail = [100, 102, 103, 104, 105, 106, 107, 108, 110, 111,
                 113, 115, 116, 117, 118, 119, 120, 121, 123, 124,
                 134, 135, 136, 137, 139, 140]
    
    zone_wing = [i-1 for i in zone_wing]
    zone_wTE  = [i-1 for i in zone_wTE ]
    zone_tail = [i-1 for i in zone_tail]
    zone_body = []
    for i in range(140):
        if (i in zone_wing) or (i in zone_tail) or (i in zone_wTE):
            continue
        zone_body.append(i)
    
    #* Rotate axis
    if True:
        
        index_xyz = [0,1,2]
        
        data, name_var, _ = read_tecplot('surface.dat')
    
        for i in range(len(data)):
            temp = data[i][:,:,:,index_xyz[1]].copy()
            data[i][:,:,:,index_xyz[1]] = data[i][:,:,:,index_xyz[2]] 
            data[i][:,:,:,index_xyz[2]] = - temp
    
        with open('surface-aircraft.dat', 'w') as f:
            f.write('Variables= X Y Z I J K U V W P T M Cp ut \n')
            for z in range(len(data)):
                f.write('zone i= %d j= %d k= %d \n'%(
                    data[z].shape[0], data[z].shape[1], data[z].shape[2]))
                for k in range(data[z].shape[2]):
                    for j in range(data[z].shape[1]):
                        for i in range(data[z].shape[0]):
                            for v in range(data[z].shape[3]):
                                f.write(' %19.9E'%(data[z][i,j,k,v]))
                            f.write('\n')
                f.write('\n')
    
    #* Wing sections
    if True:
        
        Pref = np.array([20, 0, 0])
        dir_norm = np.array([0, 0, -1])
        locations = [4, 8, 10.87, 13, 17, 21, 25, 29.3]
        sections, name_var = extract_slice(locations, Pref, dir_norm, 
                        fname='surface-aircraft.dat', zone_id=zone_wing,
                        dir_ref=np.array([-1.,0.01,0.]), arrange_method='join')

        loc_extra = [2.5, 29.4]
        
        ratio = (loc_extra[0]-locations[0])/(locations[1]-locations[0])
        sec0 = (1-ratio)*sections[0] + ratio*sections[1]
        ratio = (loc_extra[1]-locations[-1])/(locations[-2]-locations[-1])
        sec1 = (1-ratio)*sections[-1] + ratio*sections[-2]
        
        sections = [sec0] + sections + [sec1]

        sec_infos = []
        output_curve(None, fname='wing-slice.dat', append=False, name_var=name_var)
        for section in sections:
            output_curve(section, fname='wing-slice.dat', append=True, name_var=name_var)
            sec_info = extract_wing_section_info(section, n_cst=10)
            sec_infos.append(sec_info)

        print('       xLE       yLE       zLE     chord  twist(deg) rel-thick    tail      rLE')
        for p in sec_infos:
            print('  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.4f'%(
                p['xLE'], p['yLE'], -p['zLE'], p['chord'], p['twist'],
                p['rel-thick'], p['rel-tail']*p['chord'], p['rel-rLE']))
        
        print('CST_coefs:')
        for i in range(len(sec_infos)):
            print('Section %2d---------------'%(i+1))
            print(sec_infos[i]['cst_u'])
            print(sec_infos[i]['cst_l'])
    
    blwf = BLWF(ITV=0)
    
    #* Fuselage sections
    if False:

        body_secs = blwf.define_fuselage(zone_body, n_slice=51, n_point=51,
                            fname='surface-aircraft.dat', index_xyz=[0,1,2])

        blwf.write_input_file()
    
        output_curve(None, fname='curve-body.dat', append=False)
        for curve in body_secs:
            output_curve(np.array(curve), fname='curve-body.dat', append=True)
    

    #* Run BLWF
    if True:
        
        Minf = 0.85
        AoA = 2.3

        wing = Surface(n_sec=9, name='Wing', nn=73, ns=11)
        wing.read_setting('Wing.txt')
        
        for sec in wing.secs:
            sec.zLE = -sec.zLE
            
        wing.geo()
        #wing.output_tecplot('wing.dat')
        
        blwf.cst_wing(wing)
        blwf.update_input_wing()

        lines = []
        with open('blwf.in', 'r') as f:
            lines = f.readlines()
        f = open('blwf.in', 'w')
        f.write(lines[0])
        f.write(lines[1])
        f.write(' %9.4f %9.4f   0.     32000000. %9.4f %9.4f\n'%(
            Minf, AoA, 7.0053, 191.8450))
        for line in lines[3:]:
            f.write(line)
        f.close()

        os.system('blwf58.exe blwf.in')
    
    
    #* Extract slices
    if True:
        
        ratios = np.array([0.201, 0.283, 0.502, 0.603, 0.846, 0.950])
        locations = ratios*29.33
        
        sections, name_var = extract_slice(locations, Pref, dir_norm, 
                        fname='surface-aircraft.dat', zone_id=zone_wing,
                        dir_ref=np.array([-1.,0.01,0.]), arrange_method='join')

        with open('Cp-cfl3d.dat', 'w') as f:
            
            f.write('Variables= X Y Cp \n')
            
            for k in range(len(sections)):
                
                nn = sections[k].shape[0]
                
                sec = np.concatenate((sections[k][:,:3],sections[k][:,12:13]),axis=1)

                xx, yy, cp = scale_slice(sec)
                
                CL = extract_CL(xx, cp)

                f.write('zone T="Loc= %.2f CL= %.2f" i= %d \n'%(locations[k], CL, nn))
                for i in range(nn):
                    f.write('% 20.10f  %20.10f  %20.10f\n'%(xx[i], yy[i], cp[i]))

        
        sections, name_var = extract_slice(locations, Pref, -dir_norm,
                        fname='tp1.plt', 
                        dir_ref=np.array([-1.,0.01,0.]), arrange_method='join')

        with open('Cp-blwf.dat', 'w') as f:
            
            f.write('Variables= X Y Cp \n')
            
            for k in range(len(sections)):
                
                nn = sections[k].shape[0]

                xx, yy, cp = scale_slice(sections[k])
                
                CL = extract_CL(xx, cp)

                f.write('zone T="Loc= %.2f CL= %.2f" i= %d \n'%(locations[k], CL, nn))
                for i in range(nn):
                    f.write('% 20.10f  %20.10f  %20.10f\n'%(xx[i], yy[i], cp[i]))

    
    
    
    
            
