import numpy as np
import matplotlib.pyplot as plt

file_dump = './dump/'


def fitting_blade_sections():
    '''
    Transform 3D blade sections to 2D planar curves and fit with CST.
    
    For turbomachinery, the sections are a 3D curve on cylinders.
    The origin is (x,y,z)=(0,0,0), or (r,t,z)=(0,0,0). (t~theta, rad).
    '''
    from cst_modeling.basic import fromCylinder
    from cst_modeling.section import fit_curve_with_twist, find_circle_3p
    from cst_modeling.surface import OpenSurface

    n_sec = 12

    blade = OpenSurface(n_sec=n_sec, name='Blade', nn=101, ns=101, projection=False)

    CST = []
    origins = []

    print('Layout')
    with open('./files/fan-blade-sections.dat', 'r') as f:
        line = f.readline()    # Variables= X Y Z

        for i in range(n_sec):

            #* Original curve on cylinder
            line = f.readline().split()
            n_point = int(line[2])
            x = np.zeros(n_point)
            y = np.zeros(n_point)
            z = np.zeros(n_point)
            for j in range(n_point):
                line = f.readline().split()
                x[j] = float(line[0])
                y[j] = float(line[1])
                z[j] = float(line[2])
            line = f.readline() # Empty line

            if x[0]>x[-1]:
                x = np.flip(x)
                y = np.flip(y)
                z = np.flip(z)

            #* Locate Origin of cylinder
            ii = int(0.5*n_point)
            _, origin = find_circle_3p([x[0], y[0]], [x[ii], y[ii]], [x[-1], y[-1]])
            origins.append(origin)
            
            #* Convert to plane curve
            xx, yy, zz = fromCylinder(x, y, z, flip=True, origin=origin)

            XLE = xx[0]
            YLE = yy[0]
            ZLE = zz[0]

            #* CST coefficients
            cst, chord, twist, thick = fit_curve_with_twist(xx, yy, n_cst=7)
            CST.append(cst)

            print(np.array([XLE, YLE, ZLE, chord, twist, thick]))

            blade.secs[i].xLE   = XLE
            blade.secs[i].yLE   = YLE
            blade.secs[i].zLE   = ZLE
            blade.secs[i].chord = chord
            blade.secs[i].twist = twist
            blade.secs[i].thick = thick
            blade.secs[i].cst   = cst.copy()

    print()
    print('CST Parameters')
    for cst in CST:
        print(cst)
        
    print()
    print('Cylinder Origins')
    for origin in origins:
        print(origin)

    blade.update_sections()
    blade.surf_to_cylinder(flip=True, origin=origins)
    blade.output_section(file_dump+'fitting_blade_sections.dat', TwoD=False)

def fan_blade():
    
    from cst_modeling.surface import Surface
    
    blade = Surface(n_sec=6, name='Blade-simple',nn=101, ns=51, projection=False)

    blade.read_setting('./files/Fan.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.1, 0.05])

    blade.geo()

    blade.smooth(isec0=0, isec1=4)
    blade.smooth(isec0=4, isec1=5, smooth0=True)

    blade.surf_to_cylinder(flip=True)

    blade.output_tecplot(fname=file_dump+'Blade-simple.dat')

def fan_blade_split():
    
    from cst_modeling.surface import OpenSurface
    
    #* ==========================================
    #* Suction side
    blade1 = OpenSurface(n_sec=7, name='Blade-suction',nn=101, ns=51, projection=False)

    blade1.read_setting('./files/Fan.txt')

    origins = None
    origins = blade1.read_cylinder_origins('Fan.txt')

    blade1.geo()

    # Convert back to cylinder
    blade1.surf_to_cylinder(flip=True, origin=origins)

    # This outputs the surface to fname in tecplot format
    blade1.output_tecplot(fname=file_dump+'Blade-suction.dat')

    #* ==========================================
    #* Pressure side
    blade2 = OpenSurface(n_sec=7, name='Blade-pressure',nn=101, ns=51, projection=False)

    blade2.read_setting('./files/Fan.txt')

    origins = blade2.read_cylinder_origins('Fan.txt')

    blade2.geo()

    blade2.surf_to_cylinder(flip=True, origin=origins)

    blade2.output_tecplot(fname=file_dump+'Blade-pressure.dat')


if __name__ == '__main__':
    
    np.set_printoptions(formatter={'float': '{: 10.6f}'.format}, linewidth=200)
    
    fitting_blade_sections()
    
    fan_blade()
    
    fan_blade_split()
