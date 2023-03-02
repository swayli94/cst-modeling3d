
import numpy as np
import matplotlib.pyplot as plt


def transform_airfoil():
    '''
    Transform an airfoil or a curve
    '''
    from cst_modeling.basic import transform
    from cst_modeling.section import cst_foil
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])
    
    plt.figure()

    x0, yu0, yl0, _, _ = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.01)
    
    plt.plot(x0, yu0, 'k')
    plt.plot(x0, yl0, 'k')
    
    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, rot=20, projection=True)
    plt.plot(xu_, yu_, 'k--')
    plt.plot(xl_, yl_, 'k--')
    plt.text(-0.5, 0.65, 'rotation when keeping the projection length unchanged', fontsize=8, color='k')
    
    
    xLE = -0.5
    yLE = 0.2
    
    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=20, dx=xLE-0.0, dy=yLE-0.0, 
                                   x0=None, y0=None, xr=None, yr=None, projection=False)
    plt.plot(xu_, yu_, 'b')
    plt.plot(xl_, yl_, 'b')
    plt.text(-0.5, 0.60, 'scale center and rotate center are both the leading edge (%.2f, %.2f)'%(xLE, yLE), fontsize=8, color='b')


    xTE = 1.0
    yTE = 0.2

    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=-10, dx=xTE-1.0, dy=yTE-0.0, 
                                   x0=xTE, y0=yTE, xr=xTE, yr=yTE, projection=False)
    plt.plot(xu_, yu_, 'g--')
    plt.plot(xl_, yl_, 'g--')
    plt.text(-0.5, 0.55, 'scale center and rotate center are both the trailing edge (%.2f, %.2f)'%(xTE, yTE), fontsize=8, color='g')



    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=-10, dx=0.0, dy=0.0, 
                                   x0=0.5, y0=0.0, xr=None, yr=None, projection=False)
    plt.plot(xu_, yu_, 'r--')
    plt.plot(xl_, yl_, 'r--')
    plt.text(-0.5, 0.50, 'scale center is (0.5, 0.0), rotate center is the leading edge', fontsize=8, color='r')


    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('figures/transform_airfoil.jpg', dpi=300)
    plt.close()

def normalize_airfoil():
    '''
    Normalize an airfoil
    '''
    from cst_modeling.basic import transform
    from cst_modeling.section import cst_foil, normalize_foil
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])
    
    x0, yu0, yl0, t0, r0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.01)

    xu1, xl1, yu1, yl1 = transform(x0, x0, yu0, yl0, scale=1.5, rot=20, dx=-0.4, dy=0.1)
    
    xu2, yu2, xl2, yl2, twist, chord, tail = normalize_foil(xu1, yu1, xl1, yl1)

    plt.figure()
    plt.plot(x0, yu0, 'k')
    plt.plot(xu1, yu1, 'b')
    plt.plot(xu2, yu2, 'g--')
    
    plt.legend(['original', 'transformed', 'normalized'])
    
    plt.plot(x0, yl0, 'k')
    plt.plot(xl1, yl1, 'b')
    plt.plot(xl2, yl2, 'g--')
    
    plt.text(0.0, 0.10, 'chord= %.2f, twist= %.2f, tail= %.3f'%(chord, twist, tail), fontsize=12, color='k')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('figures/normalize_airfoil.jpg', dpi=300)
    plt.close()

def stretch_curve():
    '''
    Linearly stretch a curve when a certain point is fixed
    '''
    
    from cst_modeling.basic import stretch_fixed_point
    from cst_modeling.section import cst_foil
    
    cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
    cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])
    
    plt.figure()

    x0, yu0, yl0, _, _ = cst_foil(1001, cst_u, cst_l, x=None, t=0.11, tail=0.01)
    
    plt.plot(x0, yu0, 'k')
    plt.plot(x0, yl0, 'k')
    
    
    dx=0.5; dy=0.5; xm=1.0; ym=0.0; xf=0.0; yf=0.0
    xu_, yu_ = stretch_fixed_point(x0, yu0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    xl_, yl_ = stretch_fixed_point(x0, yl0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    
    plt.plot(xu_, yu_, 'b')
    plt.plot(xl_, yl_, 'b')
    plt.plot([xf], [yf], 'bo')
    plt.plot([xm, xm+dx], [ym, ym+dy], 'b--')
    plt.text(0.0, 0.45, 'stretch (%.1f, %.1f) to (%.1f, %.1f), fixed point (%.2f, %.2f)'%(xm, ym, xm+dx, ym+dy, xf, yf), fontsize=12, color='b')
    
    
    dx=0.2; dy=0.0; xm=0.0; ym=0.0; xf=x0[500]; yf=yu0[500]
    xu_, yu_ = stretch_fixed_point(x0, yu0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    xl_, yl_ = stretch_fixed_point(x0, yl0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    
    plt.plot(xu_, yu_, 'g')
    plt.plot(xl_, yl_, 'g')
    plt.plot([xf], [yf], 'go')
    plt.plot([xm, xm+dx], [ym, ym+dy], 'g--')
    plt.text(0.0, 0.40, 'stretch (%.1f, %.1f) to (%.1f, %.1f), fixed point (%.2f, %.2f)'%(xm, ym, xm+dx, ym+dy, xf, yf), fontsize=12, color='g')
    
    
    dx=0.0; dy=0.5; xm=0.0; ym=0.0; xf=0.5; yf=0.0
    xu_, yu_ = stretch_fixed_point(x0, yu0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    xl_, yl_ = stretch_fixed_point(x0, yl0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf)
    
    plt.plot(xu_, yu_, 'r')
    plt.plot(xl_, yl_, 'r')
    plt.plot([xf], [yf], 'ro')
    plt.plot([xm, xm+dx], [ym, ym+dy], 'r--')
    plt.text(0.0, 0.35, 'stretch (%.1f, %.1f) to (%.1f, %.1f), fixed point (%.2f, %.2f)'%(xm, ym, xm+dx, ym+dy, xf, yf), fontsize=12, color='r')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('figures/stretch_curve.jpg', dpi=300)
    plt.close()

def cylinder_to_plane():
    '''
    Transform a curve on a cylinder to a curve on a plane
    '''
    from cst_modeling.basic import toCylinder, fromCylinder, rotate
    from cst_modeling.section import cst_foil

    fig = plt.figure(figsize=(10,4))

    def _figure(ax, zLE, thetaLE, radius, origin):
    
        u = np.linspace(0, 2*np.pi, 101, endpoint=True)
        h = np.linspace(0, radius, 20)
        x_cylinder = np.outer(np.sin(u), np.ones_like(h)*radius)
        y_cylinder = np.outer(np.cos(u), np.ones_like(h)*radius)
        z_cylinder = np.outer(np.ones_like(u), h)
        
        
        cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
        cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

        X, YU, YL, _, _ = cst_foil(1001, cst_u, cst_l)
        
        xu, yu, zu = toCylinder(X, YU, np.ones_like(X)*radius, flip=False, origin=origin)
        xl, yl, zl = toCylinder(X, YL, np.ones_like(X)*radius, flip=False, origin=origin)
        
        xu, yu, zu = rotate(zu, yu, zu, thetaLE, axis='Z')
        xl, yl, zl = rotate(xl, yl, zl, thetaLE, axis='Z')
        zu = zu + zLE
        zl = zl + zLE
        
        XU_, YU_, ZU_ = fromCylinder(xu, yu, zu, flip=False, origin=origin)
        XL_, YL_, ZL_ = fromCylinder(xl, yl, zl, flip=False, origin=origin)

        ax.plot3D(xu, yu, zu, 'k')
        ax.plot3D(XU_, YU_, ZU_, 'r')

        ax.plot3D([-0.2, 0.2], [0, 0], [radius, radius], 'r--')
        ax.plot3D([0, 0], [-0.2, 0.2], [radius, radius], 'r--')

        ax.plot3D(xl, yl, zl, 'k')
        ax.plot3D(XL_, YL_, ZL_, 'r')

        ax.plot_surface(x_cylinder, y_cylinder, z_cylinder, 
                        alpha=0.2, rstride=1, cstride=1, color='gray')
        
        plt.legend(['cylinder curve', 'plane curve', 'origin of plane'])
        
        plt.title('r=%.1f, zLE=%.1f, thetaLE=%.1f'%(radius, zLE, thetaLE))

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim3d(origin[0]-radius*1.2, origin[0]+radius*1.2)
        ax.set_ylim3d(origin[1]-radius*1.2, origin[1]+radius*1.2)
        ax.set_zlim3d(0.0, 1.2*radius)
        

    ax = fig.add_subplot(1,2,1,projection='3d')
    
    _figure(ax, zLE=0.5, thetaLE=0.0, radius=1.0, origin=[0, 0])
    
    ax.view_init(elev=32, azim=-130)
    
    
    ax = fig.add_subplot(1,2,2,projection='3d')
    
    _figure(ax, zLE=0.5, thetaLE=0.0, radius=2.0, origin=[0, 0])
    
    ax.view_init(elev=32, azim=-130)
    
    
    plt.savefig('figures/cylinder_to_plane.jpg', dpi=300)
    plt.close()
    
    
if __name__ == '__main__':
    
    transform_airfoil()
    
    normalize_airfoil()

    stretch_curve()
    
    cylinder_to_plane()
    