import numpy as np
import matplotlib.pyplot as plt

def intersection_curve():
    '''
    Intersection between curves
    '''

    from cst_modeling.basic import intersect_index, intersect_point
    
    x1 = np.linspace(0.5, 1.0, 151)
    x2 = np.linspace(0.5, 1.0, 101)
    y1 = np.sin(x1*np.pi)
    y2 = 1.2 - x2
    
    i1, i2, points = intersect_index(x1, y1, x2, y2)
    
    p1 = np.array([x1[i1],   y1[i1]])
    p2 = np.array([x1[i1+1], y1[i1+1]])
    p3 = np.array([x2[i2],   y2[i2]])
    p4 = np.array([x2[i2+1], y2[i2+1]])
    intersection = intersect_point(p1, p2, p3, p4)
    
    
    plt.figure()
    plt.plot(x1, y1, 'k')
    plt.plot(x2, y2, 'b')
    plt.plot(points[0][0], points[0][1], 'ko')
    plt.plot(points[1][0], points[1][1], 'bo')
    plt.plot(intersection[0], intersection[1], 'r*')
    plt.text(0.5, 0.15, 'Intersection index %d (black) and %d (blue)'%(i1, i2), fontsize=12, color='k')
    plt.text(0.5, 0.05, 'Intersection point (%.3f, %.3f)'%(intersection[0], intersection[1]), fontsize=12, color='k')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    
    plt.savefig('figures/intersection_curve.jpg', dpi=300)
    plt.close()

def intersection_vector_plane():
    '''
    Intersection between vector and plane
    '''
    
    from cst_modeling.basic import intersect_vec_plane, plane_3points

    V0 = np.array([0.0, 0.0, 0.0])
    V1 = np.array([1.0, 1.0, 1.0])
    P0 = np.array([0.0, 0.0, 1.0])
    P1 = np.array([0.0, 1.0, 0.0])
    P3 = np.array([1.0, 0.0, 0.0])
    
    uu = np.linspace(0, 1, 101)
    u, v = np.meshgrid(uu, uu)
    w = plane_3points(P0, P1, P3, u, v)
        
    xi, t1, t3, rv = intersect_vec_plane(V0, V1, P0, P1, P3)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    ax.plot_surface(u, v, w, alpha=0.2, rstride=1, cstride=1, color='gray')
    
    ax.plot3D([P0[0], P1[0]], [P0[1], P1[1]], [P0[2], P1[2]], 'k')
    ax.plot3D([P0[0], P3[0]], [P0[1], P3[1]], [P0[2], P3[2]], 'k')
    ax.plot3D([V0[0], V1[0]], [V0[1], V1[1]], [V0[2], V1[2]], 'b--')
    
    ax.plot3D([xi[0]], [xi[1]], [xi[2]], 'ro')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim3d(0.0, 1.0)
    ax.set_ylim3d(0.0, 1.0)
    ax.set_zlim3d(0.0, 1.0)
    ax.view_init(elev=20, azim=120)
    
    plt.title('Ratio in V01: %.2f; P01: %.2f; P03: %.2f'%(rv, t1, t3))
    plt.savefig('figures/intersection_vector_plane.jpg', dpi=300)
    plt.close()

def intersection_surface_plane():
    '''
    Intersection between vector and plane
    '''
    
    from cst_modeling.basic import plane_3points, intersect_surface_plane

    def _cylinder(radius, height, origin):
    
        u = np.linspace(0, 2*np.pi, 101, endpoint=True)
        h = np.linspace(0, height, 20)
        cylinder = np.zeros([u.shape[0], h.shape[0], 3])
        cylinder[:,:,0] = np.outer(np.sin(u), np.ones_like(h)*radius)+origin[0]
        cylinder[:,:,1] = np.outer(np.cos(u), np.ones_like(h)*radius)+origin[1]
        cylinder[:,:,2] = np.outer(np.ones_like(u), h)
        
        return cylinder

    #* Plane
    P0 = np.array([0.0, 0.0, 1.0])
    P1 = np.array([0.0, 1.0, 0.0])
    P3 = np.array([1.0, 0.0, 0.0])
    
    uu = np.linspace(0, 1, 101)
    u, v = np.meshgrid(uu, uu)
    w = plane_3points(P0, P1, P3, u, v)
    
    #* Figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    
    ax.plot_surface(u, v, w, alpha=0.2, rstride=1, cstride=1, color='gray')
    ax.plot3D([P0[0], P1[0]], [P0[1], P1[1]], [P0[2], P1[2]], 'k')
    ax.plot3D([P0[0], P3[0]], [P0[1], P3[1]], [P0[2], P3[2]], 'k')
    
    #* Cylinder 1
    origin=[0.5,0.5,0]
    cylinder = _cylinder(radius=0.3, height=1.0, origin=origin)
    
    curve, ij_curve, xi_curve, yt_curve = intersect_surface_plane(
        cylinder, P0, P1, P3, within_bounds=True, original_order=False)
    curve = np.array(curve)
    
    ax.plot_surface(cylinder[:,:,0], cylinder[:,:,1], cylinder[:,:,2], 
                    alpha=0.2, rstride=1, cstride=1, color='blue')
    
    ax.plot3D(origin[0], origin[1], origin[2], 'bx')
    
    ax.plot3D(curve[:,0], curve[:,1], curve[:,2], 'b')
    
    #* Cylinder 2
    origin=[0.1,0.1,0]
    cylinder = _cylinder(radius=0.2, height=1.0, origin=origin)
    
    ax.plot_surface(cylinder[:,:,0], cylinder[:,:,1], cylinder[:,:,2], 
                    alpha=0.2, rstride=1, cstride=1, color='g')
    ax.plot3D(origin[0], origin[1], origin[2], 'gx')
    
    curve, ij_curve, xi_curve, yt_curve = intersect_surface_plane(
        cylinder, P0, P1, P3, within_bounds=True, original_order=False)
    curve = np.array(curve)

    ax.plot3D(curve[:,0], curve[:,1], curve[:,2], 'g')
    
    curve, ij_curve, xi_curve, yt_curve = intersect_surface_plane(
        cylinder, P0, P1, P3, within_bounds=False, original_order=False)
    curve = np.array(curve)

    ax.plot3D(curve[:,0], curve[:,1], curve[:,2], 'g--')
    
    
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim3d(0.0, 1.0)
    ax.set_ylim3d(0.0, 1.0)
    ax.set_zlim3d(0.0, 1.0)
    ax.view_init(elev=20, azim=120)
    
    plt.savefig('figures/intersection_surface_plane.jpg', dpi=300)
    plt.close()


if __name__ == '__main__':
    
    intersection_curve()

    intersection_vector_plane()

    intersection_surface_plane()

