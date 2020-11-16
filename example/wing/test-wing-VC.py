
# Run this in the directory where the folder cst_modeling is
from cst_modeling.auxiliary import WingVariableCamber


if __name__ == "__main__":

    print('This is a example for constructing a wing')

    flap_loc = [5.0, 10.0, 12.0, 25.0]
    flap_angle = [-5.0, 5.0]
    axis_xloc = [0.8, 0.7]
    axis_dy = []

    wingVC = WingVariableCamber(n_sec=6, name='Wing-tip', fname='Wing.txt',
        nn=101, ns=101, flap_loc=flap_loc, flap_trans=0.2, flap_angle=flap_angle,
        axis_xloc=axis_xloc, axis_dy=axis_dy)

    wingVC.build(split=True, one_piece=False, f_tecplot=None, f_plot3d=None)

    wingVC.bend(4, 5, leader=[[21.0, 2.1, 30.0, 1.6]], kx=[0.6983, 4.0], ky=[0.1043, 1.10], rot_x=True)

    wingVC.output_tecplot(fname='Wing-VC.dat', one_piece=False)


