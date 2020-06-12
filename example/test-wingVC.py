
# Run this in the directory where the folder cst_modeling is
from cst_modeling.auxiliary import WingVariableCamber


if __name__ == "__main__":

    print('This is a example for constructing a wing')

    flap_loc = [5.0, 10.0, 18.0, 23.0]
    flap_angle = [-10.0, 10.0]
    axis_xloc = [0.8, 0.7]
    axis_dy = []

    wingVC = WingVariableCamber(n_sec=4, name='WingVC', fname='Wing.txt',
        nn=201, ns=51, flap_loc=flap_loc, flap_trans=0.2, flap_angle=flap_angle,
        axis_xloc=axis_xloc, axis_dy=axis_dy)

    wingVC.build(split=True, one_piece=False, f_tecplot=None, f_plot3d=None)

    # This is for constructing surfaces with curved leading edge lines
    # e.g., strut of the SBW aircraft, wing-let
    wingVC.bend(isec0=wingVC.n_sec-2, isec1=wingVC.n_sec-1,
        leader=[[22.9, 1.2, 27.1, 0.75]], kx=[0.4, 1.0], ky=[0.0, 1.6], rot_x=True)

    wingVC.output_tecplot(fname='Wing.dat', one_piece=False)

