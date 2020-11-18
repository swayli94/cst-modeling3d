

from cst_modeling.surface import Surface


if __name__ == "__main__":


    wing = Surface(n_sec=6, name='Wing-tip', nn=101, ns=51)

    wing.read_setting('Wing.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.05, 0.01])

    wing.geo()

    wing.bend(4, 5, leader=[[21.0, 2.1, 30.0, 1.6]], kx=[0.6983, 4.0], ky=[0.1043, 1.10], rot_x=True)

    wing.smooth(0,2)
    wing.smooth(2,4)

    wing.output_tecplot(fname='Wing-tip.dat', one_piece=False, split=False)

