import os
import sys
sys.path.append('.')

from cst_modeling.surface import Surface


if __name__ == "__main__":
    
    path = os.path.dirname(sys.argv[0])


    wing = Surface(n_sec=4, name='Propeller', nn=101, ns=51)

    wing.read_setting(os.path.join(path, 'Propeller.txt'), tail=[0.1, 0.1, 0.1, 0.1])

    wing.geo()

    wing.bend(2, 3, leader=[[3.27, 1.3, 1.00, 1.23]], kx=[1.0, 1.0], ky=[1.0, 1.0], rot_x=True)

    #wing.smooth(0,3)
    #wing.smooth(2,3)

    wing.output_tecplot(fname=os.path.join(path, 'Propeller.dat'), one_piece=False, split=False)

