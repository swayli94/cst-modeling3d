

from cst_modeling.surface import Surface


if __name__ == "__main__":

    #* Define a [Surface] object, which is called wing
    #  It has n_sec control sections, and its name is 'Wing'
    #  Each section's upper/lower surface are constructed by CST method
    #  The chord-wise has nn points in the upper/lower surface
    #  The span-wise has ns points in each surface between two adjcent sections
    wing = Surface(n_sec=3, name='Wing-basic', nn=101, ns=21)


    #* Read settings from file 'Wing.txt'
    #  The settings of [wing] object is under its name 'Wing'
    #  First part is the parameters fro the layout
    #  Second part is the CST parameters of upper/lower surface of each section
    wing.read_setting('Wing.txt')


    #* Constructs the surfaces between sections
    #  split and showfoil are options
    wing.geo(split=False, showfoil=True)


    #* Output tecplot format to fname
    # one_piece is an option for combining all surfaces in different sections into one piece
    wing.output_tecplot(fname='Wing-basic.dat', one_piece=False)


    #* Plot the wing surface on screen
    #  This is flipping the surface, it can do turing 90deg or mirror
    wing.flip(axis='+X +Z')

    wing.plot()


