

from cst_modeling.surface import Surface, plot3d_to_igs


if __name__ == "__main__":

    #* Define a [Surface] object, which is called wing
    #  It has n_sec control sections, and its name is 'Wing'
    #  Each section's upper/lower surface are constructed by CST method
    #  The chord-wise has nn points in the upper/lower surface
    #  The span-wise has ns points in each surface between two adjcent sections
    wing = Surface(n_sec=3, name='Wing-basic', nn=101, ns=52)


    #* Read settings from file 'Wing.txt'
    #  The settings of [wing] object is under its name 'Wing'
    #  First part is the parameters for the layout
    #  Second part is the CST parameters of upper/lower surface of each section
    wing.read_setting('Wing.txt')


    #* Add auxiliary section
    wing.add_sec([5.0], axis='Z')


    #* Constructs the surfaces between sections
    wing.geo()


    #* Output tecplot format to fname
    # one_piece is an option for combining all surfaces in different sections into one piece
    # split is an option for splitting to upper and lower surfaces
    wing.output_tecplot(fname='Wing-basic.dat', one_piece=False, split=True)

    wing.output_plot3d(fname='Wing-basic.grd', split=True)

    wing.output_section(fname='Wing-basic-sec.dat', TwoD=False)

    plot3d_to_igs(fname='Wing-basic')


    #* Plot the wing surface on screen
    #  This is flipping the surface, it can do turing 90deg or mirror
    wing.flip(axis='+X +Z')

    wing.plot()

