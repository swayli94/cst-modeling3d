# cst-modeling

Surface and curve modeling via CST (class shape transformation) method

The curves, e.g., foil's upper and lower surface, are constructed via CST method.

The multi-section surface is interpolated from several control sections.

## Modules

### foil

1. Section class:

    This is a class for control sections.

    The section stores unit 2D curves of upper and lower surface (xx, yu, yl) and CST coefficients (cst_u, cst_l). The 3D curve (x, y, z) is generated from the 2D curves. It starts from the lower surface trailing edge and ends at the upper surface trailing edge. The 3D curve is still a plane curve, the leading edge location (xLE, yLE, zLE), chord length, twist angle (deg), maximum relative thickness, relative tail thickness are defined by user.

2. Supportive functions:

    There are several functions dealing with curves.

        [curve_curvature]       get curvature of a curve
        [transform]             2D curve's translation, scale, rotation in z-axis
        [rotate]                3D curve's rotation defined by given origin, angle and axis
        [stretch_fixed_point]   linearly stretch a curve with certain fixed point 
        [output_foil]           output the airfoil/curve in tecplot format

### surface

1. Surface class:

    This is a class for multi-section surface.

        [__init__]              number of control sections, cst parameters of each curve must be specified immediately
        [read_setting]          layout and cst parameters can be read in from a control file
        [geo]                   the function of generating curves and surfaces
        [add_sec]               add sections to the surface, the new sections are interpolated from current ones
        [flip]                  mirror the surface by certain plane, or rotate 90 deg by certain axis
        [bend]                  bend the surface in x-axis and y-axis, also can rotate sections in x-axis
        [smooth]                smooth the surface
        [toCylinder]            convert the plane curves to curves on cylinders (for turbomachinery)
        [fromCylinder]          convert curves on cylinder to plane curves
        [plot]                  show the surface
        [output_tecplot]        output surface in tecplot format
        [output_plot3d]         output surface in plot3d format
