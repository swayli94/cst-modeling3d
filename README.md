# cst_modeling

Surface and curve modeling via CST (class shape transformation) method

The curves, e.g., foil's upper and lower surface, are constructed via CST method.

The multi-section surface is interpolated from several control sections.

## Installation

``` bash
pip install cst-modeling3d
```

## Functions

**Note:** Details for usage can be found in code comments.

### Airfoil

1. build an airfoil

    CST coefficients of upper and lower surface can have different sizes.

    ``` python
    >>> from cst_modeling.foil import cst_foil
    >>> nn = 1001
    >>> cst_u = [ 1.0,  1.0,  1.0]
    >>> cst_l = [-1.0, -1.0, -1.0]
    >>> x, yu, yl, t0, R0 = cst_foil(nn, cst_u, cst_l)
    ```

2. find CST coefficients of an airfoil

    cst_u, cst_l are of the same size

    ``` python
    >>> from cst_modeling.foil import cst_foil_fit
    >>> cst_u, cst_l = cst_foil_fit(x, yu, x, yl, n_order=7)
    ```

3. get geometry information of an airfoil

    thickness/camber distribution, curvature of upper/lower surface

    ```python
    >>> from cst_modeling.foil import foil_tcc
    >>> thickness, curv_u, curv_l, camber = foil_tcc(x, yu, yl)
    ```

4. check whether the airfoil is reasonable or not

    ```python
    >>> from cst_modeling.foil import check_valid
    >>> rule_invalid = check_valid(x, yu, yl, RLE)
    ```

5. modify airfoil by adding a bump

    ```python
    >>> from cst_modeling.foil import foil_bump_modify
    >>> yu_new, yl_new = foil_bump_modify(x, yu, yl, xc, h, s, side)
    ```

6. modify curve by adding an incremental curve

    ```python
    >>> from cst_modeling.foil import foil_increment
    >>> yu_new, yl_new = foil_increment(x, yu, yl, cst_u, cst_l)
    ```

7. output airfoil

    ```python
    >>> from cst_modeling.foil import output_foil
    >>> output_foil(x, yu, yl, fname='airfoil.dat', ID=0, info=False)
    ```

8. 2D curve's translation, scale, and rotation in z-axis

    ```python
    >>> from cst_modeling.foil import transform
    >>> xu_, xl_, yu_, yl_ = transform(x, x, yu, yl, scale=1.0)
    ```

9. rotate the 3D curve according to origin

    ```python
    >>> from cst_modeling.foil import rotate
    >>> x_, y_, z_ = rotate(x, y, z, angle=0.0, axis='X')
    ```

10. linearly stretch a curve when certain point is fixed

    ```python
    >>> from cst_modeling.foil import stretch_fixed_point
    >>> x_, y_ = stretch_fixed_point(x, y, dx=0.0, dy=0.0)
    ```

### Wing with variable camber or flap

**WingVariableCamber** is a sub-class of Surface.

```text
+x:     flow direction (m)
+y:     upside (m)
+z:     span-wise (m)
twist:  +z direction (deg)
chord:  chord length (m)
thick:  relative maximum thickness
tail:   absolute tail thickness (m)
```

1. Initialization

    ```python
    >>> from cst_modeling.auxiliary import WingVariableCamber
    >>> flap_loc = [5.0, 10.0, 18.0, 23.0]
    >>> flap_angle = [-10.0, 10.0]
    >>> axis_xloc = [0.8, 0.7]
    >>> wingVC = WingVariableCamber(n_sec=4, name='WingVC', fname='Wing.txt', nn=201, ns=51,
            flap_loc=flap_loc, flap_trans=0.2, flap_angle=flap_angle, axis_xloc=axis_xloc)
    ```

    Parameters for setting up wing and flap.

    ```text
    n_sec:   number of control sections (2D if set to 0 or 1)
    tail:    tail thickness (m), list or float
    name:    name of the surface
    fname:   name of control file
    nn:      number of points of upper/lower section
    ns:      number of span-wise points
    project: True ~ projected chord length does not change when twisted

    flap_loc:   list [2*n_flap], z coordinates of the flap ends. [z_flap1_1, z_flap1_2, z_flap2_1, z_flap2_2, ...]
    flap_trans: float, transition length of the flap deflection
    flap_angle: list [n_flap], deflect angle of flaps
    axis_xloc:  list [n_flap], chord-wise ratio of the flap rotation axis
    axis_dy:    empty list or list [n_flap], scaled y location of the rotation axis
    ```

2. Building the geometry

    ```python
    >>> wingVC.build(split=True, one_piece=False, f_tecplot='Wing.dat', f_plot3d='Wing.grd')
    ```

    ```text
    split:      True ~ generate [surfs] as upper and lower separately
    showfoil:   True ~ output name-foil.dat of airfoils
    one_piece:  True ~ combine the span-wise sections into one piece (for tecplot format)

    f_tecplot:  file name of tecplot format file. If None, do not output.
    f_plot3d:   file name of tecplot format file. If None, do not output.
    ```

## Modules and classes

**Note:** Details for usage can be found in code comments.

### Section

**Section** is a class for control sections.

The section stores unit 2D curves of upper and lower surface (xx, yu, yl), as well as CST coefficients (cst_u, cst_l).

The 3D curve (x, y, z) is then generated from the 2D curves. It starts from the lower surface trailing edge and ends at the upper surface trailing edge.

The 3D curve is still a plane curve, the leading edge location (xLE, yLE, zLE), chord length, twist angle (deg), maximum relative thickness, relative tail thickness are defined by user.

```python
>>> from cst_modeling.foil import Section
>>> sec = Section(thick=None, chord=1.0, twist=0.0, tail=0.0)
>>> sec.foil(cst_u=cst_u, cst_l=cst_l, nn=1001)
```

### Multi-section surface

**Surface** is a class for multi-section surface.

The surface class can be used for both wings and fan blades.

```python
>>> from cst_modeling.surface import Surface
>>> wing = Surface(n_sec=3, name='Wing', nn=101, ns=51, project=True)
>>> wing.read_setting('Wing.txt', tail=[0.1, 0.1, 0.1])
```

Parameters:

```text
n_sec:   number of control sections (2D if set to 0 or 1)
tail:    tail thickness (m), list or float
name:    name of the surface
fname:   name of control file
nn:      number of points of upper/lower section
ns:      number of span-wise points
project: True ~ projected chord length does not change when twisted
```

1. construct geometry

    ```python
    >>> wing.geo(split=True, showfoil=True)
    ```

    ```text
    showfoil:   if True, output name-foil.dat of airfoils
    split:      if True, generate surfaces as upper and lower separately
    ```

2. output the surface in tecplot format

    ```python
    >>> wing.output_tecplot(fname='Wing.dat', one_piece=False)
    ```

    ```text
    one_piece:  if True, combine the span-wise sections into one piece
    ```

3. output the surface in plot3d format

    ```python
    >>> wing.output_plot3d(fname='Wing.grd')
    ```

4. plot surface by matplotlib

    ```python
    >>> wing.plot(fig_id=1, type='wireframe')
    ```

5. flipping the surface (turing 90 degrees or mirror)

    **Note:** This should be the last action before output

    ```python
    >>> wing.flip(axis='+X +Z')
    ```

6. add sections to the surface

    the new sections are interpolated from current ones

    **Note:** must run before geo() and flip()

    ```python
    >>> wing.add_sec(location=[8, 12.0], axis='Z')
    ```

7. bend surfaces by angle and leader curve

    This is for constructing surfaces with curved leading edge lines,
    e.g., strut of the SBW aircraft, wing-let.

    ```python
    >>> wing.bend(******)
    ```

8. smooth the surface between given sections

    ```python
    >>> wing.smooth(isec0, isec1)
    ```

9. convert actual coordinates x,y,z of curves on cylinder to plane curves

    The section of fan blades are usually defined by sections sliced by cylinders.

    Since the multi-section surface is constructed by plane sections,
    the actual control sections on cylinders need to be convert to plane sections.
    After the surfaces are built based on these plane sections,
    the surfaces are convert back to the actual coordinates.

    Use the following function to convert leading edge, trailing edge coordinates,
    then define the multi-section surface layout parameters of the plane sections.

    ```python
    >>> x_, y_, z_ = Surface.fromCylinder(x, y, z, flip=True)
    ```

10. convert plane sections and corresponding surfaces to cylinder coordinates

    Usually project = False for fans.

    ```python
    >>> fan = Surface(n_sec=7, name='Fan',nn=101, ns=51, project=False)
    >>> fan.read_setting('Fan.txt', tail=0.1)
    >>> fan.geo(split=True, showfoil=False)
    >>> fan.Surf2Cylinder(flip=True)
    ```
