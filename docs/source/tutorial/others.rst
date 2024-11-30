Others
========================

Input and output
-------------------

There are several formats available for I/O.

.. code-block:: python
    :linenos:
    
    output_curve(x, y, fname='curve.dat', ID=0)

    output_foil(x, yu, yl, fname='airfoil.dat', ID=0, info=False)

    output_plot3d(X, Y, Z, fname, scale=1.0)

    xs, ys = read_curves(fname='curve.dat')

    data, name_var, titles = read_tecplot(fname='tecplot.dat')

    xyz, iLine0_new = read_block_plot3d(lines, iLine0, ni, nj, nk)


Format transformation
-----------------------

Plot3D format to IGES format.

.. code-block:: python
    :linenos:
    
    output_plot3d(X, Y, Z, 'name.xyz', scale=1.0)

    plot3d_to_igs(fname='name')
