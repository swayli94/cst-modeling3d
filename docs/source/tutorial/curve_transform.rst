Curve transformation
=======================

Transform airfoil
------------------------

The following figure shows the transformation of an airfoil. It takes three steps in order:
translation (`dx`, `dy`), scale (`scale`, `x0`, `y0`), and rotation (`rot`, `xr`, `yr`, `projection`).

.. code-block:: python
    :linenos:
    
    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, rot=20, projection=True)           # black

    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=20, 
        dx=xLE-0.0, dy=yLE-0.0, x0=None, y0=None, xr=None, yr=None, projection=False)   # blue

    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=-10, 
        dx=xTE-1.0, dy=yTE-0.0, x0=xTE, y0=yTE, xr=xTE, yr=yTE, projection=False)       # green

    xu_, xl_, yu_, yl_ = transform(x0, x0, yu0, yl0, scale=1.5, rot=-10,    
        dx=0.0, dy=0.0, x0=0.5, y0=0.0, xr=None, yr=None, projection=False)             # red


.. figure:: ../../tutorial/figures/transform_airfoil.jpg
    :width: 90 %
    :align: center

    Transform airfoil


Normalize airfoil
------------------------

The following figure shows the normalization of any airfoil.

.. code-block:: python
    :linenos:
    
    xu_, yu_, xl_, yl_, twist, chord, tail = normalize_foil(xu, yu, xl, yl)

.. figure:: ../../tutorial/figures/normalize_airfoil.jpg
    :width: 70 %
    :align: center

    Normalize airfoil


Stretch curve
------------------------

The following figure shows the result of stretching a curve.

.. code-block:: python
    :linenos:
    
    xu_, yu_ = stretch_fixed_point(x0, yu0, dx=dx, dy=dy, xm=xm, ym=ym, xf=xf, yf=yf) 

.. figure:: ../../tutorial/figures/stretch_curve.jpg
    :width: 70 %
    :align: center

    Stretch curve


