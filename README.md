# CST Modeling

Surface and curve modeling via class shape function transformation (CST) method.

The CST method combines class functions and shape functions to describe an arbitrary geometry and can guarantee airfoil smoothness with comparatively fewer design variables. Usually, a sixth-order Bernstein polynomial is used as the shape function; i.e., seven CST parameters are used to describe upper and lower surfaces.

Reference: Kulfan, B. M., “Universal parametric geometry representation method,” Journal of Aircraft, vol. 45, No. 1, 2008, pp. 142-158. (doi: 10.2514/1.29958)

The curves, e.g., foil's upper and lower surfaces, are constructed via CST method. The multi-section surface is interpolated from several control sections.

## Installation

``` bash
pip install cst-modeling3d
```

## Example

### (1) Airfoil

```python
import numpy as np
from cst_modeling.foil import cst_foil

cst_u = np.array([ 0.1185,  0.1189,  0.1557,  0.1367,  0.2092,  0.1483,  0.1935])
cst_l = np.array([-0.1155, -0.1341, -0.1091, -0.2532, -0.0122, -0.1184,  0.0641])
x, yu, yl, t0, R0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)
```

<div align=center>
	<img src="example\airfoil\airfoil.png" width="400"> <br>
    Fig. A clean airfoil
</div>

### (2) Wing

```python
from cst_modeling.surface import Surface
wing = Surface(n_sec=6, name='Wing-tip', nn=101, ns=101)
wing.read_setting('Wing.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.05, 0.01])
wing.geo(split=False, showfoil=False)
wing.bend(4, 5, leader=[[21.0, 2.1, 30.0, 1.6]], kx=[0.6983, 4.0], ky=[0.1043, 1.10], rot_x=True)
wing.smooth(0,2)
wing.smooth(2,4)
wing.output_tecplot(fname='Wing-tip.dat', one_piece=False)
```

<div align=center>
	<img src="example\wing\wing-tip.jpg" width="500"> <br>
</div>


### (3) Blade

```python
from cst_modeling.surface import Surface
blade = Surface(n_sec=6, name='Blade-simple',nn=101, ns=51, projection=False)
blade.read_setting('Fan.txt', tail=[0.1, 0.1, 0.1, 0.1, 0.1, 0.05])
blade.geo()
blade.smooth(isec0=0, isec1=4)
blade.smooth(isec0=4, isec1=5, smooth0=True)
blade.surf_to_cylinder(flip=True)
blade.output_tecplot(fname='Blade-simple.dat')
```

<div align=center>
	<img src="example\blade\blade-simple-1.jpg" width="300"> <br>
</div>


### (4) Nacelle

<div align=center>
    <img src="example\nacelle\nacelle.jpg" width="280">
    <img src="example\nacelle\nacelle-frontview.jpg" width="280"> <br>
</div>

### (5) Fuselage

<div align=center>
    <img src="example\fuselage\fuselage.jpg" width="500"><br>
</div>

### (6) Fairing

<div align=center>
    <img src="example\fairing\fairing-fuselage.jpg" width="400"><br>
</div>

### (7) Delta wing

<div align=center>
    <img src="example\delta-wing\delta.jpg" width="400"><br>
</div>
