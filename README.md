# CST Modeling

Surface and curve modeling via class shape function transformation (CST) method.

The CST method combines class functions and shape functions to describe an arbitrary geometry and can guarantee airfoil smoothness with comparatively fewer design variables. Usually, a sixth-order Bernstein polynomial is used as the shape function; i.e., seven CST parameters are used to describe upper and lower surfaces.

Reference: Kulfan, B. M., “Universal parametric geometry representation method,” Journal of Aircraft, vol. 45, No. 1, 2008, pp. 142-158. (doi: 10.2514/1.29958)

The curves, e.g., foil's upper and lower surfaces, are constructed via CST method. The multi-section surface is interpolated from several control sections.

## Installation

``` bash
# Install the package from PyPI
pip install cst-modeling3d

# Install the package from the source code
git clone https://github.com/swayli94/cst-modeling3d
cd cst-modeling3d
pip install -e .
```

## Tutorial

https://cst-modeling3d.readthedocs.io/en/latest/

## Example

### (1) Airfoil

```python
import numpy as np
from cst_modeling.section import cst_foil

cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])

x1, yu1, yl1, tmax1, rLE1 = cst_foil(201, cst_u, cst_l, x=None, t=None, tail=0.0)
x2, yu2, yl2, tmax2, rLE2 = cst_foil(201, cst_u, cst_l, x=None, t=None, tail=0.004)
x3, yu3, yl3, tmax3, rLE3 = cst_foil(201, cst_u, cst_l, x=None, t=0.11, tail=0.004)
```

<div align=center>
	<img src="example\airfoil-basic\cst-airfoil-add-tail.png" width="300">
    <img src="example\airfoil-basic\cst-airfoil-add-tail-same-tmax.png" width="300"> <br>
</div>


### (2) Wing

```python
from cst_modeling.surface2 import Surface

wing = Surface(n_sec=10, name='Wing-CRM-winglet', nn=201, ns=51, 
                smooth_surface=True, smooth_sections=[(0, 2), (4, 7), (8, 9)],
                rotate_x_section=True, rotation_sections=[(8, 9)])

wing.read_setting('Wing.txt')
wing.prepare()
wing.geo()
```

<div align=center>
	<img src="example\winglet\wing-crm-winglet-rotx.png" width="500"> <br>
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
    <img src="example\nacelle\nacelle-axisymmetric.png" width="500"> <br>
    <img src="example\nacelle\nacelle-non-axisymmetric-2.png" width="500"> <br>
</div>

### (5) Fuselage

<div align=center>
    <img src="example\fuselage\fuselage.jpg" width="400"><br>
</div>

### (6) Fairing

<div align=center>
    <img src="example\fairing\fairing-fuselage.jpg" width="400"><br>
</div>

### (7) Delta wing

<div align=center>
    <img src="example\delta-wing\delta.jpg" width="400"><br>
</div>
