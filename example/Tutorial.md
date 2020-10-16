# Tutorial

## 1. CST method

The class shape function transformation (CST) method combines class functions and shape functions to describe an arbitrary geometry and can guarantee airfoil smoothness with comparatively fewer design variables. Usually, a sixth-order Bernstein polynomial is used as the shape function; i.e., seven CST parameters are used to describe upper and lower surfaces.

Reference: Kulfan, B. M., “Universal parametric geometry representation method,” Journal of Aircraft, vol. 45, No. 1, 2008, pp. 142-158. (doi: 10.2514/1.29958)



<div align=center>
	<img src="airfoil\shape-function.png" width="300"> <br>
    Fig. CST shape functions (sixth-order)
</div>

```python
from cst_modeling.foil import cst_curve
from matplotlib import pyplot as plt

n_cst = 7
plt.figure()

for i in range(n_cst):
    cst = np.zeros(n_cst)
    cst[i] = 1.0
    x, y = cst_curve(101, cst)
    plt.plot(x, y)

plt.show()
```

The airfoil upper and lower surfaces are the **linear combinations** of these shape functions.



## 2. Airfoil

### 2.1 clean airfoils

Build a clean airfoil with given CST coefficients.

```python
import numpy as np
from cst_modeling.foil import cst_foil

cst_u = np.array([ 0.118598,  0.118914,  0.155731,  0.136732,  0.209265,  0.148305,  0.193591])
cst_l = np.array([-0.115514, -0.134195, -0.109145, -0.253206, -0.012220, -0.118463,  0.064100])
x, yu, yl, t0, R0 = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)

plt.figure()
plt.plot(x, yu, 'b')
plt.plot(x, yl, 'b')
plt.show()
```

<div align=center>
	<img src="airfoil\airfoil.png" width="300"> <br>
    Fig. A clean airfoil
</div>

Settings of this function:

```text
>>> cst_foil(nn, coef_upp, coef_low, x=None, t=None, tail=0.0)
nn:         total amount of points on the upper/lower surfaces
x:          optional ndarray [nn]. If given, the points on the airfoil are placed by the given x.
t:          optional float. If given, the airfoil y is scaled to mathc the given relative maximum thickness
tail:       optional float. If given, the airfoil is stretched to have the given relative tail thickness
            Meanwhile, the relative maximum thickness is kept unchanged if t is specified.
            Otherwise, the thickness will increase, when the tail is added.
```

<div align=center>
	<img src="airfoil\airfoil-tail-1.png" width="300"> <img src="airfoil\airfoil-tail-2.png" width="300"> <br>
    Fig. Adding tail to an airfoil (left: t=None, right: t=0.11)
</div>



### 2.2 CST coefficients of certain airfoil

Find the CST coefficients of a given airfoil by least square method.

This function allows the airfoil has non-zero tail thickness.

Also allows the airfoil chord length not equals to one.

```pytho
from cst_modeling.foil import cst_foil_fit
cst_u, cst_l = cst_foil_fit(x, yu, x, yl, n_order=n_cst)
# The cst_u, cst_l are ndarrays [n_cst]
```



### 2.3 output airfoils

Output the curves of airfoil upper and lower surfaces in Tecplot format.

```python
from cst_modeling.foil import output_foil
output_foil(x, yu, yl, fname='airfoil.dat', ID=0, info=False)
```

Settings of this function:

```text
fname:  otuput file name
ID:     ID of this airfoil, also the zone name
        When ID = 0, outputs the file header
info:   if True, then outputs the curvature, thickness and camber distribution as well
```



### 2.4 airfoil modification

Add bump function to the airfoil.

```python
from cst_modeling.foil import foil_bump_modify
xc = 0.7
h  = 0.02
s  = 0.6
yu_, yl_ = foil_bump_modify(x, yu, yl, xc, h, s, side=1)
```

Settings of this function:

```text
xc:         x location of the bump center
h:          bump height (can be either positive or negetive)
s:          bump width
side:       1/-1, modification to the upper/lower surface
n_order:    if specified (>0), then use CST to fit the new foil
```

<div align=center>
	<img src="airfoil\airfoil-bump.png" width="300"> <br>
    Fig. Bump modification to airfoil upper surface
</div>



### 2.5 airfoil increment

Add an incremental CST curve to the airfoil.

```python
from cst_modeling.foil import foil_increment
cst_u_ = np.zeros(16)
cst_l_ = np.zeros(16)
cst_u_[12] = 0.05
yu_, yl_ = foil_increment(x, yu, yl, cst_u_, cst_l_, t=t0)
```

Settings of this function:

```text
cst_u_:     CST coefficients of the incremental curve on the upper surface
cst_l_:     CST coefficients of the incremental curve on the lower surface
t:          optional float. If given, the airfoil relative maximum thickness is kept unchanged
```

<div align=center>
	<img src="airfoil\airfoil-increment.png" width="300"> <br>
    Fig. Increment to airfoil upper surface
</div>

The green dash curve is the incremental CST curve.



## 3. Wing





## 4. Fan blade



