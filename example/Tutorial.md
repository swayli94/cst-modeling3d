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

cst_u = np.array([ 0.699328, 0.701191, 0.918286, 0.806252, 1.233955, 0.874498, 1.141528])
cst_l = np.array([-0.681143,-0.791295,-0.643583,-1.493055,-0.072058,-0.698530, 0.377973])
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
            (Meanwhile, the relative maximum thickness is kept the same)
```



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
>>> from cst_modeling.foil import output_foil
>>> output_foil(x, yu, yl, fname='airfoil.dat', ID=0, info=False)
```

Settings of this function:

```text
fname:  otuput file name
ID:     ID of this airfoil, also the zone name
        When ID = 0, outputs the file header
info:   if True, then outputs the curvature, thickness and camber distribution as well
```



### 2.4 airfoil modification





## 3. Wing





## 4. Fan blade



