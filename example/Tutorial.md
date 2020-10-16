# Tutorial

## 1. CST method

The class shape function transformation (CST) method combines class functions and shape functions to describe an arbitrary geometry and can guarantee airfoil smoothness with comparatively fewer design variables. Usually, a sixth-order Bernstein polynomial is used as the shape function; i.e., seven CST parameters are used to describe upper and lower surfaces.

Reference: Kulfan, B. M., “Universal parametric geometry representation method,” Journal of Aircraft, vol. 45, No. 1, 2008, pp. 142-158. (doi: 10.2514/1.29958)



<div align=center>
	<img src="airfoil\shape-function.png" width="300"> <br>
    Fig. CST shape functions (sixth-order)
</div>

Code:

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

The airfoil upper and lower surfaces are the linear combinations of these shape functions.



## 2. Airfoil

### 2.1 clean airfoils

Build a clean airfoil with given CST coefficients.









## 3. Wing





## 4. Fan blade



