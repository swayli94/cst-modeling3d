Introduction
=====================

CST method
---------------------------

Surface and curve modeling via class shape function transformation (CST) method.

The CST method combines class functions and shape functions to describe an arbitrary geometry 
and can guarantee airfoil smoothness with comparatively fewer design variables. 
Usually, a sixth-order Bernstein polynomial is used as the shape function; i.e., seven CST parameters are used to describe upper and lower surfaces.

.. seealso::

    Kulfan, B. M., “Universal parametric geometry representation method,” Journal of Aircraft, vol. 45, No. 1, 2008, pp. 142-158. (doi: 10.2514/1.29958)

The curves, e.g., foil's upper and lower surfaces, are constructed via CST method. The multi-section surface is interpolated from several control sections.


Project structure
---------------------------

The following is an overview of the project's structure:

.. code-block:: text

    cst-modeling3d/
    ├── cst_modeling/
    │   ├── __init__.py
    │   ├── basic.py
    │   ├── section.py
    │   ├── surface.py
    │   ├── foil.py
    │   └── tools/
    │       ├── __init__.py
    │       ├── naca.py
    │       ├── xfoil.py
    │       ├── blwf.py
    │       └── auxiliary.py
    ├── docs/
    ├── example/
    ├── .gitignore
    ├── readthedocs.yaml
    ├── setup.py
    ├── LICENSE
    ├── README.md
    └── requirements.txt


Directories and Files
---------------------

- **cst_modeling/**: Contains the source code.

  - **__init__.py**: Initialization file for the package.
  - **basic.py**: Basic classes for sections and surfaces, and fundamental functions.
  - **section.py**: Classes and functions for construction two-dimensional sections.
  - **surface.py**: Classes and functions for construction three-dimensional surfaces.
  - **foil.py**: Classes and functions for airfoil geometric feature extraction and modification.
  - **tools/**: Contains auxiliary tools.
  
    - **__init__.py**: Initialization file for the tools package.
    - **naca.py**: Functions for NACA airfoil generation.
    - **xfoil.py**: Functions for interfacing with XFOIL.
    - **blwf.py**: Functions for interfacing with BLWF.
    - **auxiliary.py**: Auxiliary functions.

- **docs/**: Contains the Sphinx documentation source files.
- **example/**: Contains example scripts.
- **.gitignore**: Git ignore file specifying files and directories to ignore.
- **readthedocs.yaml**: Configuration file for Read the Docs.
- **requirements.txt**: File specifying project dependencies.
- **setup.py**: Setup script for the project.
- **LICENSE**: License file.
- **README.md**: Project README file.


Installation
---------------------------

To install the package, run the following command:

.. code-block:: bash

    # Install the package from PyPI
    pip install cst-modeling3d

    # Install the package from the source code
    git clone https://github.com/swayli94/cst-modeling3d
    cd cst-modeling3d
    pip install -e .

Please note that the package in the PyPI repository may not be the latest version.

