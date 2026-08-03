'''
CST modeling of airfoils and 3D surfaces.

Information of modeling
=======================

1. Coordinates

    The flow direction is in +X direction;
    The symmetry plane is XY plane, the up direction is +Y;
    The span-wise direction is +Z;

2. Curve orientation

    The 2D airfoil is built by the upper and lower surfaces in the XY plane,
    i.e., xu, yu, xl, yl.

    The 3D airfoil is a single curve combining the upper and lower surfaces,
    i.e., x, y, z, the curve starts from the trailing edge of the lower surface.

    Without 3D rotation, the closed 3D curve is in the clockwise direction in XY plane;
    Without 3D rotation, the open 3D curve is in the +X direction y=f(x),
    +Y direction z=f(y), or +Z direction x=f(z);
    In other words, the direction of curves is based on the left-hand rule.
'''

__version__ = "0.2.8"

#* Public API (new lofting-based architecture).
#* The classic-architecture `Surface`/`OpenSurface` remain importable via `cst_modeling.basic`.

from .surface import Surface, OpenSurface, BasicSurface
from .operation import GuideCurve, Lofting, Lofting_2Profile, Lofting_Revolution
from .section import BasicSection, Section, OpenSection, RoundTipSection, cst_foil, cst_foil_fit
from .math import cst_curve, dist_clustcos
from .foil import FoilGeoFeatures, FoilModification

__all__ = [
    'Surface', 'OpenSurface', 'BasicSurface',
    'GuideCurve', 'Lofting', 'Lofting_2Profile', 'Lofting_Revolution',
    'BasicSection', 'Section', 'OpenSection', 'RoundTipSection', 'cst_foil', 'cst_foil_fit',
    'cst_curve', 'dist_clustcos',
    'FoilGeoFeatures', 'FoilModification',
]
