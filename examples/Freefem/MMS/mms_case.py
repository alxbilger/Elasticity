"""Cross-dim base class for manufactured-solution cases.

Dim-specific bases (`MMSCase1D`, `MMSCase2D`, `MMSCase3D`) inherit from
`MMSCase` and declare the abstract solution/derivative/source/BC methods
plus their `source_quadrature_*` placeholders. The case files
(`cubic.py`, `trigonometric.py`, `sinus_neumann.py`, …) inherit from the
dim-specific bases and override the abstract methods explicitly.
"""

from abc import ABC


class MMSCase(ABC):
    name       = None  # case identifier (must match the params.json key)
    plot_label = None  # LaTeX label for the exact solution
