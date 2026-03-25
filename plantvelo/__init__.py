"""plantvelo
~~~~~~~~~~~
Plant-specific RNA velocity analysis.

Extends velocyto with a three-state splicing model:
    unspliced (U) → intron_retained (IR) → spliced (S)

Key additions over velocyto:
- ``PlantPermissive10X`` and ``PlantValidated10X`` logic classes that emit an
  ``intron_retained`` layer in addition to spliced / unspliced / ambiguous.
- ``plantvelo run`` CLI command with ``--ir-flanking`` parameter.
"""

from plantvelo._version import __version__
from plantvelo.logic import (
    Logic,            # re-export for convenience
    PlantPermissive10X,
    PlantValidated10X,
    PLANT_LOGICS,
)

__all__ = [
    "__version__",
    "Logic",
    "PlantPermissive10X",
    "PlantValidated10X",
    "PLANT_LOGICS",
]
