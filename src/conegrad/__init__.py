"""Python wrapper for the conegrad / conegrad_vd 3D Cones MRI trajectory designer.

The underlying C code (originally part of a GE EPIC pulse sequence) takes
scan-time parameters (resolution, FOV, max gradient / slew, etc.) and produces
gradient waveforms on the x, y, z axes for a family of conical interleaves
that together tile k-space.

Two entry points are exposed:

- :func:`conegrad`     - constant FOV design.
- :func:`conegrad_vd`  - variable-density FOV (FOV varies along the k-radius).

Both return a :class:`ConegradResult` with numpy arrays for the gradient
waveforms and ancillary outputs. Use :meth:`ConegradResult.expand_trajectory`
to materialize every shot of the acquisition from the per-cone base waveforms.
:func:`make_vd_fov` builds the variable-density FOV table for ``conegrad_vd``.
"""
from __future__ import annotations

from ._design import conegrad, conegrad_vd, make_vd_fov
from ._results import ConegradResult, ExpandedTrajectory

__all__ = [
    "ConegradResult",
    "ExpandedTrajectory",
    "conegrad",
    "conegrad_vd",
    "make_vd_fov",
]
