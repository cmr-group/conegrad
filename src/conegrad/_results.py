"""Result dataclasses for the conegrad designer.

Kept in a separate module from the C-binding wrappers so that downstream code
can import the result types without pulling in the compiled extension.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ExpandedTrajectory:
    """Full per-shot gradient waveforms for a cones acquisition.

    Produced by :meth:`ConegradResult.expand_trajectory`. Each row is one
    interleaf (one TR's worth of gradient play), already polar-scaled and
    azimuthally rotated from the base cones in the parent result.

    Attributes
    ----------
    gx, gy, gz : np.ndarray[float32], shape (n_shots, ngradact)
        Per-shot gradient waveforms in physical units -- G/cm or T/m,
        matching the ``units`` mode of the parent ``ConegradResult``.
    band_index : np.ndarray[int32], shape (n_shots,)
        Polar band index (0 .. 2*ntheta-2) for each shot.
    interleaf_index : np.ndarray[int32], shape (n_shots,)
        Azimuthal interleaf index within the polar band
        (0 .. ``nintpc[band_index]-1``).
    traj_length : np.ndarray[int32], shape (n_shots,)
        Per-shot copy of the parent ConegradResult's ``traj_length`` for the
        base cone this shot was rotated from. ``gx[s, :traj_length[s]]``
        slices to the spiral.
    rampdown_end : np.ndarray[int32], shape (n_shots,)
        Per-shot copy of the parent ``rampdown_end``. ``gx[s, :rampdown_end[s]]``
        slices to the spiral plus its ramp-to-zero (excludes the rewinder).
    """

    gx: np.ndarray
    gy: np.ndarray
    gz: np.ndarray
    band_index: np.ndarray
    interleaf_index: np.ndarray
    traj_length: np.ndarray
    rampdown_end: np.ndarray


@dataclass
class ConegradResult:
    """Outputs of a cones trajectory design.

    Attributes
    ----------
    gx, gy, gz : np.ndarray[float32], shape (numcones, ngradact)
        Per-cone base gradient waveforms in physical units -- G/cm (if called
        with ``units="cgs"``) or T/m (if ``units="SI"``).
    ntheta : int
        Number of elevation angles (theta bands) above the equator.
        Total bands designed is ``2 * ntheta - 1`` (mirror-symmetric design).
    nintpc : np.ndarray[int32], shape (2*ntheta-1,)
        Number of interleaves required per theta band to sample the cone.
    rspthetas : np.ndarray[float32], shape (2*ntheta-1,)
        Elevation angle (radians) of each theta band.
    max_gradient : float
        Maximum gradient amplitude actually used in the design.
    snr_efficiency : float
        SNR efficiency of the design. Only populated when ``ktraj_out_flag=1``
        (which also writes k-space / density / nintpc files to disk); zero
        otherwise. Generally not needed and left disabled by default.
    ngradact : int
        Actual waveform length per cone (readout + rewinder). This is the
        max over all cones; per-cone values live in :attr:`traj_length` plus
        a rewinder tail.
    traj_length : np.ndarray[int32], shape (numcones,)
        Per-cone *end of spiral* in gradient samples -- the last sample of the
        actual k-space-traversing trajectory. At this point the gradient is
        still near peak. ``gx[i, :traj_length[i]]`` slices to just the spiral.
        Varies cone-to-cone because the cones designer's binary search lands a
        few samples short or over ``grad_points``.
    rampdown_end : np.ndarray[int32], shape (numcones,)
        Per-cone *end of ramp-to-zero* -- the first sample at or after the
        spiral end where the gradient on the slowest axis has finished
        slewing back to zero. The rewinder lobe begins at the next sample.
        ``gx[i, :rampdown_end[i]]`` includes the spiral plus its
        ramp-to-zero tail but excludes the rewinder.
    """

    gx: np.ndarray
    gy: np.ndarray
    gz: np.ndarray
    ntheta: int
    nintpc: np.ndarray
    rspthetas: np.ndarray
    max_gradient: float
    snr_efficiency: float
    ngradact: int
    traj_length: np.ndarray
    rampdown_end: np.ndarray

    def expand_trajectory(self, rot_flag: float = 1.0) -> ExpandedTrajectory:
        """Materialize every shot of the cones acquisition.

        The base waveforms in ``self.gx/gy/gz`` cover only the upper
        hemisphere (``numcones`` polar bands from 0 to pi/2). This method
        reproduces the EPIC-side replication that builds the full sphere:

        - For each of the ``2*ntheta-1`` polar bands at angle
          ``theta = rspthetas[t]``, picks base cone
          ``ci = ceil(|theta| / (pi/2) * numcones) - 1`` and scales:

          .. code-block:: text

              xyscale = cos(theta) / cos(theta_down)
              zscale  = sin(theta) / sin(theta_up)

          ``zscale`` is automatically negative below the equator, flipping z
          for the southern hemisphere. ``theta_down``/``theta_up`` are the
          band edges that base cone ``ci`` was designed between.

        - Within each polar band, fires ``nintpc[t]`` interleaves rotated
          around z by ``phi = k * rot_flag * 2*pi / nintpc[t]``:

          .. code-block:: text

              gx_shot = xyscale * (cos(phi)*gx_base - sin(phi)*gy_base)
              gy_shot = xyscale * (sin(phi)*gx_base + cos(phi)*gy_base)
              gz_shot = zscale  * gz_base

        Within each (polar band, interleaf), all three axes come from the
        *same* base cone index -- never mix ``gx[a]`` with ``gy[b]`` for
        ``a != b`` (the cone spiral shape is jointly designed across axes).

        Parameters
        ----------
        rot_flag : float
            +1 or -1; direction of azimuthal rotation. Pass the same value
            you used (or will use) for the ``conegrad`` call.

        Returns
        -------
        ExpandedTrajectory
        """
        numcones = int(self.gx.shape[0])
        n_bands = int(self.rspthetas.size)
        # 0-indexed base cone selection per polar band, mirroring the EPIC
        # mapping at mm4dflow.e:17050.
        conenum = np.ceil(
            np.abs(self.rspthetas, dtype=np.float64) / (np.pi / 2) * numcones
        ).astype(np.int32)
        conenum = np.clip(conenum, 1, numcones) - 1

        n_shots = int(self.nintpc.sum())
        gx_out = np.empty((n_shots, self.ngradact), dtype=np.float32)
        gy_out = np.empty((n_shots, self.ngradact), dtype=np.float32)
        gz_out = np.empty((n_shots, self.ngradact), dtype=np.float32)
        band_idx = np.empty(n_shots, dtype=np.int32)
        leaf_idx = np.empty(n_shots, dtype=np.int32)
        traj_len_out = np.empty(n_shots, dtype=np.int32)
        rampdown_out = np.empty(n_shots, dtype=np.int32)

        half_pi_per_cone = (np.pi / 2) / numcones
        shot = 0
        for t in range(n_bands):
            theta = float(self.rspthetas[t])
            ci = int(conenum[t])
            theta_down = ci * half_pi_per_cone
            theta_up = (ci + 1) * half_pi_per_cone
            xyscale = np.float32(np.cos(theta) / np.cos(theta_down))
            zscale = np.float32(np.sin(theta) / np.sin(theta_up))

            gx_b = self.gx[ci].astype(np.float32)
            gy_b = self.gy[ci].astype(np.float32)
            gz_b_scaled = self.gz[ci].astype(np.float32) * zscale

            npc = int(self.nintpc[t])
            phi = np.arange(npc, dtype=np.float64) * rot_flag * 2 * np.pi / npc
            c = np.cos(phi).astype(np.float32)[:, None]
            s = np.sin(phi).astype(np.float32)[:, None]

            sl = slice(shot, shot + npc)
            gx_out[sl] = xyscale * (c * gx_b - s * gy_b)
            gy_out[sl] = xyscale * (s * gx_b + c * gy_b)
            gz_out[sl] = gz_b_scaled  # broadcasts across npc rows
            band_idx[sl] = t
            leaf_idx[sl] = np.arange(npc, dtype=np.int32)
            # Per-shot timing boundaries are identical to the base cone's,
            # since rotation+scaling preserves the time axis.
            traj_len_out[sl] = int(self.traj_length[ci])
            rampdown_out[sl] = int(self.rampdown_end[ci])
            shot += npc

        return ExpandedTrajectory(
            gx=gx_out, gy=gy_out, gz=gz_out,
            band_index=band_idx, interleaf_index=leaf_idx,
            traj_length=traj_len_out, rampdown_end=rampdown_out,
        )
