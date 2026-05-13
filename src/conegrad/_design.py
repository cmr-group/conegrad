"""Designer-side Python wrappers around the C extension.

This module is responsible for unit conversion at the C boundary, output
buffer scaling, and the variable-density FOV table helper. It depends on
the compiled nanobind extension :mod:`conegrad._conegrad`.
"""
from __future__ import annotations

import numpy as np

from . import _conegrad as _ext
from ._results import ConegradResult

# Conversion factors between the C designer's native CGS units (G, cm) and SI
# (T, m). These are exact: 1 G = 1e-4 T and 1 cm = 1e-2 m, so 1 G/cm = 1e-2 T/m
# (and the same factor of 100 covers G/cm -> T/m, G/cm/s -> T/m/s, etc.).
_GRAD_CGS_PER_SI = 100.0   # multiply T/m by this to get G/cm
_FOV_CGS_PER_SI = 100.0    # multiply m by this to get cm
_RES_CGS_PER_SI = 1000.0   # multiply m by this to get mm
_KZ_CGS_PER_SI = 0.01      # multiply 1/m by this to get 1/cm (kz units)

# MAX_PG_WAMP in scones_design.c; the C scales each physical gradient by
# sysGMAX / 32766 before storing it as an int instruction amp.
_MAX_PG_WAMP = 32766.0

# Number of FOV-vs-kr samples expected by the C designer in the variable-density
# FOV table. The C source hardcodes 100 across three concatenated channels.
_VD_NKR = 100


def _as_float32(arr, name):
    a = np.ascontiguousarray(arr, dtype=np.float32)
    if a.ndim != 1:
        raise ValueError(f"{name} must be 1-D")
    return a


def _compute_rampdown_end(gx, gy, gz, traj_length, ngradact):
    """Per-cone index of the first sample at or after the spiral end where the
    slowest gradient axis has finished slewing back to zero.

    The C designer's rampdown writes values descending to exactly 0 at sample
    ``traj_length + nrampdowna_axis - 1`` for each axis, and ``nrampdowna``
    varies per axis (it depends on each axis's gradient amplitude at the spiral
    end). The rewinder lobe starts at the very next sample, so we find the
    first per-axis zero after ``traj_length`` and take the max across axes.
    """
    n_cones = gx.shape[0]
    out = np.empty(n_cones, dtype=np.int32)
    for i in range(n_cones):
        ts = int(traj_length[i])
        latest = ts
        for axis in (gx[i, ts:], gy[i, ts:], gz[i, ts:]):
            zeros = np.flatnonzero(axis == 0)
            if zeros.size:
                latest = max(latest, ts + int(zeros[0]))
            else:
                # No zero found -- treat the end of the buffer as the boundary.
                latest = max(latest, ngradact)
        out[i] = latest
    return out


def _call(cfunc, *, res, fov, expected_fov_len,
          numcones, grad_points, gsp, read_points, tsp,
          precision, dcf, oversample,
          smax, max_s_rewind, gmax, sys_gmax, sys_max_rewind_g,
          mindens, rewind_flag,
          t_fracx, t_fracy, t_fracz,
          acq_mode, rot_flag, sym_flag, slab_kz,
          output_grad, ktraj_flag, ktraj_out_flag, endian_flag,
          rhkacq_uid, cones_plot_flag, verbose,
          output_grad_scale):
    res = _as_float32(res, "res")
    fov = _as_float32(fov, "fov")
    if res.size != 2:
        raise ValueError("res must have exactly 2 elements: [xy_res, z_res]")
    if fov.size != expected_fov_len:
        raise ValueError(
            f"fov must have exactly {expected_fov_len} elements "
            f"(got {fov.size})"
        )

    out = cfunc(
        res, fov,
        int(numcones), int(grad_points), float(gsp),
        int(read_points), float(tsp),
        float(precision), float(dcf), int(oversample),
        float(smax), float(max_s_rewind),
        float(gmax), float(sys_gmax), float(sys_max_rewind_g),
        int(output_grad), int(ktraj_flag), int(ktraj_out_flag),
        int(endian_flag), int(rewind_flag), float(mindens),
        float(t_fracx), float(t_fracy), float(t_fracz),
        int(acq_mode), float(rot_flag), int(sym_flag),
        float(slab_kz), int(rhkacq_uid), int(cones_plot_flag),
        bool(verbose),
    )
    (gx, gy, gz, ntheta, nintpc, rspthetas, max_gradient, snr_eff,
     ngradact, traj_length) = out
    rampdown_end = _compute_rampdown_end(gx, gy, gz, traj_length, ngradact)
    # Convert the int32 instruction amplitudes back to physical gradient units.
    # The C designer scales each float gradient by sysGMAX/MAX_PG_WAMP=32766; we
    # invert that here using the caller's user-units sys_gmax so the returned
    # waveforms are in G/cm (cgs mode) or T/m (SI mode).
    scale = np.float32(output_grad_scale)
    return ConegradResult(
        gx=gx.astype(np.float32) * scale,
        gy=gy.astype(np.float32) * scale,
        gz=gz.astype(np.float32) * scale,
        ntheta=int(ntheta),
        nintpc=nintpc,
        rspthetas=rspthetas,
        max_gradient=float(max_gradient),
        snr_efficiency=float(snr_eff),
        ngradact=int(ngradact),
        traj_length=traj_length,
        rampdown_end=rampdown_end,
    )


def _validate_units(units):
    if units not in ("cgs", "SI"):
        raise ValueError(f"units must be 'cgs' or 'SI', got {units!r}")


def conegrad(
    res,
    fov,
    *,
    numcones: int = 32,
    grad_points: int,
    gsp: float,
    read_points: int | None = None,
    tsp: float | None = None,
    smax: float,
    gmax: float,
    max_s_rewind: float | None = None,
    sys_gmax: float | None = None,
    sys_max_rewind_g: float | None = None,
    precision: float = 0.5,
    dcf: float = 0.0,
    oversample: int = 1,
    mindens: float = 0.0,
    rewind_flag: int = 1,
    t_fracx: float = 0.0,
    t_fracy: float = 0.0,
    t_fracz: float = 0.0,
    acq_mode: int = 1,
    rot_flag: float = 1.0,
    sym_flag: int = 1,
    slab_kz: float = 0.0,
    output_grad: int = 0,
    ktraj_flag: int = 1,
    ktraj_out_flag: int = 0,
    endian_flag: int = 0,
    rhkacq_uid: int = 0,
    cones_plot_flag: int = 0,
    verbose: bool = False,
    units: str = "cgs",
) -> ConegradResult:
    """Design a constant-FOV 3D Cones gradient trajectory.

    Parameters
    ----------
    res : array-like, shape (2,)
        Spatial resolution ``[xy, z]``, in mm (cgs) or m (SI), each scaled by
        the gyromagnetic ratio relative to proton.
    fov : array-like, shape (2,)
        Field of view ``[xy, z]``, in cm (cgs) or m (SI).
    numcones : int
        Number of cone bands to design (1..100). Default 32.
    grad_points : int
        Number of gradient samples during the readout window.
    gsp : float
        Gradient sample period in seconds (e.g. 4e-6 for 4 us).
    read_points : int, optional
        Number of (oversampled) readout sample points. Only used when
        ``ktraj_out_flag=1`` (writes the k-space trajectory file at the ADC
        sample rate). If ``None`` (default), set equal to ``grad_points``.
    tsp : float, optional
        Readout sample period in seconds. Only used when ``ktraj_out_flag=1``.
        If ``None`` (default), set equal to ``gsp``.
    smax : float
        Maximum slew rate during readout, in G/cm/s (cgs) or T/m/s (SI).
    gmax : float
        Maximum gradient amplitude during readout, in G/cm (cgs) or T/m (SI).
    max_s_rewind : float, optional
        Maximum slew rate during rewinder. Defaults to ``smax / sqrt(3)`` --
        the multi-axis scaling that ensures the vector slew never exceeds
        ``smax`` when all three axes ramp simultaneously during the rewinder.
    sys_max_rewind_g : float, optional
        Maximum gradient amplitude during rewinder. Defaults to
        ``gmax / sqrt(3)`` for the same multi-axis reason.
    sys_gmax : float, optional
        Hardware gradient ceiling, used only to scale the C designer's
        internal integer instruction amplitudes (the wrapper inverts this on
        the way out so it has no user-visible effect). Defaults to ``gmax``.
    precision : float
        Cones design precision/tolerance (0..1). Default 0.5.
    dcf : float
        Density compensation factor (0..1). Default 0.
    oversample : int
        Trajectory oversampling factor. Default 1.
    mindens : float
        Minimum sampling density (0..10). Default 0.
    rewind_flag : int
        Rewinder strategy (0..4). Default 1.
    t_fracx, t_fracy, t_fracz : float
        Per-axis fractional gradient delays (timing calibration).
    acq_mode : int
        1 = Cones, 2 = Projection Reconstruction. Default 1.
    rot_flag : float
        Rotation direction between interleaves (+/- 1). Default 1.
    sym_flag : int
        1 = k-space symmetric, 2 = asymmetric. Default 1.
    slab_kz : float
        kz phase across the slab. Default 0.
    output_grad, ktraj_flag, ktraj_out_flag, endian_flag : int
        File-output flags (write waveforms / k-space trajectory to disk in the
        current working directory). All default to 0 in the Python wrapper.
    rhkacq_uid : int
        UID used in output filenames when file-output flags are set.
    cones_plot_flag : int
        Gnuplot-driven plot output (0..3). Default 0.
    verbose : bool
        If True, the C code's debug printf output is sent to stdout.
    units : {"cgs", "SI"}
        Unit system for the dimensioned inputs and outputs. ``"cgs"`` (default)
        keeps the C designer's native units (G/cm, G/cm/s, cm, mm). ``"SI"``
        accepts T/m, T/m/s, m, and converts at the boundary; ``result.max_gradient``
        and the gradient arrays are returned in T/m.

    Returns
    -------
    ConegradResult
    """
    _validate_units(units)
    if read_points is None:
        read_points = grad_points
    if tsp is None:
        tsp = gsp
    # Derive rewinder limits and hardware ceiling from the readout limits if
    # not supplied. The 1/sqrt(3) multi-axis factor matches the EPIC default
    # at mm4dflow.e:3815-3828; sys_gmax falls back to gmax (its only effect
    # is the internal int-amp quantization step, cancelled on output).
    if max_s_rewind is None:
        max_s_rewind = smax / np.sqrt(3)
    if sys_max_rewind_g is None:
        sys_max_rewind_g = gmax / np.sqrt(3)
    if sys_gmax is None:
        sys_gmax = gmax
    sys_gmax_user = float(sys_gmax)
    if units == "SI":
        res = np.asarray(res, dtype=np.float32) * _RES_CGS_PER_SI
        fov = np.asarray(fov, dtype=np.float32) * _FOV_CGS_PER_SI
        gmax = gmax * _GRAD_CGS_PER_SI
        sys_gmax = sys_gmax * _GRAD_CGS_PER_SI
        sys_max_rewind_g = sys_max_rewind_g * _GRAD_CGS_PER_SI
        smax = smax * _GRAD_CGS_PER_SI
        max_s_rewind = max_s_rewind * _GRAD_CGS_PER_SI
        slab_kz = slab_kz * _KZ_CGS_PER_SI

    result = _call(
        _ext.conegrad,
        res=res, fov=fov, expected_fov_len=2,
        numcones=numcones, grad_points=grad_points, gsp=gsp,
        read_points=read_points, tsp=tsp,
        precision=precision, dcf=dcf, oversample=oversample,
        smax=smax, max_s_rewind=max_s_rewind,
        gmax=gmax, sys_gmax=sys_gmax, sys_max_rewind_g=sys_max_rewind_g,
        mindens=mindens, rewind_flag=rewind_flag,
        t_fracx=t_fracx, t_fracy=t_fracy, t_fracz=t_fracz,
        acq_mode=acq_mode, rot_flag=rot_flag, sym_flag=sym_flag,
        slab_kz=slab_kz,
        output_grad=output_grad, ktraj_flag=ktraj_flag,
        ktraj_out_flag=ktraj_out_flag, endian_flag=endian_flag,
        rhkacq_uid=rhkacq_uid, cones_plot_flag=cones_plot_flag,
        verbose=verbose,
        output_grad_scale=sys_gmax_user / _MAX_PG_WAMP,
    )
    if units == "SI":
        result.max_gradient = result.max_gradient / _GRAD_CGS_PER_SI
    return result


def conegrad_vd(
    res,
    fov,
    *,
    numcones: int = 32,
    grad_points: int,
    gsp: float,
    read_points: int | None = None,
    tsp: float | None = None,
    smax: float,
    gmax: float,
    max_s_rewind: float | None = None,
    sys_gmax: float | None = None,
    sys_max_rewind_g: float | None = None,
    precision: float = 0.5,
    dcf: float = 0.0,
    oversample: int = 1,
    mindens: float = 0.0,
    rewind_flag: int = 1,
    t_fracx: float = 0.0,
    t_fracy: float = 0.0,
    t_fracz: float = 0.0,
    acq_mode: int = 1,
    rot_flag: float = 1.0,
    sym_flag: int = 1,
    slab_kz: float = 0.0,
    output_grad: int = 0,
    ktraj_flag: int = 1,
    ktraj_out_flag: int = 0,
    endian_flag: int = 0,
    rhkacq_uid: int = 0,
    cones_plot_flag: int = 0,
    verbose: bool = False,
    units: str = "cgs",
) -> ConegradResult:
    """Design a variable-density 3D Cones gradient trajectory.

    Identical to :func:`conegrad` except that ``fov`` is a length-300 array
    laid out as three concatenated channels of 100 samples each:

    - ``fov[0:100]``    -- FOV_xy as a function of normalized kr
    - ``fov[100:200]``  -- FOV_z  as a function of normalized kr
    - ``fov[200:300]``  -- normalized kr (in [0, 1], unitless in both unit systems)

    Use :func:`make_vd_fov` to build this table with the standard EPIC profile.
    See :func:`conegrad` for full parameter documentation, including the
    ``units`` kwarg.
    """
    _validate_units(units)
    if read_points is None:
        read_points = grad_points
    if tsp is None:
        tsp = gsp
    if max_s_rewind is None:
        max_s_rewind = smax / np.sqrt(3)
    if sys_max_rewind_g is None:
        sys_max_rewind_g = gmax / np.sqrt(3)
    if sys_gmax is None:
        sys_gmax = gmax
    sys_gmax_user = float(sys_gmax)
    if units == "SI":
        res = np.asarray(res, dtype=np.float32) * _RES_CGS_PER_SI
        # FOV layout: first two 100-sample channels are FOV_xy and FOV_z (m
        # in SI); the third channel is normalized kr and stays unitless.
        fov = np.asarray(fov, dtype=np.float32).copy()
        fov[0:200] *= _FOV_CGS_PER_SI
        gmax = gmax * _GRAD_CGS_PER_SI
        sys_gmax = sys_gmax * _GRAD_CGS_PER_SI
        sys_max_rewind_g = sys_max_rewind_g * _GRAD_CGS_PER_SI
        smax = smax * _GRAD_CGS_PER_SI
        max_s_rewind = max_s_rewind * _GRAD_CGS_PER_SI
        slab_kz = slab_kz * _KZ_CGS_PER_SI

    result = _call(
        _ext.conegrad_vd,
        res=res, fov=fov, expected_fov_len=300,
        numcones=numcones, grad_points=grad_points, gsp=gsp,
        read_points=read_points, tsp=tsp,
        precision=precision, dcf=dcf, oversample=oversample,
        smax=smax, max_s_rewind=max_s_rewind,
        gmax=gmax, sys_gmax=sys_gmax, sys_max_rewind_g=sys_max_rewind_g,
        mindens=mindens, rewind_flag=rewind_flag,
        t_fracx=t_fracx, t_fracy=t_fracy, t_fracz=t_fracz,
        acq_mode=acq_mode, rot_flag=rot_flag, sym_flag=sym_flag,
        slab_kz=slab_kz,
        output_grad=output_grad, ktraj_flag=ktraj_flag,
        ktraj_out_flag=ktraj_out_flag, endian_flag=endian_flag,
        rhkacq_uid=rhkacq_uid, cones_plot_flag=cones_plot_flag,
        verbose=verbose,
        output_grad_scale=sys_gmax_user / _MAX_PG_WAMP,
    )
    if units == "SI":
        result.max_gradient = result.max_gradient / _GRAD_CGS_PER_SI
    return result


def make_vd_fov(fov_full, fov_under=0.5, n_full=20):
    """Build the 300-element FOV table that :func:`conegrad_vd` expects.

    Layout matches the EPIC convention at ``mm4dflow.e:3792-3803``: a flat
    ``inner | outer`` profile where the k-space center is sampled at the full
    FOV (Nyquist) and the periphery ramps linearly down to ``fov_under * fov_full``.

    Parameters
    ----------
    fov_full : array-like, shape (2,)
        ``[FOV_xy, FOV_z]`` -- the true (Nyquist) field of view. Same units as
        you pass to :func:`conegrad_vd` (cm in cgs mode, m in SI mode).
    fov_under : float or array-like, shape (2,)
        Undersampling factor(s). Scalar applies to both axes; pass ``[ux, uz]``
        for asymmetric undersampling. ``1.0`` = fully sampled; ``0.5`` = 2x
        undersampled at the edge of k-space, etc. Default ``0.5``.
    n_full : int
        Number of inner samples (out of 100) held at the full FOV before the
        ramp begins. Default ``20`` matches the EPIC reference.

    Returns
    -------
    np.ndarray, shape (300,), float32
        Concatenated ``[FOV_xy(kr) | FOV_z(kr) | kr]`` with 100 samples each.
        ``kr`` is normalized to ``[0, 1]``.
    """
    fov_full = np.asarray(fov_full, dtype=np.float32)
    if fov_full.size != 2:
        raise ValueError("fov_full must have 2 elements: [FOV_xy, FOV_z]")
    fov_under_arr = np.broadcast_to(
        np.asarray(fov_under, dtype=np.float32), (2,)
    ).astype(np.float32)
    if not (0 <= n_full <= _VD_NKR):
        raise ValueError(f"n_full must be in [0, {_VD_NKR}], got {n_full}")

    n_ramp = _VD_NKR - n_full
    fov_edge = fov_full * fov_under_arr  # FOV at kr=1 (most undersampled)

    out = np.empty(3 * _VD_NKR, dtype=np.float32)
    for ch, (full, edge) in enumerate(zip(fov_full, fov_edge)):
        base = ch * _VD_NKR
        out[base : base + n_full] = full
        if n_ramp > 0:
            # linspace(full, edge, n_ramp+1)[1:] reproduces the C arithmetic:
            #   FOV[n_full + i] = full + (edge - full) * (i + 1) / n_ramp.
            out[base + n_full : base + _VD_NKR] = np.linspace(
                full, edge, n_ramp + 1, dtype=np.float32
            )[1:]
    out[2 * _VD_NKR : 3 * _VD_NKR] = np.linspace(0, 1, _VD_NKR, dtype=np.float32)
    return out
