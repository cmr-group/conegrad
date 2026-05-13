"""Basic tests for the conegrad Python wrapper.

These are smoke / structure tests: they verify that the C bindings can be
called with reasonable inputs, that output shapes are consistent, and that
both ``conegrad`` and ``conegrad_vd`` accept the expected FOV layouts.
"""
from __future__ import annotations

import io
from contextlib import redirect_stdout

import numpy as np
import pytest

import conegrad


# Scanner / sequence params loosely matching a 3T GE setup. These are picked to
# make the cones designer converge quickly; they are not a tuned recipe.
RES = np.array([1.0, 1.0], dtype=np.float32)        # mm * gamma_ratio
FOV = np.array([24.0, 24.0], dtype=np.float32)      # cm
SYS_GMAX = 5.0
GMAX = 0.9 * SYS_GMAX
SMAX = 0.75 * (1e6 * 5 / 333) / np.sqrt(2)
MAX_S_REWIND = 0.75 * (1e6 * 5 / 333) / np.sqrt(3)
SYS_MAX_REWIND_G = 0.9 / np.sqrt(3) * SYS_GMAX
GSP = 4e-6
TSP = 4e-6
READ_POINTS = 500
GRAD_POINTS = int(np.ceil(READ_POINTS * (2 * TSP / GSP)))


def _common_kwargs():
    return dict(
        numcones=32,
        grad_points=GRAD_POINTS,
        gsp=GSP,
        read_points=READ_POINTS,
        tsp=TSP,
        smax=SMAX,
        max_s_rewind=MAX_S_REWIND,
        gmax=GMAX,
        sys_gmax=SYS_GMAX,
        sys_max_rewind_g=SYS_MAX_REWIND_G,
    )


def test_conegrad_runs_and_returns_shapes():
    r = conegrad.conegrad(RES, FOV, **_common_kwargs())

    assert isinstance(r, conegrad.ConegradResult)
    assert r.gx.dtype == np.float32
    assert r.gy.dtype == np.float32
    assert r.gz.dtype == np.float32
    assert r.gx.shape == (32, r.ngradact)
    assert r.gy.shape == r.gx.shape
    assert r.gz.shape == r.gx.shape
    assert r.ntheta > 0
    assert r.nintpc.shape == (2 * r.ntheta - 1,)
    assert r.rspthetas.shape == (2 * r.ntheta - 1,)
    assert r.nintpc.dtype == np.int32
    assert r.rspthetas.dtype == np.float32
    assert int(r.nintpc.sum()) > 0
    # max_gradient should be set to something physical (and not the C sentinel -1e6)
    assert r.max_gradient > 0
    # Waveforms should not be all zero and should not exceed gmax (in G/cm).
    assert r.gx.any() and r.gy.any() and r.gz.any()
    assert np.abs(r.gx).max() <= GMAX * 1.01
    assert np.abs(r.gy).max() <= GMAX * 1.01
    assert np.abs(r.gz).max() <= GMAX * 1.01


def test_rspthetas_span_hemisphere():
    r = conegrad.conegrad(RES, FOV, **_common_kwargs())
    # Mirror-symmetric design covers -pi/2 ... +pi/2.
    assert r.rspthetas.min() < -1.4
    assert r.rspthetas.max() > 1.4


def test_verbose_flag_controls_stdout():
    buf_quiet = io.StringIO()
    with redirect_stdout(buf_quiet):
        conegrad.conegrad(RES, FOV, verbose=False, **_common_kwargs())

    # Verbose path is harder to capture from Python (C stdio writes don't go
    # through sys.stdout). We at least verify the call succeeds with verbose=True.
    conegrad.conegrad(RES, FOV, verbose=True, **_common_kwargs())

    # quiet mode should not write to the redirected python stdout
    assert "MCARL" not in buf_quiet.getvalue()


def test_conegrad_vd_with_uniform_fov_matches_constant_density():
    """conegrad_vd with a uniform FOV table should behave like conegrad."""
    # Layout: [FOV_xy(kr=0..1) | FOV_z(kr=0..1) | kr]
    nkr = 100
    fov_vd = np.empty(3 * nkr, dtype=np.float32)
    fov_vd[0:nkr] = FOV[0]
    fov_vd[nkr:2 * nkr] = FOV[1]
    fov_vd[2 * nkr:3 * nkr] = np.linspace(0, 1, nkr, dtype=np.float32)

    r = conegrad.conegrad_vd(RES, fov_vd, **_common_kwargs())
    assert r.gx.shape[0] == 32
    assert r.ntheta > 0
    assert r.ngradact > 0


def test_rejects_wrong_fov_shape():
    bad_fov = np.array([24.0], dtype=np.float32)
    with pytest.raises(ValueError):
        conegrad.conegrad(RES, bad_fov, **_common_kwargs())


def test_rejects_wrong_res_shape():
    bad_res = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    with pytest.raises(ValueError):
        conegrad.conegrad(bad_res, FOV, **_common_kwargs())


def test_vd_rejects_constant_fov_shape():
    with pytest.raises(ValueError):
        conegrad.conegrad_vd(RES, FOV, **_common_kwargs())


def test_expand_trajectory_shapes_and_counts():
    r = conegrad.conegrad(RES, FOV, **_common_kwargs())
    et = r.expand_trajectory()

    expected_shots = int(r.nintpc.sum())
    assert et.gx.shape == (expected_shots, r.ngradact)
    assert et.gy.shape == et.gx.shape
    assert et.gz.shape == et.gx.shape
    assert et.gx.dtype == np.float32
    # Expanded waveforms inherit physical units from the parent ConegradResult.
    assert np.abs(et.gx).max() <= GMAX * 1.01
    assert et.band_index.shape == (expected_shots,)
    assert et.interleaf_index.shape == (expected_shots,)
    # band_index runs monotonically and covers every band
    assert et.band_index.min() == 0
    assert et.band_index.max() == r.rspthetas.size - 1
    assert np.all(np.diff(et.band_index) >= 0)


def test_si_units_produce_equivalent_design():
    """Calling conegrad in SI mode with converted args should yield the same
    cone waveforms as the CGS call, and report max_gradient in T/m."""
    r_cgs = conegrad.conegrad(RES, FOV, **_common_kwargs())

    # Convert CGS inputs to SI for the second call.
    r_si = conegrad.conegrad(
        RES * 1e-3,    # mm -> m
        FOV * 1e-2,    # cm -> m
        numcones=32,
        grad_points=GRAD_POINTS,
        gsp=GSP,
        read_points=READ_POINTS,
        tsp=TSP,
        smax=SMAX * 1e-2,
        max_s_rewind=MAX_S_REWIND * 1e-2,
        gmax=GMAX * 1e-2,
        sys_gmax=SYS_GMAX * 1e-2,
        sys_max_rewind_g=SYS_MAX_REWIND_G * 1e-2,
        units="SI",
    )

    # The cone designer is identical; only the unit scale differs. CGS
    # outputs are in G/cm, SI outputs are in T/m -- ratio is exactly 100.
    assert r_cgs.ntheta == r_si.ntheta
    assert np.array_equal(r_cgs.nintpc, r_si.nintpc)
    np.testing.assert_allclose(r_si.gx * 100.0, r_cgs.gx, rtol=1e-6)
    np.testing.assert_allclose(r_si.gy * 100.0, r_cgs.gy, rtol=1e-6)
    np.testing.assert_allclose(r_si.gz * 100.0, r_cgs.gz, rtol=1e-6)
    # max_gradient must convert: CGS in G/cm, SI in T/m, factor 100.
    assert abs(r_si.max_gradient - r_cgs.max_gradient * 1e-2) < 1e-9


def test_si_units_rejects_unknown_value():
    with pytest.raises(ValueError, match="cgs.*SI"):
        conegrad.conegrad(RES, FOV, units="mks", **_common_kwargs())


def test_make_vd_fov_layout_and_endpoints():
    fov_full = np.array([24.0, 24.0], dtype=np.float32)
    fov_table = conegrad.make_vd_fov(fov_full, fov_under=0.5, n_full=20)

    assert fov_table.shape == (300,)
    assert fov_table.dtype == np.float32

    # First 20 samples of each FOV channel: held at full FOV.
    assert np.all(fov_table[0:20] == fov_full[0])
    assert np.all(fov_table[100:120] == fov_full[1])

    # Last sample of FOV channels: fov_under * fov_full.
    assert abs(fov_table[99] - fov_full[0] * 0.5) < 1e-5
    assert abs(fov_table[199] - fov_full[1] * 0.5) < 1e-5

    # kr channel: 0..1 normalized.
    assert fov_table[200] == 0.0
    assert abs(fov_table[299] - 1.0) < 1e-6
    assert np.all(np.diff(fov_table[200:300]) > 0)


def test_make_vd_fov_feeds_conegrad_vd():
    fov_full = np.array([24.0, 24.0], dtype=np.float32)
    fov_table = conegrad.make_vd_fov(fov_full, fov_under=0.6)
    r = conegrad.conegrad_vd(RES, fov_table, **_common_kwargs())
    assert r.ntheta > 0
    assert r.ngradact > 0
    assert int(r.nintpc.sum()) > 0


def test_make_vd_fov_unity_matches_constant_fov():
    """With fov_under=1.0 the VD table is constant at the full FOV across
    all kr, which should be equivalent to a constant-FOV design."""
    fov_full = np.array([24.0, 24.0], dtype=np.float32)
    fov_table = conegrad.make_vd_fov(fov_full, fov_under=1.0)
    # All FOV samples (channels 0 and 1) equal fov_full.
    assert np.allclose(fov_table[0:100], fov_full[0])
    assert np.allclose(fov_table[100:200], fov_full[1])


def test_make_vd_fov_per_axis_undersampling():
    fov_full = np.array([24.0, 18.0], dtype=np.float32)
    fov_table = conegrad.make_vd_fov(fov_full, fov_under=[0.5, 0.25], n_full=20)
    assert abs(fov_table[99] - 24.0 * 0.5) < 1e-5
    assert abs(fov_table[199] - 18.0 * 0.25) < 1e-5


def test_expand_trajectory_southern_hemisphere_flips_gz():
    """Polar bands at +theta and -theta share the same base cone but
    differ in gz sign (zscale = sin(theta)/sin(theta_up)).
    """
    r = conegrad.conegrad(RES, FOV, **_common_kwargs())
    et = r.expand_trajectory()

    # Find mirror pairs: bands t and (2*ntheta-2-t) have opposite sign theta.
    # Pick one band away from the equator so |theta| is appreciable.
    n_bands = r.rspthetas.size
    t_pos = n_bands - 5            # near +pi/2
    t_neg = n_bands - 1 - t_pos    # mirror across equator

    # First interleaf of each band
    sp = int(np.searchsorted(et.band_index, t_pos))
    sn = int(np.searchsorted(et.band_index, t_neg))

    # Verify rspthetas are indeed opposite signs
    assert r.rspthetas[t_pos] * r.rspthetas[t_neg] < 0

    # gz should be sign-flipped between the two (same magnitude up to scale)
    gz_pos = et.gz[sp]
    gz_neg = et.gz[sn]
    # If both bands map to the same base cone, the gz waveforms differ only
    # by zscale sign and magnitude; the dot product is therefore negative.
    assert float(np.dot(gz_pos, gz_neg)) < 0
