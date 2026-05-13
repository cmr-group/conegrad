# conegrad

<img src="assets/traj3.gif" width="200" />

Python wrapper around the `conegrad` and `conegrad_vd` C functions for designing
3D Cones MRI k-space trajectories.

Gurney, P.T., Hargreaves, B.A. and Nishimura, D.G. (2006), Design and analysis of a practical 3D cones trajectory. Magn. Reson. Med., 55: 575-582. https://doi.org/10.1002/mrm.20796

## Installation

```bash
pip install conegrad
```

## Quick start

```python
import numpy as np
import conegrad

res = [1e-3, 1e-3]  # [m], this is XY, and Z
fov = np.array([0.24, 0.24])  # [m], this is XY, and Z
gmax = 0.08  # [T/m]
smax = 80  # [T/m/s]
# Note: units can be in CGS, just change the `units` argument

dt = 4e-6  # [s]
t_readout = 3e-3 # [s]

rewind_flag = 1  # Do rewinding lobe

n_pts = int(t_readout / dt)
out = conegrad.conegrad(res, fov, 
                        gmax=gmax, smax=smax, 
                        gsp=dt, grad_points=n_pts, 
                        units='SI', rewind_flag=rewind_flag)

# out is now the base cones
print(result.gx.shape, result.gx.dtype)  # (numcones, ngradact) float32, G/cm
print(result.ntheta, result.nintpc)      # number of theta bands, interleaves per band

# To exand the base cones to all cones in the trajectory:
exp = out.expand_trajectory()

# Now exp contains the full gradient waveforms
print(exp.gx.shape)
```

## Arguments

The unit columns below are written as **cgs / SI** — both unit systems are
supported via the `units="cgs"` (default) or `units="SI"` keyword argument.

### Scan geometry (required positional)

| Arg | Type | Shape | Units (cgs / SI) | Meaning |
|---|---|---|---|---|
| `res` | float32 | `(2,)` | mm / m, × gamma_ratio | Spatial resolution `[xy, z]`. For proton, gamma_ratio = 1; for other nuclei, scale by `gamma_nucleus / gamma_proton`. |
| `fov` | float32 | `(2,)` for `conegrad`; `(300,)` for `conegrad_vd` | cm / m | Field of view. See "Variable-density FOV layout" below for `conegrad_vd`. |

### Sequence timing (required, keyword-only)

| Arg | Type | Units | Typical | Meaning |
|---|---|---|---|---|
| `gsp` | float | seconds | `4e-6` | Gradient sample period (i.e. 4 µs grid). |
| `grad_points` | int | samples | a few hundred to several thousand | Number of gradient samples during the readout window. |
| `tsp` | float, optional | seconds | `4e-6` | Readout (ADC) sample period. Defaults to `gsp` if omitted. |
| `read_points` | int, optional | samples | up to a few thousand | Number of oversampled readout sample points. Defaults to `grad_points` if omitted. |

**Note:** `grad_points` drives the cone design; `read_points`/`tsp` are only
used when writing the k-space trajectory file (`ktraj_out_flag=1`). Omit them
unless you need the file output.

### Scanner gradient limits (keyword-only)

| Arg | Type | Units (cgs / SI) | Typical 3T (cgs / SI) | Meaning |
|---|---|---|---|---|
| `gmax` | float | G/cm / T/m | 4.5 / 0.045 | Max gradient amplitude during readout. **Required.** |
| `smax` | float | G/cm/s / T/m/s | ~7960 / ~79.6 | Max slew rate during readout. **Required.** |
| `sys_max_rewind_g` | float, optional | G/cm / T/m | 2.6 / 0.026 | Max gradient during rewinder. Defaults to `gmax / sqrt(3)`. |
| `max_s_rewind` | float, optional | G/cm/s / T/m/s | ~4600 / ~46 | Max slew rate during rewinder. Defaults to `smax / sqrt(3)`. |
| `sys_gmax` | float, optional | G/cm / T/m | — | Hardware gradient ceiling, used only for internal int-amp scaling; the wrapper inverts it on output so it has no user-visible effect. Defaults to `gmax`. |

### Cone design (optional, keyword-only)

| Arg | Type | Default | Range | Meaning |
|---|---|---|---|---|
| `numcones` | int | 32 | 1..100 | Number of base cone waveforms designed (covers polar angles 0..π/2). |
| `precision` | float | 0.5 | 0..1 | Cones design tolerance. |
| `dcf` | float | 0.0 | 0..1 | Density compensation factor. |
| `oversample` | int | 1 | 1..1024 | Trajectory oversampling factor. |
| `mindens` | float | 0.0 | 0..10 | Minimum sampling density. |
| `rewind_flag` | int | 1 | 0 or 1 | `1` = full rewinder (return k to 0); anything else = ramp grad to 0 only. The original EPIC sequence also defined values 2/3/4 for alternative rewinder strategies, but those were handled outside `scones_design.c` and have no effect here. |
| `sym_flag` | int | 1 | 1, 2 | 1 = k-space symmetric (mirrors design across equator), 2 = asymmetric. |
| `acq_mode` | int | 1 | 1, 2 | 1 = Cones, 2 = Projection Reconstruction. |
| `rot_flag` | float | 1.0 | -1 or +1 | Direction of azimuthal rotation between interleaves. |
| `slab_kz` | float | 0.0 | — | kz phase across the slab (added to all kz values). |
| `t_fracx`, `t_fracy`, `t_fracz` | float | 0.0 | — | Per-axis fractional gradient delays (timing calibration). |

### File output (optional, keyword-only — all default to 0)

These trigger file writes into the current working directory. Defaults are off
since most users only want the in-memory arrays.

| Arg | Type | Effect |
|---|---|---|
| `output_grad` | int | 1 → write `kacq_cone_g{x,y,z}.<uid>` (per-cone gradient binaries). |
| `ktraj_flag` | int | 1 → run k-space-trajectory + density computation (needed for `nintpc`, `rspthetas`, `max_gradient`). **Defaults to 1** in this wrapper. |
| `ktraj_out_flag` | int | 1 → write `kacq_cone_{ktraj,dens,npc}.<uid>` AND populate `snr_efficiency`. |
| `endian_flag` | int | Output file endianness flag (currently inert — endian conversion code is commented in scones_design.c). |
| `rhkacq_uid` | int | UID used in output filenames. |
| `cones_plot_flag` | int | 0..3 — generate gnuplot scripts of gradient (1) / k-space (2) / both (3). |
| `verbose` | bool | True → forward C-side `printf` debug output to stdout. |
| `units` | str | `"cgs"` (default) or `"SI"`. Sets units for both inputs (m, T/m, T/m/s in SI mode; cm, G/cm, G/cm/s in cgs) and the returned `gx/gy/gz/max_gradient` (T/m vs G/cm). |

### Variable-density FOV layout (`conegrad_vd` only)

`fov` is a length-300 array laid out as three concatenated 100-sample channels:

| Slice | Meaning |
|---|---|
| `fov[0:100]` | FOV_xy as a function of normalized k-radius |
| `fov[100:200]` | FOV_z as a function of normalized k-radius |
| `fov[200:300]` | Normalized k-radius (typically `linspace(0, 1, 100)`) |

Use `conegrad.make_vd_fov(fov_full, fov_under=0.5, n_full=20)` to build this
table with the standard EPIC profile: inner `n_full` samples at the full FOV,
outer samples linearly ramping to `fov_under * fov_full`.

## Outputs

Both `conegrad()` and `conegrad_vd()` return a `ConegradResult`:

| Field | Type | Shape | Units | Meaning |
|---|---|---|---|---|
| `gx`, `gy`, `gz` | float32 | `(numcones, ngradact)` | G/cm (cgs) or T/m (SI) | Per-cone base waveforms in physical units. `gx[i]`, `gy[i]`, `gz[i]` are a **coupled triplet** for cone `i` — do not mix indices across axes. |
| `ntheta` | int | scalar | — | Number of polar bands designed in the upper hemisphere. Total bands in the full sphere acquisition = `2 * ntheta - 1`. |
| `nintpc` | int32 | `(2*ntheta-1,)` | shots | Number of azimuthal interleaves per polar band. `sum(nintpc)` is total shots in the full acquisition. |
| `rspthetas` | float32 | `(2*ntheta-1,)` | radians | Polar angle θ of each band, in `[-π/2, π/2]` (full hemisphere range). |
| `max_gradient` | float | scalar | G/cm | Max gradient amplitude actually used in the design. |
| `snr_efficiency` | float | scalar | percent | Only populated when `ktraj_out_flag=1` (which also writes files). 0 otherwise. |
| `ngradact` | int | scalar | samples | Per-cone waveform length (readout + rewinder). |

### Expanding to full per-shot trajectory

The 32 base waveforms are designed for the upper hemisphere only; the full
sphere is built at scan time by polar-scaling and azimuthally rotating each
base cone. Call `result.expand_trajectory(rot_flag=...)` to do this in Python:

```python
et = result.expand_trajectory(rot_flag=1.0)
et.gx.shape  # (sum(nintpc), ngradact)  — every shot in the acquisition
et.band_index, et.interleaf_index  # which polar band / which interleaf each shot is
```

Returns an `ExpandedTrajectory`:

| Field | Type | Shape | Meaning |
|---|---|---|---|
| `gx`, `gy`, `gz` | float32 | `(n_shots, ngradact)` | Per-shot gradient waveforms in the same physical units as the parent `ConegradResult` (G/cm in cgs mode, T/m in SI mode), already polar-scaled and z-rotated. |
| `band_index` | int32 | `(n_shots,)` | Polar band each shot belongs to (`0..2*ntheta-2`). |
| `interleaf_index` | int32 | `(n_shots,)` | Azimuthal interleaf within the band (`0..nintpc[band]-1`). |

The polar scaling and z-rotation applied per shot match the EPIC scan-time
code at `mm4dflow.e:17050+`. Specifically: for polar band `t` and interleaf
`k`, base cone `ci = ceil(|rspthetas[t]| / (π/2) · numcones) - 1`:

```
xyscale = cos(theta) / cos(theta_down)
zscale  = sin(theta) / sin(theta_up)         # automatically negative below the equator
phi     = k · rot_flag · 2π / nintpc[t]

gx_shot = xyscale · (cos(phi)·gx_base − sin(phi)·gy_base)
gy_shot = xyscale · (sin(phi)·gx_base + cos(phi)·gy_base)
gz_shot = zscale  · gz_base
```

## Building from source

```bash
pip install .
```

Requires a C and C++17 compiler, CMake >= 3.18, and Python >= 3.11. nanobind is
fetched automatically as a build dependency.
