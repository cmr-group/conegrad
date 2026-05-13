"""Quick smoke test: call conegrad with reasonable MRI defaults and print
the shape/summary of the result. Mirrors the parameter regime from the
mm4dflow.e usage so the C designer has a feasible problem to solve.
"""
import numpy as np

import conegrad


def main():
    # 1 mm isotropic, 24 cm FOV (gamma_ratio = 1.0 for proton imaging)
    res = np.array([1.0, 1.0], dtype=np.float32)
    fov = np.array([24.0, 24.0], dtype=np.float32)

    # Typical 3T GE scanner: 5 G/cm peak, 333 us rise time to peak.
    sys_gmax = 5.0
    gmax = 0.9 * sys_gmax           # GmaxScale = 0.9
    # SlewScale=0.75, 1e6*xfs/xrt = 1e6*5/333 = 15015 G/cm/s, then /sqrt(2)/(sqrt(3))
    smax = 0.75 * 15015 / np.sqrt(2)
    max_s_rewind = 0.75 * 15015 / np.sqrt(3)
    sys_max_rewind_g = 0.9 / np.sqrt(3) * sys_gmax

    # Sample timing
    gsp = 4e-6
    tsp = 4e-6
    read_points = 500
    grad_points = int(np.ceil(read_points * (2 * tsp / gsp)))

    print(f"Inputs: res={res}, fov={fov}, gmax={gmax:.2f}, smax={smax:.0f}, "
          f"grad_points={grad_points}")

    r = conegrad.conegrad(
        res, fov,
        numcones=32,
        grad_points=grad_points,
        gsp=gsp,
        read_points=read_points,
        tsp=tsp,
        smax=smax, max_s_rewind=max_s_rewind,
        gmax=gmax, sys_gmax=sys_gmax, sys_max_rewind_g=sys_max_rewind_g,
        verbose=False,
    )

    print("\n=== Result ===")
    print(f"gx shape:           {r.gx.shape}, dtype={r.gx.dtype}")
    print(f"gy shape:           {r.gy.shape}, dtype={r.gy.dtype}")
    print(f"gz shape:           {r.gz.shape}, dtype={r.gz.dtype}")
    print(f"ntheta:             {r.ntheta}")
    print(f"nintpc shape:       {r.nintpc.shape}, total interleaves={int(r.nintpc.sum())}")
    print(f"rspthetas range:    [{r.rspthetas.min():.3f}, {r.rspthetas.max():.3f}] rad")
    print(f"max_gradient:       {r.max_gradient:.4f}")
    print(f"snr_efficiency:     {r.snr_efficiency:.4f}")
    print(f"ngradact:           {r.ngradact}")

    # Sanity checks: arrays should not be all zeros
    assert r.gx.any(), "gx is all zero"
    assert r.gy.any(), "gy is all zero"
    assert r.ntheta > 0, f"ntheta should be positive, got {r.ntheta}"
    print("\nSmoke test passed.")


if __name__ == "__main__":
    main()
