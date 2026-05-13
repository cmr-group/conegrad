// nanobind bindings for the conegrad / conegrad_vd C functions.
//
// The underlying C code designs 3D Cones MRI k-space trajectories: for a given
// scanner / scan configuration it produces integer-amplitude gradient waveforms
// on three axes for many conical interleaves.
//
// This binding owns the output buffers, calls the C entry points, and returns
// the results as numpy arrays packaged inside a ConegradResult NamedTuple.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;

extern "C" {

// Verbose flag defined in scones_design.c. The Python wrapper toggles this
// to silence the ~167 printf debug statements scattered through the C code.
extern int conegrad_verbose;

int conegrad(int *nintpc, float *rspthetas, float *RES, float *FOV,
             int *gxcones_host, int *gycones_host, int *gzcones_host,
             int NUMCONES, int GRAD_POINTS, float GSP, int READ_POINTS,
             float TSP, float PRECISION, float DCF, int OVERSAMPLE,
             float SMAX, float maxSRewind, float GMAX, float sysGMAX,
             float SysMaxRewindG, int output_grad, int ktraj_flag,
             int ktraj_out_flag, int endian_flag, int rewind_flag,
             float MINDENS, float *maxgr2, float *snr_eff, int *ngradact,
             float t_fracx, float t_fracy, float t_fracz, int acq_mode,
             float rot_flag, int SymFlag, float SlabKz, int rhkacq_uid,
             int Cones_Plot_Flag, int *traj_length);

int conegrad_vd(int *nintpc, float *rspthetas, float *RES, float *FOV,
                int *gxcones_host, int *gycones_host, int *gzcones_host,
                int NUMCONES, int GRAD_POINTS, float GSP, int READ_POINTS,
                float TSP, float PRECISION, float DCF, int OVERSAMPLE,
                float SMAX, float maxSRewind, float GMAX, float sysGMAX,
                float SysMaxRewindG, int output_grad, int ktraj_flag,
                int ktraj_out_flag, int endian_flag, int rewind_flag,
                float MINDENS, float *maxgr2, float *snr_eff, int *ngradact,
                float t_fracx, float t_fracy, float t_fracz, int acq_mode,
                float rot_flag, int SymFlag, float SlabKz, int rhkacq_uid,
                int Cones_Plot_Flag, int *traj_length);

}  // extern "C"

// MAX_THETAS matches scones_design.h. The C functions write up to 2*ntheta-1
// entries into nintpc / rspthetas where ntheta <= MAX_THETAS.
constexpr int CG_MAX_THETAS = 1000;

// Safety margin added on top of GRAD_POINTS for the per-cone waveform buffer.
// ngradact = GRAD_POINTS + nramp, where nramp is bounded by the formula at
// scones_design.c:~3997. 5000 covers any realistic Gmax/Smax/GSP combination.
constexpr int CG_BUFFER_MARGIN = 5000;

using FloatArr = nb::ndarray<float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// Wrap a heap-allocated raw buffer as a numpy array. The capsule frees the
// buffer when Python garbage-collects the array.
template <typename T>
nb::ndarray<nb::numpy, T> make_owned_1d(T *data, size_t n) {
    nb::capsule owner(data, [](void *p) noexcept {
        delete[] static_cast<T *>(p);
    });
    size_t shape[1] = {n};
    return nb::ndarray<nb::numpy, T>(data, 1, shape, owner);
}

template <typename T>
nb::ndarray<nb::numpy, T> make_owned_2d(T *data, size_t rows, size_t cols) {
    nb::capsule owner(data, [](void *p) noexcept {
        delete[] static_cast<T *>(p);
    });
    size_t shape[2] = {rows, cols};
    return nb::ndarray<nb::numpy, T>(data, 2, shape, owner);
}

// Validate the input arrays and run the chosen C entry point. Shared between
// conegrad() and conegrad_vd() — the only difference is the expected FOV shape
// and which C function is called.
struct CallResult {
    nb::object gx;
    nb::object gy;
    nb::object gz;
    int ntheta;
    nb::object nintpc;
    nb::object rspthetas;
    float max_gradient;
    float snr_efficiency;
    int ngradact;
    nb::object traj_length;
};

using CFunc = int (*)(int *, float *, float *, float *, int *, int *, int *,
                      int, int, float, int, float, float, float, int, float,
                      float, float, float, float, int, int, int, int, int,
                      float, float *, float *, int *, float, float, float,
                      int, float, int, float, int, int, int *);

static CallResult call_conegrad_impl(
    CFunc cfunc,
    const char *cfunc_name,
    FloatArr res,
    FloatArr fov,
    size_t expected_fov_size,
    int numcones, int grad_points, float gsp,
    int read_points, float tsp,
    float precision, float dcf, int oversample,
    float smax, float max_s_rewind,
    float gmax, float sys_gmax, float sys_max_rewind_g,
    int output_grad, int ktraj_flag, int ktraj_out_flag,
    int endian_flag, int rewind_flag, float mindens,
    float t_fracx, float t_fracy, float t_fracz,
    int acq_mode, float rot_flag, int sym_flag,
    float slab_kz, int rhkacq_uid, int cones_plot_flag,
    bool verbose)
{
    if (res.shape(0) != 2) {
        throw std::invalid_argument("res must be a 1-D float32 array of length 2 (xy, z)");
    }
    if (fov.shape(0) != expected_fov_size) {
        throw std::invalid_argument(std::string(cfunc_name) +
            ": fov must be a 1-D float32 array of length " +
            std::to_string(expected_fov_size));
    }
    if (numcones < 1 || numcones > 100) {
        throw std::invalid_argument("numcones must be in [1, 100]");
    }
    if (grad_points < 1) {
        throw std::invalid_argument("grad_points must be >= 1");
    }

    conegrad_verbose = verbose ? 1 : 0;

    int max_pts_per_cone = grad_points + CG_BUFFER_MARGIN;
    int total_grad_size = numcones * max_pts_per_cone;
    int max_theta_pts = 2 * CG_MAX_THETAS;

    // Use zeroed allocations — easier to spot bugs if a region is left untouched.
    int32_t *gx = new int32_t[total_grad_size]();
    int32_t *gy = new int32_t[total_grad_size]();
    int32_t *gz = new int32_t[total_grad_size]();
    int32_t *nintpc_buf = new int32_t[max_theta_pts]();
    float   *rspthetas_buf = new float[max_theta_pts]();
    int32_t *traj_length_buf = new int32_t[numcones]();

    float maxgr = 0.0f;
    float snr_eff = 0.0f;
    int ngradact = 0;
    int ntheta = 0;

    try {
        ntheta = cfunc(
            nintpc_buf, rspthetas_buf,
            res.data(), fov.data(),
            gx, gy, gz,
            numcones, grad_points, gsp, read_points,
            tsp, precision, dcf, oversample,
            smax, max_s_rewind, gmax, sys_gmax, sys_max_rewind_g,
            output_grad, ktraj_flag, ktraj_out_flag, endian_flag, rewind_flag,
            mindens, &maxgr, &snr_eff, &ngradact,
            t_fracx, t_fracy, t_fracz, acq_mode, rot_flag, sym_flag,
            slab_kz, rhkacq_uid, cones_plot_flag, traj_length_buf);
    } catch (...) {
        delete[] gx; delete[] gy; delete[] gz;
        delete[] nintpc_buf; delete[] rspthetas_buf;
        delete[] traj_length_buf;
        throw;
    }

    if (ngradact > max_pts_per_cone) {
        delete[] gx; delete[] gy; delete[] gz;
        delete[] nintpc_buf; delete[] rspthetas_buf;
        delete[] traj_length_buf;
        throw std::runtime_error(
            std::string(cfunc_name) + ": ngradact (" + std::to_string(ngradact) +
            ") exceeded allocated buffer (" + std::to_string(max_pts_per_cone) +
            "). Increase grad_points or CG_BUFFER_MARGIN.");
    }
    if (ntheta < 0 || ntheta > CG_MAX_THETAS) {
        delete[] gx; delete[] gy; delete[] gz;
        delete[] nintpc_buf; delete[] rspthetas_buf;
        delete[] traj_length_buf;
        throw std::runtime_error(
            std::string(cfunc_name) + ": returned ntheta=" + std::to_string(ntheta) +
            " out of range [0, " + std::to_string(CG_MAX_THETAS) + "]");
    }

    // The C code writes gxcones_host[j + i * (*ngradact)], so the per-cone
    // stride is ngradact (the output, not the allocated buffer width). The
    // data is already densely packed in the first numcones * ngradact slots;
    // we simply present that slice as a (numcones, ngradact) ndarray below.
    // Anything past numcones * ngradact is untouched zero-init padding.

    size_t theta_n = (ntheta > 0) ? static_cast<size_t>(2 * ntheta - 1) : 0;

    CallResult r;
    r.gx = nb::cast(make_owned_2d<int32_t>(gx, numcones, ngradact));
    r.gy = nb::cast(make_owned_2d<int32_t>(gy, numcones, ngradact));
    r.gz = nb::cast(make_owned_2d<int32_t>(gz, numcones, ngradact));
    r.ntheta = ntheta;
    // For nintpc / rspthetas we hand back arrays sized 2*ntheta-1 but the
    // underlying buffer is still max_theta_pts; the capsule frees the full
    // allocation regardless.
    r.nintpc = nb::cast(make_owned_1d<int32_t>(nintpc_buf, theta_n));
    r.rspthetas = nb::cast(make_owned_1d<float>(rspthetas_buf, theta_n));
    r.max_gradient = maxgr;
    r.snr_efficiency = snr_eff;
    r.ngradact = ngradact;
    r.traj_length = nb::cast(make_owned_1d<int32_t>(
        traj_length_buf, static_cast<size_t>(numcones)));
    return r;
}

NB_MODULE(_conegrad, m) {
    m.doc() = "Low-level nanobind bindings for the conegrad / conegrad_vd C functions.";

    // The Python wrapper builds the ConegradResult NamedTuple from these
    // primitive returns; the binding itself returns a plain tuple to keep
    // nanobind dependencies minimal.
    auto pack = [](const CallResult &r) {
        return nb::make_tuple(r.gx, r.gy, r.gz,
                              r.ntheta, r.nintpc, r.rspthetas,
                              r.max_gradient, r.snr_efficiency, r.ngradact,
                              r.traj_length);
    };

    m.def("conegrad",
          [pack](FloatArr res, FloatArr fov,
                 int numcones, int grad_points, float gsp,
                 int read_points, float tsp,
                 float precision, float dcf, int oversample,
                 float smax, float max_s_rewind,
                 float gmax, float sys_gmax, float sys_max_rewind_g,
                 int output_grad, int ktraj_flag, int ktraj_out_flag,
                 int endian_flag, int rewind_flag, float mindens,
                 float t_fracx, float t_fracy, float t_fracz,
                 int acq_mode, float rot_flag, int sym_flag,
                 float slab_kz, int rhkacq_uid, int cones_plot_flag,
                 bool verbose) {
            auto r = call_conegrad_impl(
                &conegrad, "conegrad",
                res, fov, /*expected_fov_size=*/2,
                numcones, grad_points, gsp, read_points, tsp,
                precision, dcf, oversample,
                smax, max_s_rewind, gmax, sys_gmax, sys_max_rewind_g,
                output_grad, ktraj_flag, ktraj_out_flag, endian_flag, rewind_flag,
                mindens, t_fracx, t_fracy, t_fracz,
                acq_mode, rot_flag, sym_flag, slab_kz, rhkacq_uid, cones_plot_flag,
                verbose);
            return pack(r);
          },
          nb::arg("res"), nb::arg("fov"),
          nb::arg("numcones"), nb::arg("grad_points"), nb::arg("gsp"),
          nb::arg("read_points"), nb::arg("tsp"),
          nb::arg("precision"), nb::arg("dcf"), nb::arg("oversample"),
          nb::arg("smax"), nb::arg("max_s_rewind"),
          nb::arg("gmax"), nb::arg("sys_gmax"), nb::arg("sys_max_rewind_g"),
          nb::arg("output_grad"), nb::arg("ktraj_flag"), nb::arg("ktraj_out_flag"),
          nb::arg("endian_flag"), nb::arg("rewind_flag"), nb::arg("mindens"),
          nb::arg("t_fracx"), nb::arg("t_fracy"), nb::arg("t_fracz"),
          nb::arg("acq_mode"), nb::arg("rot_flag"), nb::arg("sym_flag"),
          nb::arg("slab_kz"), nb::arg("rhkacq_uid"), nb::arg("cones_plot_flag"),
          nb::arg("verbose"));

    // conegrad_vd takes a variable-density FOV table laid out as three blocks of
    // length Nkr=100: [FOVxy(kr), FOVz(kr), kr_normalized]. Expected length 300.
    m.def("conegrad_vd",
          [pack](FloatArr res, FloatArr fov,
                 int numcones, int grad_points, float gsp,
                 int read_points, float tsp,
                 float precision, float dcf, int oversample,
                 float smax, float max_s_rewind,
                 float gmax, float sys_gmax, float sys_max_rewind_g,
                 int output_grad, int ktraj_flag, int ktraj_out_flag,
                 int endian_flag, int rewind_flag, float mindens,
                 float t_fracx, float t_fracy, float t_fracz,
                 int acq_mode, float rot_flag, int sym_flag,
                 float slab_kz, int rhkacq_uid, int cones_plot_flag,
                 bool verbose) {
            auto r = call_conegrad_impl(
                &conegrad_vd, "conegrad_vd",
                res, fov, /*expected_fov_size=*/300,
                numcones, grad_points, gsp, read_points, tsp,
                precision, dcf, oversample,
                smax, max_s_rewind, gmax, sys_gmax, sys_max_rewind_g,
                output_grad, ktraj_flag, ktraj_out_flag, endian_flag, rewind_flag,
                mindens, t_fracx, t_fracy, t_fracz,
                acq_mode, rot_flag, sym_flag, slab_kz, rhkacq_uid, cones_plot_flag,
                verbose);
            return pack(r);
          },
          nb::arg("res"), nb::arg("fov"),
          nb::arg("numcones"), nb::arg("grad_points"), nb::arg("gsp"),
          nb::arg("read_points"), nb::arg("tsp"),
          nb::arg("precision"), nb::arg("dcf"), nb::arg("oversample"),
          nb::arg("smax"), nb::arg("max_s_rewind"),
          nb::arg("gmax"), nb::arg("sys_gmax"), nb::arg("sys_max_rewind_g"),
          nb::arg("output_grad"), nb::arg("ktraj_flag"), nb::arg("ktraj_out_flag"),
          nb::arg("endian_flag"), nb::arg("rewind_flag"), nb::arg("mindens"),
          nb::arg("t_fracx"), nb::arg("t_fracy"), nb::arg("t_fracz"),
          nb::arg("acq_mode"), nb::arg("rot_flag"), nb::arg("sym_flag"),
          nb::arg("slab_kz"), nb::arg("rhkacq_uid"), nb::arg("cones_plot_flag"),
          nb::arg("verbose"));
}
