// Generated by dust (version 0.10.0) - do not edit
#include <cpp11.hpp>
[[cpp11::register]]
SEXP dust_model_alloc(cpp11::list r_pars, bool pars_multi, size_t step,
                         cpp11::sexp r_n_particles, size_t n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp device_config);

[[cpp11::register]]
SEXP dust_model_run(SEXP ptr, size_t step_end, bool device);

[[cpp11::register]]
SEXP dust_model_simulate(SEXP ptr, cpp11::sexp step_end, bool device);

[[cpp11::register]]
SEXP dust_model_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                SEXP r_step, SEXP r_set_initial_state);

[[cpp11::register]]
SEXP dust_model_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
size_t dust_model_step(SEXP ptr);

[[cpp11::register]]
void dust_model_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_model_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_model_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_model_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_model_set_data(SEXP ptr, cpp11::list data);

[[cpp11::register]]
SEXP dust_model_compare_data(SEXP ptr, bool device);

[[cpp11::register]]
SEXP dust_model_filter(SEXP ptr, bool save_trajectories,
                          cpp11::sexp step_snapshot,
                          bool device);

[[cpp11::register]]
cpp11::sexp dust_model_capabilities();

[[cpp11::register]]
void dust_model_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_model_n_state(SEXP ptr);

[[cpp11::register]]
cpp11::sexp dust_model_device_info();

#include <dust/dust.hpp>
#include <dust/interface/dust.hpp>

// Generated by odin.dust (version 0.2.11) - do not edit
template <typename real_type, typename T, typename U>
HOSTDEVICE real_type fmodr(T x, U y) {
  real_type tmp = std::fmod(static_cast<real_type>(x),
                            static_cast<real_type>(y));
  if (tmp * y < 0) {
    tmp += y;
  }
  return tmp;
}

// These exist to support the model on the gpu, as in C++14 std::min
// and std::max are constexpr and error without --expt-relaxed-constexpr
template <typename T>
HOSTDEVICE T odin_min(T x, T y) {
  return x < y ? x : y;
}

template <typename T>
HOSTDEVICE T odin_max(T x, T y) {
  return x > y ? x : y;
}
// [[dust::class(model)]]
// [[dust::param(alpha, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_A, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_E, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_F, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_I, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_R, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(delta_S, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(omega, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_F, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_I, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_R, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_S, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(sigma, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_s, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(A0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(F0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S10, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S20, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(U0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t0, has_default = TRUE, default_value = 0L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class model {
public:
  typedef double real_type;
  typedef dust::random::generator<real_type> rng_state_type;
  typedef dust::no_data data_type;
  struct shared_type {
    real_type A0;
    real_type E0;
    real_type F0;
    real_type I0;
    real_type R0;
    real_type S10;
    real_type S20;
    real_type U0;
    real_type alpha;
    real_type beta;
    real_type delta_A;
    real_type delta_E;
    real_type delta_F;
    real_type delta_I;
    real_type delta_R;
    real_type delta_S;
    real_type dt;
    real_type initial_A;
    real_type initial_E;
    real_type initial_F;
    real_type initial_I;
    real_type initial_N;
    real_type initial_R;
    real_type initial_S1;
    real_type initial_S2;
    real_type initial_U;
    real_type initial_igas_inc;
    real_type initial_pharyngitis_inc;
    real_type initial_scarlet_fever_inc;
    real_type initial_time;
    real_type omega;
    real_type p_F;
    real_type p_I;
    real_type p_R;
    real_type p_S;
    real_type pi;
    real_type r_A;
    real_type r_AR;
    real_type r_AU;
    real_type r_E;
    real_type r_EI;
    real_type r_ES;
    real_type r_F;
    real_type r_FR;
    real_type r_I;
    real_type r_IR;
    real_type r_R;
    real_type r_RU;
    real_type r_S1;
    real_type r_S2;
    real_type r_SF;
    real_type r_SR;
    real_type r_SS;
    real_type sigma;
    real_type steps_per_week;
    real_type t0;
    real_type t_s;
  };
  struct internal_type {
  };
  model(const dust::pars_type<model>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() {
    return 13;
  }
  std::vector<real_type> initial(size_t step) {
    std::vector<real_type> state(13);
    state[0] = shared->initial_time;
    state[1] = shared->initial_U;
    state[2] = shared->initial_A;
    state[3] = shared->initial_E;
    state[4] = shared->initial_I;
    state[5] = shared->initial_S1;
    state[6] = shared->initial_S2;
    state[7] = shared->initial_F;
    state[8] = shared->initial_R;
    state[9] = shared->initial_N;
    state[10] = shared->initial_pharyngitis_inc;
    state[11] = shared->initial_scarlet_fever_inc;
    state[12] = shared->initial_igas_inc;
    return state;
  }
  HOST void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type time = state[0];
    const real_type U = state[1];
    const real_type E = state[3];
    const real_type A = state[2];
    const real_type S1 = state[5];
    const real_type S2 = state[6];
    const real_type F = state[7];
    const real_type I = state[4];
    const real_type R = state[8];
    const real_type N = state[9];
    const real_type pharyngitis_inc = state[10];
    const real_type scarlet_fever_inc = state[11];
    const real_type igas_inc = state[12];
    real_type seasonality = 1 + shared->sigma * std::cos(2 * shared->pi * (shared->t0 + time - shared->t_s) / (real_type) 365.25);
    real_type lambda = shared->beta * seasonality * (A + S1 + S2) / (real_type) N;
    real_type n_xU = dust::random::poisson<real_type>(rng_state, shared->alpha * shared->dt);
    state_next[0] = (step + 1) * shared->dt;
    real_type n_A = dust::random::binomial<real_type>(rng_state, A, 1 - std::exp(- shared->r_A * shared->dt));
    real_type n_E = dust::random::binomial<real_type>(rng_state, E, 1 - std::exp(- shared->r_E * shared->dt));
    real_type n_F = dust::random::binomial<real_type>(rng_state, F, 1 - std::exp(- shared->r_F * shared->dt));
    real_type n_I = dust::random::binomial<real_type>(rng_state, I, 1 - std::exp(- shared->r_I * shared->dt));
    real_type n_R = dust::random::binomial<real_type>(rng_state, R, 1 - std::exp(- shared->r_R * shared->dt));
    real_type n_S1 = dust::random::binomial<real_type>(rng_state, S1, 1 - std::exp(- shared->r_S1 * shared->dt));
    real_type n_S2 = dust::random::binomial<real_type>(rng_state, S2, 1 - std::exp(- shared->r_S2 * shared->dt));
    real_type r_UA = (1 - shared->p_S) * lambda;
    real_type r_UE = shared->p_S * lambda;
    real_type n_Ax = dust::random::binomial<real_type>(rng_state, n_A, shared->omega / (real_type) shared->r_A);
    real_type n_Ex = dust::random::binomial<real_type>(rng_state, n_E, shared->omega / (real_type) shared->r_E);
    real_type n_Fx = dust::random::binomial<real_type>(rng_state, n_F, shared->omega / (real_type) shared->r_F);
    real_type n_Ix = dust::random::binomial<real_type>(rng_state, n_I, shared->omega / (real_type) shared->r_I);
    real_type n_Rx = dust::random::binomial<real_type>(rng_state, n_R, shared->omega / (real_type) shared->r_R);
    real_type n_S1x = dust::random::binomial<real_type>(rng_state, n_S1, shared->omega / (real_type) shared->r_S1);
    real_type n_S2x = dust::random::binomial<real_type>(rng_state, n_S2, shared->omega / (real_type) shared->r_S2);
    real_type r_U = r_UE + r_UA + shared->omega;
    real_type n_AR = dust::random::binomial<real_type>(rng_state, n_A - n_Ax, shared->p_R);
    real_type n_EI = dust::random::binomial<real_type>(rng_state, n_E - n_Ex, shared->p_I);
    real_type n_FR = n_F - n_Fx;
    real_type n_IR = n_I - n_Ix;
    real_type n_RU = n_R - n_Rx;
    real_type n_SF = dust::random::binomial<real_type>(rng_state, n_S1 - n_S1x, shared->p_F);
    real_type n_SR = n_S2 - n_S2x;
    real_type n_U = dust::random::binomial<real_type>(rng_state, U, 1 - std::exp(- r_U * shared->dt));
    real_type n_AU = n_A - n_Ax - n_AR;
    real_type n_ES = n_E - n_Ex - n_EI;
    real_type n_SS = n_S1 - n_S1x - n_SF;
    real_type n_Ux = dust::random::binomial<real_type>(rng_state, n_U, shared->omega / (real_type) r_U);
    state_next[7] = F + n_SF - n_FR - n_Fx;
    state_next[4] = I + n_EI - n_IR - n_Ix;
    state_next[8] = R + n_AR + n_SR + n_FR + n_IR - n_RU - n_Rx;
    state_next[12] = ((fmodr<real_type>(step, shared->steps_per_week) == 0 ? n_EI : igas_inc + n_EI));
    state_next[11] = ((fmodr<real_type>(step, shared->steps_per_week) == 0 ? n_SF : scarlet_fever_inc + n_SF));
    real_type n_Nx = n_Ux + n_Ex + n_Ax + n_S1x + n_S2x + n_Fx + n_Ix + n_Rx;
    real_type n_UE = dust::random::binomial<real_type>(rng_state, n_U - n_Ux, shared->p_S);
    state_next[5] = S1 + n_ES - n_SF - n_SS - n_S1x;
    state_next[6] = S2 + n_SS - n_SR - n_S2x;
    state_next[10] = ((fmodr<real_type>(step, shared->steps_per_week) == 0 ? n_SS + n_SF : pharyngitis_inc + n_SS + n_SF));
    real_type n_UA = n_U - n_Ux - n_UE;
    state_next[3] = E + n_UE - n_ES - n_EI - n_Ex;
    state_next[9] = N + n_xU - n_Nx;
    state_next[2] = A + n_UA - n_AU - n_AR - n_Ax;
    state_next[1] = U + n_xU - n_UE - n_UA + n_AU + n_RU - n_Ux;
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
HOSTDEVICE real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace dust {
template<>
dust::pars_type<model> dust_pars<model>(cpp11::list user) {
  typedef typename model::real_type real_type;
  auto shared = std::make_shared<model::shared_type>();
  model::internal_type internal;
  shared->initial_igas_inc = 0;
  shared->initial_pharyngitis_inc = 0;
  shared->initial_scarlet_fever_inc = 0;
  shared->initial_time = 0;
  shared->pi = 3.14159265358979;
  shared->steps_per_week = 7;
  shared->dt = 1 / (real_type) shared->steps_per_week;
  shared->alpha = NA_REAL;
  shared->beta = NA_REAL;
  shared->delta_A = NA_REAL;
  shared->delta_E = NA_REAL;
  shared->delta_F = NA_REAL;
  shared->delta_I = NA_REAL;
  shared->delta_R = NA_REAL;
  shared->delta_S = NA_REAL;
  shared->omega = NA_REAL;
  shared->p_F = NA_REAL;
  shared->p_I = NA_REAL;
  shared->p_R = NA_REAL;
  shared->p_S = NA_REAL;
  shared->sigma = NA_REAL;
  shared->t_s = NA_REAL;
  shared->A0 = 0;
  shared->E0 = 0;
  shared->F0 = 0;
  shared->I0 = 0;
  shared->R0 = 0;
  shared->S10 = 0;
  shared->S20 = 0;
  shared->U0 = 0;
  shared->t0 = 0;
  shared->A0 = user_get_scalar<real_type>(user, "A0", shared->A0, NA_REAL, NA_REAL);
  shared->E0 = user_get_scalar<real_type>(user, "E0", shared->E0, NA_REAL, NA_REAL);
  shared->F0 = user_get_scalar<real_type>(user, "F0", shared->F0, NA_REAL, NA_REAL);
  shared->I0 = user_get_scalar<real_type>(user, "I0", shared->I0, NA_REAL, NA_REAL);
  shared->R0 = user_get_scalar<real_type>(user, "R0", shared->R0, NA_REAL, NA_REAL);
  shared->S10 = user_get_scalar<real_type>(user, "S10", shared->S10, NA_REAL, NA_REAL);
  shared->S20 = user_get_scalar<real_type>(user, "S20", shared->S20, NA_REAL, NA_REAL);
  shared->U0 = user_get_scalar<real_type>(user, "U0", shared->U0, NA_REAL, NA_REAL);
  shared->alpha = user_get_scalar<real_type>(user, "alpha", shared->alpha, NA_REAL, NA_REAL);
  shared->beta = user_get_scalar<real_type>(user, "beta", shared->beta, NA_REAL, NA_REAL);
  shared->delta_A = user_get_scalar<real_type>(user, "delta_A", shared->delta_A, NA_REAL, NA_REAL);
  shared->delta_E = user_get_scalar<real_type>(user, "delta_E", shared->delta_E, NA_REAL, NA_REAL);
  shared->delta_F = user_get_scalar<real_type>(user, "delta_F", shared->delta_F, NA_REAL, NA_REAL);
  shared->delta_I = user_get_scalar<real_type>(user, "delta_I", shared->delta_I, NA_REAL, NA_REAL);
  shared->delta_R = user_get_scalar<real_type>(user, "delta_R", shared->delta_R, NA_REAL, NA_REAL);
  shared->delta_S = user_get_scalar<real_type>(user, "delta_S", shared->delta_S, NA_REAL, NA_REAL);
  shared->omega = user_get_scalar<real_type>(user, "omega", shared->omega, NA_REAL, NA_REAL);
  shared->p_F = user_get_scalar<real_type>(user, "p_F", shared->p_F, NA_REAL, NA_REAL);
  shared->p_I = user_get_scalar<real_type>(user, "p_I", shared->p_I, NA_REAL, NA_REAL);
  shared->p_R = user_get_scalar<real_type>(user, "p_R", shared->p_R, NA_REAL, NA_REAL);
  shared->p_S = user_get_scalar<real_type>(user, "p_S", shared->p_S, NA_REAL, NA_REAL);
  shared->sigma = user_get_scalar<real_type>(user, "sigma", shared->sigma, NA_REAL, NA_REAL);
  shared->t0 = user_get_scalar<real_type>(user, "t0", shared->t0, NA_REAL, NA_REAL);
  shared->t_s = user_get_scalar<real_type>(user, "t_s", shared->t_s, NA_REAL, NA_REAL);
  shared->initial_A = shared->A0;
  shared->initial_E = shared->E0;
  shared->initial_F = shared->F0;
  shared->initial_I = shared->I0;
  shared->initial_N = shared->U0 + shared->A0 + shared->E0 + shared->I0 + shared->S10 + shared->S20 + shared->F0 + shared->R0;
  shared->initial_R = shared->R0;
  shared->initial_S1 = shared->S10;
  shared->initial_S2 = shared->S20;
  shared->initial_U = shared->U0;
  shared->r_AR = shared->p_R / (real_type) shared->delta_A;
  shared->r_AU = (1 - shared->p_R) / (real_type) shared->delta_A;
  shared->r_EI = shared->p_I / (real_type) shared->delta_E;
  shared->r_ES = (1 - shared->p_I) / (real_type) shared->delta_E;
  shared->r_FR = 1 / (real_type) shared->delta_F;
  shared->r_IR = 1 / (real_type) shared->delta_I;
  shared->r_RU = 1 / (real_type) shared->delta_R;
  shared->r_SF = shared->p_F / (real_type) shared->delta_S;
  shared->r_SR = 1 / (real_type) shared->delta_S;
  shared->r_SS = (1 - shared->p_F) / (real_type) shared->delta_S;
  shared->r_A = shared->r_AU + shared->r_AR + shared->omega;
  shared->r_E = shared->r_EI + shared->r_ES + shared->omega;
  shared->r_F = shared->r_FR + shared->omega;
  shared->r_I = shared->r_IR + shared->omega;
  shared->r_R = shared->r_RU + shared->omega;
  shared->r_S1 = shared->r_SF + shared->r_SS + shared->omega;
  shared->r_S2 = shared->r_SR + shared->omega;
  return dust::pars_type<model>(shared, internal);
}
template <>
cpp11::sexp dust_info<model>(const dust::pars_type<model>& pars) {
  const model::internal_type internal = pars.internal;
  const std::shared_ptr<const model::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "U", "A", "E", "I", "S1", "S2", "F", "R", "N", "pharyngitis_inc", "scarlet_fever_inc", "igas_inc"});
  cpp11::writable::list dim(13);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({1});
  dim[6] = cpp11::writable::integers({1});
  dim[7] = cpp11::writable::integers({1});
  dim[8] = cpp11::writable::integers({1});
  dim[9] = cpp11::writable::integers({1});
  dim[10] = cpp11::writable::integers({1});
  dim[11] = cpp11::writable::integers({1});
  dim[12] = cpp11::writable::integers({1});
  dim.names() = nms;
  cpp11::writable::list index(13);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = cpp11::writable::integers({6});
  index[6] = cpp11::writable::integers({7});
  index[7] = cpp11::writable::integers({8});
  index[8] = cpp11::writable::integers({9});
  index[9] = cpp11::writable::integers({10});
  index[10] = cpp11::writable::integers({11});
  index[11] = cpp11::writable::integers({12});
  index[12] = cpp11::writable::integers({13});
  index.names() = nms;
  size_t len = 13;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

SEXP dust_model_alloc(cpp11::list r_pars, bool pars_multi, size_t step,
                         cpp11::sexp r_n_particles, size_t n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp device_config) {
  return dust::r::dust_alloc<model>(r_pars, pars_multi, step, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        device_config);
}

SEXP dust_model_run(SEXP ptr, size_t step_end, bool device) {
  return dust::r::dust_run<model>(ptr, step_end, device);
}

SEXP dust_model_simulate(SEXP ptr, cpp11::sexp step_end, bool device) {
  return dust::r::dust_simulate<model>(ptr, step_end, device);
}

SEXP dust_model_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                SEXP r_step, SEXP r_set_initial_state) {
  return dust::r::dust_update_state<model>(ptr, r_pars, r_state, r_step,
                                               r_set_initial_state);
}

SEXP dust_model_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model>(ptr, r_index);
}

size_t dust_model_step(SEXP ptr) {
  return dust::r::dust_step<model>(ptr);
}

void dust_model_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model>(ptr, r_index);
}

SEXP dust_model_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model>(ptr, r_weights);
}

SEXP dust_model_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model>(ptr, first_only, last_only);
}

SEXP dust_model_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_model_set_data(SEXP ptr, cpp11::list data) {
  dust::r::dust_set_data<model>(ptr, data);
  return R_NilValue;
}

SEXP dust_model_compare_data(SEXP ptr, bool device) {
  return dust::r::dust_compare_data<model>(ptr, device);
}

SEXP dust_model_filter(SEXP ptr, bool save_trajectories,
                          cpp11::sexp step_snapshot,
                          bool device) {
  return dust::r::dust_filter<model>(ptr, save_trajectories, step_snapshot, device);
}

cpp11::sexp dust_model_capabilities() {
  return dust::r::dust_capabilities<model>();
}

void dust_model_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model>(ptr, n_threads);
}

int dust_model_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model>(ptr);
}

cpp11::sexp dust_model_device_info() {
  return dust::cuda::device_info<model::real_type>();
}
