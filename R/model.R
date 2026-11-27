##' @name model
##' @title Stochastic model of GAS
##' @description A `dust2` system generator for the stochastic compartmental
##' model of group A Streptococcus transmission, compiled from the odin code
##' in `inst/odin/model.R`.
##'
##' The model runs in discrete time. Time and durations are measured in days,
##' and the default time step is one day; transmission and demographic rates
##' are given per week. Flows such as `igas_inc` accumulate over each week,
##' and are reset every seven days. The model includes a comparison to weekly
##' surveillance data, so it can be used with [dust2::dust_filter_create()].
##'
##' Use [model_parameters()] to create a suitable list of parameters.
##' @export model
##' @examples
##' pars <- example_parameters()
##' sys <- dust2::dust_system_create(model, pars, n_particles = 5)
##' dust2::dust_system_set_state_initial(sys)
##' y <- dust2::dust_system_simulate(sys, c(0, 7, 14))
##' dust2::dust_unpack_state(sys, y)$igas_inc
NULL

##' @title Simulate the model
##' @description Simulate the model forwards from day zero, returning the
##' state at each of a set of times.
##' @param pars A list of model parameters, as created by
##' [model_parameters()]
##' @param times A vector of increasing times (in days) at which to return
##' the state
##' @param n_particles The number of particles to simulate
##' @param seed The random number seed, see [dust2::dust_system_create()]
##' @param deterministic Logical, indicating if the model should be run
##' deterministically, replacing random draws by their expectations
##' @param n_threads The number of threads to use
##' @return An array with dimensions state x particle x time. The rows are
##' named after the model states, with a numeric suffix for those that are
##' arrays (e.g. `prev_A1`, `prev_A2`).
##' @export
##' @examples
##' y <- model_simulate(example_parameters(), c(0, 7, 14), n_particles = 5)
##' y["igas_inc", , ]
model_simulate <- function(pars, times, n_particles = 1, seed = NULL,
                           deterministic = FALSE, n_threads = 1) {
  sys <- dust2::dust_system_create(model, pars, n_particles = n_particles,
                                   seed = seed, deterministic = deterministic,
                                   n_threads = n_threads,
                                   preserve_particle_dimension = TRUE)
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, times)
  rownames(y) <- model_state_names(sys)
  y
}

model_state_names <- function(sys) {
  names(unlist(dust2::dust_unpack_index(sys)))
}
