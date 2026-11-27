##' @name hydrogen_fitted_states
##' @title Hydrogen fitted states
##' @return A character vector of the model states that are fitted to data
##' @export
hydrogen_fitted_states <- function() {
  c("daily_pharyngitis_rate", "daily_scarlet_fever_rate",
    "scarlet_fever_cases", "igas_inc")
}

##' @name hydrogen_saved_states
##' @title Hydrogen saved states
##' @return A character vector of the model states that are of interest when
##' fitting the hydrogen model: those compared with data, followed by some
##' additional outputs. This can be passed as `save_trajectories` to
##' [dust2::dust_likelihood_monty()].
##' @export
##' @examples
##' hydrogen_saved_states()
hydrogen_saved_states <- function() {
  c("daily_gas_pharyngitis_rate", "daily_scarlet_fever_rate",
    "scarlet_fever_cases", "igas_inc", "etiologic_fraction",
    "prev_R", "prev_A", "N", "births_inc", "net_leavers_inc",
    "infections_inc", "gas_pharyngitis_inc", "scarlet_fever_inc")
}

##' Compare observed and modelled data from the hydrogen model. This is a
##' reference implementation in R of the comparison that is compiled into
##' the [model], which is what is used by [hydrogen_filter()].
##' @title Compare observed and modelled data for the hydrogen model
##' @param state A matrix of model states at the end of the current week,
##' with one column per particle and rows named
##' `daily_gas_pharyngitis_rate`, `daily_scarlet_fever_rate`,
##' `scarlet_fever_cases` and `igas_inc`.
##' @param observed Observed data. This will be a list with elements
##' `daily_pharyngitis_rate`, `daily_scarlet_fever_rate`, `scarlet_fever_cases`
##' and `igas_inc`
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##' @export
##' @examples
##' state <- rbind(daily_gas_pharyngitis_rate = 10:15,
##'                daily_scarlet_fever_rate = 0:5,
##'                scarlet_fever_cases = 100:105,
##'                igas_inc = 90:95)
##' observed <- list(daily_pharyngitis_rate = 13,
##'                  daily_scarlet_fever_rate = 3,
##'                  scarlet_fever_cases = 103,
##'                  igas_inc = 93)
##' pars <- example_parameters(1)
##' hydrogen_compare(state, observed, pars)
##' hydrogen_compare(state * 5, observed, pars)
hydrogen_compare <- function(state, observed, pars) {

  stopifnot(pars$n_group == 1)
  if (!all(hydrogen_fitted_states() %in% names(observed))) {
    stop("missing or misnamed data")
  }

  ll_pharyngitis <- ll_norm(observed$daily_pharyngitis_rate * pars$phi_S,
                            state["daily_gas_pharyngitis_rate", ], pars$k_gp)

  ll_sf_cases <- ll_nbinom(observed$scarlet_fever_cases,
                           state["scarlet_fever_cases", ],
                           pars$k_hpr, pars$exp_noise)

  ll_sf_rate <- ll_norm(observed$daily_scarlet_fever_rate,
                        state["daily_scarlet_fever_rate", ], pars$k_gp)

  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ],
                       pars$k_hpr, pars$exp_noise)

  ll_pharyngitis + ll_sf_cases + ll_sf_rate + ll_igas
}

##' @title Prepare particle filter data for the hydrogen model
##' @param data The data set to be used for the particle filter,
##' This is essentially a [data.frame()] with at least columns `model_week`,
##' along with the data used in [hydrogen_compare()].
##' @return A `dust_filter_data` [data.frame()] with a column `time` (in
##' days), and one column for each data stream known to the [model], prefixed
##' with `obs_`. Data streams that the hydrogen model is not fitted to are
##' set to `NA`.
##' @export
hydrogen_prepare_data <- function(data) {
  prepare_data(data, hydrogen_fitted_states())
}

##' @title Create particle filter for hydrogen model
##' @inheritParams hydrogen_prepare_data
##' @param n_particles The number of particles to simulate
##' @param n_threads The number of threads to use
##' @param seed The random number seed, see [dust2::dust_filter_create()]
##' @param deterministic Logical, indicating if a deterministic likelihood
##' should be created with [dust2::dust_unfilter_create()], in which case
##' `n_particles` and `seed` are ignored
##' @return a particle filter for the hydrogen model, which can be run with
##' [dust2::dust_likelihood_run()]
##' @export
hydrogen_filter <- function(data, n_particles, n_threads = 1, seed = NULL,
                            deterministic = FALSE) {
  data <- hydrogen_prepare_data(data)
  create_filter(data, n_particles, n_threads, seed, deterministic)
}

##' @title Create transform function for hydrogen model
##' @param demographic_pars A list of demographic parameters containing elements
##' `N0`, `alpha`, `omega`, `m`, `r_age` to be loaded into the transform
##' @return A function that converts a list of fitted parameters into the
##' list of parameters required by the [model]
##' @export
hydrogen_create_transform <- function(demographic_pars) {
  ## evaluate now, so that the transform does not depend on the calling
  ## environment (e.g. when it is sent to another process)
  force(demographic_pars)
  transform <- function(pars) {
    pars <- as.list(pars)

    ## time-varying rate of reporting SF
    pars$q_F <- reporting_probability(pars$q_F_2014, pars$q_F_2020)

    model_parameters(pars, demographic_pars = demographic_pars)
  }
  transform
}

##' @title Create parameter packer for hydrogen model
##' @param names A character vector of the names of the fitted parameters
##' @inheritParams hydrogen_create_transform
##' @param fixed An optional named list of parameters that are not fitted
##' @return A [monty::monty_packer()] that unpacks a vector of fitted
##' parameters into the list of parameters required by the [model], for use
##' with [dust2::dust_likelihood_monty()]
##' @export
hydrogen_packer <- function(names, demographic_pars, fixed = NULL) {
  create_packer(names, hydrogen_create_transform(demographic_pars), fixed)
}
