##' @name hydrogen_fitted_states
##' @title Hydrogen fitted states
##' @return A character vector of the model states that are fitted to data
##' @export
hydrogen_fitted_states <- function() {
  c("daily_pharyngitis_rate", "daily_scarlet_fever_rate",
    "scarlet_fever_cases", "igas_inc")
}

##' @name hydrogen_index
##' @title Hydrogen fitted states
##' @param info The result of running the `$info()` method on an
##' initialised model
##' @return A list with elements `run`, indicating the locations of the
##' compartments used in fitting the model, and `state` indicating the locations
##' of additional interesting model outputs
##' @export
##' @examples
##' p <- example_parameters(1)
##' mod <- model$new(p, 0, 10)
##' hydrogen_index(mod$info())
hydrogen_index <- function(info) {
  stopifnot(info$dim$N == 1)
  run <-  c("daily_gas_pharyngitis_rate", "daily_scarlet_fever_rate",
            "scarlet_fever_cases", "igas_inc")
  save <- c("prev_R", "prev_A", "N", "births_inc", "net_leavers_inc",
            "infections_inc", "gas_pharyngitis_inc", "scarlet_fever_inc")

  list(run = unlist(info$index[run]),
       state = unlist(info$index[c(run, save)]))
}

##' Compare observed and modelled data from the hydrogen model. This
##' conforms to the mcstate interface.
##' @title Compare observed and modelled data for the hydrogen model
##' @param state State vector for the end of the current week. This is
##' assumed to be filtered following [hydrogen_index()] so contains
##' rows corresponding to average daily pharyngitis rate per 100,000,
##' weekly number of scarlet fever cases and weekly number of iGAS cases.
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
                            state["daily_gas_pharyngitis_rate", ] * pars$p_T,
                            pars$k_gp)
  ll_scarlet_fever <- ll_nbinom(observed$scarlet_fever_cases,
                                state["scarlet_fever_cases", ],
                                pars$k_hpr, pars$exp_noise) +
    ll_norm(observed$daily_scarlet_fever_rate,
            state["daily_scarlet_fever_rate", ], pars$k_gp)

  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ],
                       pars$k_hpr, pars$exp_noise)

  ll_pharyngitis + ll_scarlet_fever + ll_igas
}

##' @title Prepare particle filter data for the hydrogen model
##' @param data The data set to be used for the particle filter,
##' This is essentially a [data.frame()] with at least columns `model_week`,
##' along with the data used in [hydrogen_compare()].
##' @return A a [data.frame()] with columns `step_start`, `step_end`,
##' along with the data used in [hydrogen_compare()].
##' @importFrom mcstate particle_filter_data
hydrogen_prepare_data <- function(data) {
  time <- "model_week"
  states <- hydrogen_fitted_states()
  fitted <- data[, c(time, states)]
  mcstate::particle_filter_data(fitted, time = time, rate = 7)
}

##' @title Create particle filter for hydrogen model
##' @inheritParams hydrogen_prepare_data
##' @param n_particles The number of particles to simulate
##' @return a particle filter for the hydrogen model
##' @export
##' @importFrom mcstate particle_filter
hydrogen_filter <- function(data, n_particles) {
  data <- hydrogen_prepare_data(data)
  mcstate::particle_filter$new(data, model, n_particles,
                               compare = hydrogen_compare,
                               index = hydrogen_index)
}

##' @title Create transform function for hydrogen model
##' @param demographic_pars A list of demographic parameters containing elements
##' `N0`, `alpha`, `omega`, `m`, `r_age` to be loaded into the transform
##' @return A transform function for use in [`mcstate::pmcmc()`]
##' @export
hydrogen_create_transform <- function(demographic_pars) {
  transform <- function(pars) {
    pars <- as.list(pars)

    ## time-varying rate of reporting SF
    pars$q_F_t <- seq(model_day("2014-01-01"), model_day("2020-01-01"))
    pars$q_F <- seq(pars$q_F_2014, pars$q_F_2020,
                    length.out = length(pars$q_F_t))

    model_parameters(pars, demographic_pars = demographic_pars)
  }
  transform
}
