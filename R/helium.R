##' @name helium_age_groups
##' @title Helium model age groups
##' @description Helium model age groups
##' @return A list of named model parameters
##' @export
helium_age_groups <- function() {
  age <- c(seq(0, 75, 5), 100)

  # map to UKHSA age groupings
  idx <- list(
    "04"    = 1,     #[0,5)
    "05_14" = 2:3,   #[5,10), [10,15)
    "15_44" = 4:9,   #[15,20), [20,25), [25,30), [30,35), [35,40), [40,45)
    "45_64" = 10:13, #[45,50), [50,55), [55,60), [60,65)
    "65_74" = 14:15, #[65,70), [70,75)
    "75"    = 16)    #[75)

  list(age_start = age[-length(age)],
       age_end = age[-1],
       n_group = length(age) - 1,
       idx_ukhsa = idx)
}

ukhsa_age_groups <- function() {
  age <- c(0, 5, 15, 45, 65, 75, 100)
  list(age_start = age[-length(age)],
       age_end = age[-1] - 1L,
       n_group = length(age) - 1)
}

example_helium_parameters <- function() {
  groups <- helium_age_groups()
  pars <- example_gas_parameters(groups$n_group)
  spline_pars <- list(b0_phi_S = 5, b1_phi_S = 2, b2_phi_S = 0.2,
                      b0_prev_A = 10, b1_prev_A = 3, b2_prev_A = 0.3,
                      b0_prev_R = 75, b1_prev_R = 2, b2_prev_R = 0.9)
  pars$prev_A <- pars$prev_R <- pars$phi_S <-  NULL
  pars$q_F_2014 <- 0.3
  pars$q_F_2020 <- 0.6
  c(pars, spline_pars)
}

##' @name helium_fitted_states
##' @title Helium fitted states
##' @return A character vector of data states to which the model is fitted
##' @export
helium_fitted_states <- function() {
  groups <- helium_age_groups()
  c("scarlet_fever_cases", "igas_inc",
    paste0("daily_pharyngitis_rate_", names(groups$idx_ukhsa)),
    paste0("daily_scarlet_fever_rate_", names(groups$idx_ukhsa)))
}

##' @name helium_index
##' @title Helium fitted states
##' @inheritParams hydrogen_index
##' @return A list with elements `run`, indicating the locations of the
##' compartments used in fitting the model, and `state` indicating the locations
##' of additional interesting model outputs
##' @export
##' @examples
##' p <- example_parameters(16)
##' mod <- model$new(p, 0, 10)
##' helium_index(mod$info())
helium_index <- function(info) {
  groups <- helium_age_groups()
  stopifnot(info$dim$N == groups$n_group)
  run <- c("scarlet_fever_cases", "igas_inc", "daily_pharyngitis_rate",
           paste0("daily_pharyngitis_rate_", names(groups$idx_ukhsa)),
           paste0("daily_scarlet_fever_rate_", names(groups$idx_ukhsa)))
  save <- c("prev_R", "prev_A",  "births_inc", "net_leavers_inc",
            "infections_inc", "gas_pharyngitis_inc", "scarllet_fever_inc",
            "daily_pharyngitis_rate", "daily_scarlet_fever_rate")

  list(run = unlist(info$index[run]),
       state = unlist(info$index[c(run, save)]))
}

##' Compare observed and modelled data from the helium model. This
##' conforms to the mcstate interface.
##' @title Compare observed and modelled data for the helium model
##' @param state State vector for the end of the current week. This is
##' assumed to be filtered following [helium_index()] so contains
##' rows corresponding to average daily pharyngitis rate per 100,000,
##' weekly number of scarlet fever cases and weekly number of iGAS cases.
##' @param observed Observed data. This will be a list with elements:
##' `scarlet_fever_cases`, `igas_inc`,
##' `daily_pharyngitis_rate_04`, `daily_pharyngitis_rate_05_14`,
##' `daily_pharyngitis_rate_15_44`, `daily_pharyngitis_rate_45_64`,
##' `daily_pharyngitis_rate_65_74`, `daily_pharyngitis_rate_75`,
##' `daily_scarlet_fever_rate_04`, `daily_scarlet_fever_rate_05_14`,
##' `daily_scarlet_fever_rate_15_44`, `daily_scarlet_fever_rate_45_64`,
##' `daily_scarlet_fever_rate_65_74`, `daily_scarlet_fever_rate_75`
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##' @export
helium_compare <- function(state, observed, pars) {
  groups <- helium_age_groups()

  stopifnot(pars$n_group == groups$n_group)
  if (!all(helium_fitted_states() %in% names(observed))) {
    stop("missing or misnamed data")
  }

  age_groups <- names(groups$idx_ukhsa)

  obs_sf_rate <- unlist(
    observed[paste0("daily_scarlet_fever_rate_", age_groups)])
  obs_pharyngitis_rate <- unlist(
    observed[paste0("daily_pharyngitis_rate_", age_groups)])

  model_sf_rate <- state[sprintf("daily_scarlet_fever_rate_%s", age_groups), ]
  model_pharyngitis_rate <-
    state[sprintf("daily_pharyngitis_rate_%s", age_groups), ]

  ## age-split only available between W38 2016-2018 W16
  ## need to fit to overall rate in otherwise
  if (all(is.na(obs_pharyngitis_rate))) {
    ll_pharyngitis <- ll_norm(observed$daily_pharyngitis_rate,
                             state["daily_pharyngitis_rate", ], pars$k_gp)
  } else {
    ll_pharyngitis <- colSums(ll_norm(obs_pharyngitis_rate,
                              model_pharyngitis_rate, pars$k_gp))
  }

  ll_sf_cases <- ll_nbinom(observed$scarlet_fever_cases,
                                state["scarlet_fever_cases", ],
                                pars$k_hpr, pars$exp_noise)

  ll_sf_rate <- colSums(ll_norm(obs_sf_rate, model_sf_rate, pars$k_gp))

  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ],
                       pars$k_hpr, pars$exp_noise)

  ll <- ll_pharyngitis + ll_sf_rate + ll_sf_cases + ll_igas
  ll
}

##' @title Prepare particle filter data for the helium model
##' @param data The data set to be used for the particle filter,
##' This is essentially a [data.frame()] with at least columns `model_week`,
##' along with the data used in [helium_compare()].
##' @return A a [data.frame()] with columns `step_start`, `step_end`,
##' along with the data used in [helium_compare()].
##' @importFrom mcstate particle_filter_data
helium_prepare_data <- function(data) {
  time <- "model_week"
  states <- helium_fitted_states()
  fitted <- data[, c(time, states)]
  mcstate::particle_filter_data(fitted, time = time, rate = 7)
}

##' @title Create particle filter for helium model
##' @inheritParams helium_prepare_data
##' @param constant_data a list with entries `etiologic_fraction` and
##' `asymptomatic_carriage` each containing a data.frame with entries
##' `age_start`, `age_end`, `N` and `n`.
##' @param n_particles The number of particles to simulate
##' @return a particle filter for the helium model
##' @export
helium_filter <- function(data, constant_data, n_particles) {

  data <- helium_prepare_data(data)
  constant_ll <- create_constant_log_likelihood(constant_data)
  mcstate::particle_filter$new(data, model, n_particles,
                               compare = helium_compare,
                               index = helium_index,
                               constant_log_likelihood = constant_ll)
}


create_constant_log_likelihood <- function(data) {

  constant_log_likelihood <- function(pars) {

    coefs <- helium_get_spline_coefficients(c("phi_S", "prev_A"), pars)

    loglik_phi_S <- with(data$etiologic_fraction, {
      dbinom(n, N, mean_age_spline(age_start, age_end, coefs$phi_S,
                                   age_spline_gamma), TRUE)})
    loglik_prev_A <- with(data$asymptomatic_carriage, {
      dbinom(n, N, mean_age_spline(age_start, age_end, coefs$prev_A,
                                   age_spline_gamma), TRUE)})

    sum(loglik_phi_S) + sum(loglik_prev_A)
  }

  constant_log_likelihood
}


helium_get_spline_coefficients <- function(name, pars, prefix = "b") {
  ret <- lapply(name, function(nm) {
    nms <- grep(sprintf("%s[0-9+]_%s", prefix, nm), names(pars))
    unlist(pars[nms])
  })
  names(ret) <- name
  ret
}

##' @title Create transform function for helium model
##' @inheritParams hydrogen_create_transform
##' @return A transform function for use in [`mcstate::pmcmc()`]
##' @export
helium_create_transform <- function(demographic_pars) {

  transform <- function(pars) {
    pars <- as.list(pars)

    groups <- helium_age_groups()
    pars$n_group <- groups$n_group
    ukhsa <- ukhsa_age_groups()

    # calculate age splines using gamma method
    nms <- c("phi_S", "prev_A", "prev_R")
    coefs <- helium_get_spline_coefficients(nms, pars)
    pars[nms] <- lapply(coefs, mean_age_spline, age_start = groups$age_start,
                        age_end = groups$age_end - 1, spline = age_spline_gamma)

    ## time-varying rate of reporting SF
    pars$q_F_t <- seq(model_day("2014-01-01"), model_day("2020-01-01"))
    pars$q_F <- seq(pars$q_F_2014, pars$q_F_2020,
                    length.out = length(pars$q_F_t))

    model_parameters(pars, demographic_pars = demographic_pars)
  }

  transform
}
