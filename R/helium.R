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
  c("scarlet_fever_cases", "igas_inc", "daily_pharyngitis_rate",
    paste0("daily_pharyngitis_rate_", names(groups$idx_ukhsa)),
    paste0("daily_scarlet_fever_rate_", names(groups$idx_ukhsa)))
}

##' @name helium_saved_states
##' @title Helium saved states
##' @return A character vector of the model states that are of interest when
##' fitting the helium model: those compared with data, followed by some
##' additional outputs. This can be passed as `save_trajectories` to
##' [dust2::dust_likelihood_monty()].
##' @export
##' @examples
##' helium_saved_states()
helium_saved_states <- function() {
  groups <- names(helium_age_groups()$idx_ukhsa)
  c("scarlet_fever_cases", "igas_inc", "daily_gas_pharyngitis_rate",
    "etiologic_fraction",
    paste0("etiologic_fraction_", groups),
    paste0("daily_gas_pharyngitis_rate_", groups),
    paste0("daily_scarlet_fever_rate_", groups),
    "prev_R", "prev_A", "N", "births_inc", "net_leavers_inc",
    "infections_inc", "gas_pharyngitis_inc", "scarlet_fever_inc",
    "daily_scarlet_fever_rate")
}

##' Compare observed and modelled data from the helium model. This is a
##' reference implementation in R of the comparison that is compiled into
##' the [model], which is what is used by [helium_filter()].
##' @title Compare observed and modelled data for the helium model
##' @param state A matrix of model states at the end of the current week,
##' with one column per particle and rows named after (at least) the fitted
##' states in [helium_saved_states()].
##' @param observed Observed data. This will be a list with elements:
##' `scarlet_fever_cases`, `igas_inc`, `daily_pharyngitis_rate`,
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
  obs_pharyngitis_rate <-
    unlist(observed[paste0("daily_pharyngitis_rate_", age_groups)])

  model_sf_rate <- state[sprintf("daily_scarlet_fever_rate_%s", age_groups), ]
  model_gas_pharyngitis_rate <-
    state[sprintf("daily_gas_pharyngitis_rate_%s", age_groups), ]

  if (all(is.na(obs_pharyngitis_rate))) {
    ## where age-disaggregation is not avaliable, use total rate
    phi_S <- state["etiologic_fraction", ]
    ll_pharyngitis <- ll_norm(observed$daily_pharyngitis_rate * phi_S,
                              state["daily_gas_pharyngitis_rate", ], pars$k_gp)
  } else {

    phi_S <- state[sprintf("etiologic_fraction_%s", age_groups), ]
    ll_pharyngitis <- colSums(ll_norm(obs_pharyngitis_rate * phi_S,
                                      model_gas_pharyngitis_rate, pars$k_gp))
  }

  # if q_F is time-varying then need to do this in odin
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
##' @return A `dust_filter_data` [data.frame()] with a column `time` (in
##' days), and one column for each data stream known to the [model], prefixed
##' with `obs_`. Data streams that the helium model is not fitted to are set
##' to `NA`, as is the all-age pharyngitis rate in any week where rates by
##' age are available.
##' @export
helium_prepare_data <- function(data) {
  ret <- prepare_data(data, helium_fitted_states())
  ## where age-disaggregation is not available, use total rate
  nms <- sprintf("obs_daily_pharyngitis_rate_%s", ukhsa_age_group_names())
  by_age <- rowSums(!is.na(ret[nms])) > 0
  ret$obs_daily_pharyngitis_rate[by_age] <- NA_real_
  ret
}

##' @title Create particle filter for helium model
##' @inheritParams helium_prepare_data
##' @inheritParams hydrogen_filter
##' @return a particle filter for the helium model, which can be run with
##' [dust2::dust_likelihood_run()]. Note that this does not include the
##' likelihood of the time-invariant data, see
##' [helium_constant_likelihood()].
##' @export
helium_filter <- function(data, n_particles, n_threads = 1, seed = NULL,
                          deterministic = FALSE) {
  data <- helium_prepare_data(data)
  create_filter(data, n_particles, n_threads, seed, deterministic)
}

##' @title Likelihood of time-invariant data for the helium model
##' @description The helium model is fitted to cross-sectional studies of the
##' etiologic fraction and of asymptomatic carriage by age, in addition to
##' the surveillance time series. This creates the likelihood of those
##' studies, to be combined with the particle filter likelihood.
##' @param constant_data a list with entries `etiologic_fraction` and
##' `asymptomatic_carriage` each containing a data.frame with entries
##' `age_start`, `age_end`, `N` and `n`.
##' @param packer A packer for the fitted parameters, as created by
##' [helium_packer()]
##' @return A [monty::monty_model()] over the fitted parameters, which can be
##' added to the model returned by [dust2::dust_likelihood_monty()]
##' @export
helium_constant_likelihood <- function(constant_data, packer) {
  force(packer)
  constant_log_likelihood <- create_constant_log_likelihood(constant_data)
  density <- function(x) {
    constant_log_likelihood(packer$unpack(x))
  }
  monty::monty_model(list(parameters = packer$names(), density = density))
}


##' @importFrom stats dbinom
create_constant_log_likelihood <- function(data) {
  force(data)

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
##' @return A function that converts a list of fitted parameters into the
##' list of parameters required by the [model]
##' @export
helium_create_transform <- function(demographic_pars) {
  ## evaluate now, so that the transform does not depend on the calling
  ## environment (e.g. when it is sent to another process)
  force(demographic_pars)

  transform <- function(pars) {
    pars <- as.list(pars)

    groups <- helium_age_groups()
    pars$n_group <- groups$n_group

    # calculate age splines using gamma method
    nms <- c("phi_S", "prev_A", "prev_R")
    coefs <- helium_get_spline_coefficients(nms, pars)
    pars[nms] <- lapply(coefs, mean_age_spline, age_start = groups$age_start,
                        age_end = groups$age_end - 1, spline = age_spline_gamma)

    ## time-varying rate of reporting SF
    pars$q_F <- reporting_probability(pars$q_F_2014, pars$q_F_2020)

    model_parameters(pars, demographic_pars = demographic_pars)
  }

  transform
}

##' @title Create parameter packer for helium model
##' @inheritParams hydrogen_packer
##' @return A [monty::monty_packer()] that unpacks a vector of fitted
##' parameters into the list of parameters required by the [model], for use
##' with [dust2::dust_likelihood_monty()]
##' @export
helium_packer <- function(names, demographic_pars, fixed = NULL) {
  create_packer(names, helium_create_transform(demographic_pars), fixed)
}
