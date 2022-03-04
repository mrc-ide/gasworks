##' @name helium_fitted_states
##' @title Helium fitted states
##' @return A character vector of data states to which the model is fitted
##' @export
helium_fitted_states <- function() {
  groups <- helium_age_groups()
  c("scarlet_fever_inc", "igas_inc", "daily_pharyngitis_rate",
    paste0("daily_pharyngitis_rate_", names(groups$idx_ukhsa)),
    paste0("daily_scarlet_fever_rate_", names(groups$idx_ukhsa)))
}

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

##' @name helium_index
##' @title Helium fitted states
##' @inheritParams helium_index
##' @return A list with elements `run`, indicating the locations of the
##' compartments used in fitting the model, and `state` indicating the locations
##' of additional interesting model outputs
##' @export
##' @examples
##' p <- example_parameters(19)
##' mod <- model$new(p, 0, 10)
##' helium_index(mod$info())
helium_index <- function(info) {
  groups <- helium_age_groups()
  stopifnot(info$dim$N == groups$n_group)
  run <- c("scarlet_fever_inc", "igas_inc", "daily_pharyngitis_rate", "N",
           paste0("pharyngitis_prop_", names(groups$idx_ukhsa)),
           paste0("scarlet_fever_prop_", names(groups$idx_ukhsa)))
  save <- c("prev_R", "prev_A",  "births_inc", "net_leavers_inc",
            "infections_inc", "gas_pharyngitis_inc",
            "daily_scarlet_fever_rate")

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
##' `scarlet_fever_inc`, `igas_inc`, `daily_pharyngitis_rate`,
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
  pharyngitis_prop <- state[sprintf("pharyngitis_prop_%s", age_groups), ]
  scarlet_fever_prop <- state[sprintf("scarlet_fever_prop_%s", age_groups), ]

  cases <- helium_convert_incidence_to_cases(state, observed)

  ll_pharyngitis <- ll_norm(observed$daily_pharyngitis_rate,
                            state["daily_pharyngitis_rate", ], pars$k_gp) +
    ll_multinom(cases$pharyngitis, pharyngitis_prop, 1e-10)

  ll_scarlet_fever <- ll_nbinom(observed$scarlet_fever_inc,
                                state["scarlet_fever_inc", ],
                                pars$k_hpr, pars$exp_noise) +
    ll_multinom(cases$scarlet_fever, scarlet_fever_prop, 1e-10)

  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ], pars$k_hpr,
                       pars$exp_noise)

  ll_pharyngitis + ll_scarlet_fever + ll_igas
}


helium_convert_incidence_to_cases <- function(state, observed) {
  # state and observed are same as supplied to compare function
  groups <- helium_age_groups()
  # can take a single particle as all will be the same due to deterministic
  # demography
  N <- state[paste0("N", seq_len(groups$n_group)), 1]
  # convert to UKHSA age groups
  pop <- vapply(groups$idx_ukhsa, function(idx) sum(N[idx]), numeric(1))

  nms <- names(groups$idx_ukhsa)
  # extract rates by age from data
  rates <- list(
    pharyngitis = observed[paste0("daily_pharyngitis_rate_", nms)],
    scarlet_fever = observed[paste0("daily_scarlet_fever_rate_", nms)])

  # convert daily -> weekly
  # per 100,000 -> per person
  # per person -> whole pop
  cases <- lapply(rates, function(x) round(unlist(x) * pop * 7 / 1e5))
  cases
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

    # calculate age splines using gamma method
    nms <- c("phi_S", "prev_A", "prev_R")
    coefs <- helium_get_spline_coefficients(nms, pars)
    pars[nms] <- lapply(coefs, mean_age_spline, age_start = groups$age_start,
                        age_end = groups$age_end - 1, spline = age_spline_gamma)

    model_parameters(pars, demographic_pars = demographic_pars)
  }

  transform
}
