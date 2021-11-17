
##' @title Calculate the log likelihood of the data given the parameters
##' @description Calculate the log likelihood of the data given the parameters
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [model_index()] so contains
##'   3 rows corresponding to pharyngitis, scarlet_fever and igas flows
##' @param observed Observed data.
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return a single log likelihood
compare <- function(state, observed, pars) {
  idx <- model_index()

  # GP surveillance data are per 100,000 population
  pharyngitis <- state[idx$pharyngitis, ] * 1e5 / state[idx$N, ]
  scarlet_fever <- state[idx$scarlet_fever, ]
  igas <- state[idx$igas, ]

  ## TODO: incorporate % of pharyngitis attributable to GAS
  ll_pharyngitis <- dnbinom(observed$pharyngitis * pars[["phi_gas"]],
                            size = 1 / pars[["k_gp"]],
                            mu = pharyngitis, log = TRUE)
  ll_scarlet_fever <- dnbinom(observed$scarlet_fever, size = 1 / pars[["k_gp"]],
                              mu = scarlet_fever, log = TRUE)
  ll_igas <- dnbinom(observed$igas, size = 1 / pars[["k_hpr"]],
                     mu = igas, log = TRUE)

  ll_pharyngitis + ll_scarlet_fever + ll_igas
}
