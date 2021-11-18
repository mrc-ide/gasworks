
##' @title Calculate the log likelihood of the data given the parameters
##' @description Calculate the log likelihood of the data given the parameters
##' @param state State vector for the end of the current day. This is
##'   assumed to be filtered following [model_index()] so contains
##'   3 rows corresponding to pharyngitis, scarlet_fever and igas flows
##' @param observed Observed data.
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return a single log likelihood
compare <- function(state, observed, pars) {
  exp_noise <- 1e6
  idx <- model_index()

  pharyngitis <- calculate_pharyngitis_incidence(state, idx, pars)
  scarlet_fever <- state[idx$scarlet_fever_inc, ]
  igas <- state[idx$igas_inc, ]

  ## continuous dist - need to use a normal, relate variance to mean
  ll_pharyngitis <- ll_norm(observed$pharyngitis ,
                            model_mean = pharyngitis,
                            model_sd = sqrt(pharyngitis) * pars[["k_gp"]])
  ll_scarlet_fever <- ll_nbinom(observed$scarlet_fever, model = scarlet_fever,
                                kappa = 1 / pars[["k_gp"]], exp_noise)
  ll_igas <- ll_nbinom(observed$igas, model = igas, kappa = 1 / pars[["k_hpr"]],
                       exp_noise)

  ll_pharyngitis + ll_scarlet_fever + ll_igas
}

calculate_pharyngitis_incidence <- function(state, idx, pars) {
  # GP surveillance data are per 100,000 population allowing for misattribution
  # with prob phi_S
  state[idx$pharyngitis_inc, ] * 1e5 / state[idx$N, ] / pars[["phi_S"]]
}
