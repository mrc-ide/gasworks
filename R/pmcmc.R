
##' @title Calculate the log likelihood of the data given the parameters
##' @description Calculate the log likelihood of the data given the parameters
##' @param state State vector for the end of the current day. This is
##' assumed to be filtered following [index()] so contains 3 rows corresponding
##' to daily_pharyngitis_scarlet_fever_rate, scarlet_fever and igas flows
##' @param observed Observed data containing entries
##'  `daily_pharyngitis_scarlet_fever_rate`, `scarlet_fever_inc` and `igas_inc`.
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return a single log likelihood
##' @export
compare <- function(state, observed, pars) {

  if (!all(rownames(state) %in% names(observed))) {
    stop("missing or misnamed data")
  }

  ## use k_gp for GP surveillance data
  ## rate is a continuous dist - need to use a normal
  ll_pharyngitis_scarlet_fever <-
    ll_norm(observed$daily_pharyngitis_scarlet_fever_rate,
            state["daily_pharyngitis_scarlet_fever_rate", ],
            pars$k_gp, pars$exp_noise)

  ## use k_hpr for Health Protection Report data

  ll_scarlet_fever <- ll_nbinom(observed$scarlet_fever_inc,
                                state["scarlet_fever_inc", ],
                                pars$k_hpr, pars$exp_noise)

  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ],
                       pars$k_hpr, pars$exp_noise)

  ll_pharyngitis_scarlet_fever + ll_scarlet_fever + ll_igas
}
