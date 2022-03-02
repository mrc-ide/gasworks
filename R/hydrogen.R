##' @name hydrogen_fitted_states
##' @title Hydrogen fitted states
##' @return A character vector of the model states that are fitted to data
##' @export
hydrogen_fitted_states <- function() {
  c("daily_pharyngitis_rate", "scarlet_fever_inc", "igas_inc")
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
##' p <- example_gas_parameters(1)
##' mod <- model$new(p, 0, 10)
##' hydrogen_index(mod$info())
hydrogen_index <- function(info) {
  stopifnot(info$dim$N == 1)
  run <- hydrogen_fitted_states()
  save <- c("prev_R", "prev_A", "N", "births_inc", "net_leavers_inc",
            "infections_inc", "gas_pharyngitis_inc",
            "daily_pharyngitis_scarlet_fever_rate",
            "daily_scarlet_fever_rate")

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
##' `daily_pharyngitis_rate`, `scarlet_fever_inc` and `igas_inc`
##' @param pars A list of parameters, as created by [model_parameters()]
##' @return A vector of log likelihoods, the same length as the number
##'   of particles (the number of columns in the modelled state)
##' @export
##' @examples
##' state <- rbind(daily_pharyngitis_rate = 10:15,
##'                scarlet_fever_inc = 100:105,
##'                igas_inc = 90:95)
##' observed <- list(daily_pharyngitis_rate = 13, scarlet_fever_inc = 103,
##'                  igas_inc = 93)
##' pars <- example_gas_parameters(1)
##' hydrogen_compare(state, observed, pars)
##' hydrogen_compare(state * 5, observed, pars)
hydrogen_compare <- function(state, observed, pars) {

  stopifnot(pars$n_group == 1)
  if (!all(rownames(state) %in% names(observed))) {
    stop("missing or misnamed data")
  }

  ll_pharyngitis <- ll_norm(observed$daily_pharyngitis_rate,
                            state["daily_pharyngitis_rate", ],
                            pars$k_gp, pars$exp_noise)
  ll_scarlet_fever <- ll_nbinom(observed$scarlet_fever_inc,
                                state["scarlet_fever_inc", ],
                                pars$k_hpr, pars$exp_noise)
  ll_igas <- ll_nbinom(observed$igas_inc, state["igas_inc", ],
                       pars$k_hpr, pars$exp_noise)
  ll_pharyngitis + ll_scarlet_fever + ll_igas
}
