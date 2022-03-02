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
