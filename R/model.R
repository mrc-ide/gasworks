##' @name model
##' @title Stochastic model of GAS
##'This is an odin.dust model.
##' @export model
NULL

##' Index of interesting elements for the model. This function
##' conforms to the mcstate interface.
##'@name index
##'@title Index of the model
##'@param info The result of running the `$info()` method on an initialised
##' model
##'@return A list with element `run`, indicating the locations of (in order)
##' (1) pharyngitis_scarlet_fever_rate, (2) scarlet_fever_inc,
##' (3) iGAS_inc and with element `state` containing the same values,
##' followed by U, R, entrants_inc, leavers_inc, infections_inc and
##' pharyngitis_inc.
##' @export
index <- function(info) {
  run <- c("daily_pharyngitis_scarlet_fever_rate",
           "scarlet_fever_inc", "igas_inc")
  save <- c("U", "R", "births_inc", "net_leavers_inc", "infections_inc",
            "pharyngitis_inc")

  list(run = unlist(info$index[run]),
       state = unlist(info$index[c(run, save)]))
}
