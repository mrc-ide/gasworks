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
##'@return A list with element `run`, indicating the locations of (in
##' order) (1) pharyngitis_rate, (2) scarlet_fever_inc, (3) iGAS_inc
##' and with element `state` containing the same values followed by the
##' U, R, entrants_inc and leavers_inc, infections_inc and pharyngitis_inc.
##' @export
index <- function(info) {
  index_run <- c(pharyngitis_rate = info$index$pharyngitis_rate,
                 scarlet_fever_inc = info$index$scarlet_fever_inc,
                 igas_inc = info$index$igas_inc)

  index_save <- c(U = info$index$U,
                  R = info$index$R,
                  entrants_inc = info$index$entrants_inc,
                  leavers_inc = info$index$leavers_inc,
                  infections_inc = info$index$infections_inc,
                  pharyngitis_inc = info$index$infections_inc)

  list(run = index_run,
       state = c(index_run, index_save))
}
