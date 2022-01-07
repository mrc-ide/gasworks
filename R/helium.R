
##' @name helium_age_groups
##' @title Helium model age groups
##' @description Helium model age groups
##' @return A list of named model parameters
##' @export
helium_age_groups <- function() {
  age <- c(seq(0, 90, 5), 100)
  list(age_start = age[-length(age)],
       age_end = age[-1],
       n_group = length(age) - 1)
}
