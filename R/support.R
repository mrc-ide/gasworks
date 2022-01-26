
##' @name model_day
##' @title Convert date to model_day
##' @description Convert date to model_day
##' @param date A vector of dates or character strings containing dates
##' @return A vector of model_days (i.e. integers),
##' where model_day 1 = 1 Jan 2014
##' @examples
##' model_day("2014-01-01")
##' model_day("2014-01-07")
##' @export
model_day <- function(date) {
  start_date <- as.Date("2013-12-31")
  as.numeric(as.Date(date) - start_date)
}

##' @name model_date
##' @title Convert model_day to date
##' @description Convert model_day to date
##' @param day A vector of model_days (i.e. integers),
##'  where model_day 1 = 1 Jan 2014
##' @return A vector of dates
##' @examples
##' model_date(1)
##' model_date(7)
##' @export
model_date <- function(day) {
  as.Date("2013-12-31") + day
}

##' @name model_week
##' @title Convert date to model_week
##' @description Convert date to model_week
##' @param date A vector of dates or character strings containing dates
##' @return A vector of model_weeks (i.e. integers),
##'  where model_week 1 = the week ending 7 Jan 2014
##' @examples
##' model_week("2014-01-01")
##' model_week("2014-01-07")
##' @export
model_week <- function(date) {
  ceiling(model_day(date) / 7)
}

##' @name model_week_date
##' @title Convert model_week to date on which the week ends
##' @description Convert model_week to date on which the week ends
##' @param week A vector of model_weeks (i.e. integers),
##'  where model_week 1 ends on 7 Jan 2014
##' @return A vector of dates giving the date on which the week ends
##' @examples
##' model_week_date(1)
##' model_week_date(7)
##' @export
model_week_date <- function(week) {
  model_date(week * 7)
}

##'@name model_compartments
##'@title Names of model compartments
##'@description Names of model compartments
##'@return Names of model compartments
model_compartments <- function() {
  c("U", "A", "E", "S", "P", "F", "R")
}

##'@name model_index
##'@title Named list of model output indices
##'@description Named list of model output indices
##'@param n_group Integer number of age groups
##'@return Named list of model output indices
model_index <- function(n_group = 1) {
  pars <- example_gas_parameters(n_group)
  mod <- model$new(pars, 1, 1)
  idx <- mod$info()$index
  ret <- unlist(idx)

  if (n_group > 1) {
    suffix <- sprintf("_%02d", unlist(lapply(idx, seq_along)))
    suffix[which(lengths(idx) == 1)] <- ""
    nms <- paste0(rep(names(idx), lengths(idx)), suffix)
    names(ret) <- nms
  }

  ret
}

##' @importFrom stats dnbinom rexp
ll_nbinom <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  mu <- model + rexp(length(model), rate = exp_noise)
  # var(x) =  mu + mu ^ 2 * kappa, so kappa = 0 is Poisson,
  # large kappa is over-dispersed
  dnbinom(data, size = 1 / kappa, mu = mu, log = TRUE)
}


##' @importFrom stats dnorm
ll_norm <- function(data, model, kappa, exp_noise) {
  if (is.na(data)) {
    return(numeric(length(model)))
  }
  sd <- sqrt(model * (1 + kappa)) + rexp(length(model), exp_noise)
  # var(x) =  mu * (1 + kappa), so kappa = 0 is Poisson,
  # large kappa is over-dispersed
  dnorm(data, model, sd, log = TRUE)
}
