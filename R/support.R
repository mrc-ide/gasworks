
model_start_date <- function() {
  ## PHE/UKHSA report data ending for w/e on a Sunday (i.e. isoweeks)
  ## Model week 1 should end on 2014-01-05 i.e. first Sunday of 2014
  start_date <- as.Date("2014-01-05") - 7
  start_date
}

##' @name model_day
##' @title Convert date to model_day
##' @description Convert date to model_day
##' @param date A vector of dates or character strings containing dates
##' @return A vector of model_days (i.e. integers),
##' where model_day 1 = 30 Dec 2013, to align with the PHE/UKHSA definition of
##' epi weeks, whereby week 1 of 2014 ends on Sunday 5 Jan 2014
##' @examples
##' model_day("2013-12-30")
##' model_day("2014-01-05")
##' @export
model_day <- function(date) {
  start_date <- model_start_date()
  as.numeric(as.Date(date) - start_date)
}

##' @name model_date
##' @title Convert model_day to date
##' @description Convert model_day to date
##' @param day A vector of model_days (i.e. integers),
##' where model_day 1 = 30 Dec 2013, to align with the PHE/UKHSA definition of
##' epi weeks, whereby week 1 of 2014 ends on Sunday 5 Jan 2014
##' @examples
##' model_date(1)
##' model_date(7)
##' @export
model_date <- function(day) {
  model_start_date() + day
}

##' @name model_week
##' @title Convert date to model_week
##' @description Convert date to model_week
##' @param date A vector of dates or character strings containing dates
##' @return A vector of model_weeks (i.e. integers),
##' where model_week 1 = ends on Sunday 5 Jan 2014
##' to align with the PHE/UKHSA definition of epi weeks
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
##' where model_week 1 = ends on Sunday 5 Jan 2014
##' to align with the PHE/UKHSA definition of epi weeks
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
  c("U", "A", "E", "S1", "S2", "P", "F1", "F2", "R")
}

##'@name model_index
##'@title Named list of model output indices
##'@description Named list of model output indices
##'@param n_group Integer number of age groups
##'@return Named list of model output indices
model_index <- function(n_group = 1) {
  pars <- example_parameters(n_group)
  mod <- model$new(pars, 1, 1)
  idx <- mod$info()$index
  ret <- unlist(idx)
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
ll_norm <- function(data, model, sd) {
  ret <- dnorm(data, model, sd, log = TRUE)
  ret[is.na(ret)] <- 0
  ret
}

##' @title ll_multinom
##' @importFrom stats dmultinom
##' @param data a vector containing observations from a multinomial distribution
##' @param prob a matrix containing sets of probabilities for the multinomial
##' distribution, with one set per column. Values will scale automatically to
##' so that each column sums to 1. the number of rows should be the same length
##' as the data.
##' @param noise exponential noise for case when all probabilities are 0.
ll_multinom <- function(data, prob, noise) {
  stopifnot(nrow(prob) == length(data))
  if (any(is.na(data))) {
    return(numeric(ncol(prob)))
  }
  prob[is.na(prob)] <- 0
  ## need to return -Inf when p all NA / 0
  apply(prob, 2, function(p) dmultinom(data, prob = p + noise, log = TRUE))
}

##' @title ll_dirichlet
##' @importFrom extraDistr ddirichlet
##' @param data a vector containing observations from a Dirichlet distribution
##' i.e. a vector of probabilities that sum to 1
##' @param state a matrix containing sets of shape parameters, with one row
##' per parameter set
ll_dirichlet <- function(data, state, exp_noise) {
  stopifnot(ncol(state) == length(data))
  if (any(is.na(data))) {
    return(numeric(nrow(state)))
  }
  stopifnot(all.equal(sum(data), 1))
  alpha <- state + rexp(length(state), rate = exp_noise)
  ## need to return -Inf when p all NA / 0
  extraDistr::ddirichlet(data, alpha, TRUE)
}

##' @title age_spline_polynomial
##' @param age a vector of ages
##' @param pars a vector of length three giving the polynomial coefficients
##' @return polynomial spline at input ages
##' @export
age_spline_polynomial <- function(age, pars) {
  assert_length(pars, 3)
  assert_nonnegative(age)
  exp(pars[1] + log(age + 1) * pars[2] + log(age + 1) ^ 2 * pars[3])
}

##' @title age_spline_lognormal
##' @param age a vector of ages
##' @param pars a vector of length three giving the lognormal coefficients
##' mode, sd, peak (i.e. height of peak at the mode), mode and sd should be > 0,
##' peak should be in unit interval.
##' @return lognormal-based spline at input ages
##' @export
##' @importFrom stats dlnorm
age_spline_lognormal <- function(age, pars) {
  assert_length(pars, 3)
  assert_nonnegative(age)

  mode  <- pars[1]
  sigma <- pars[2]
  peak  <- pars[3]

  assert_positive(mode)
  assert_positive(sigma)
  assert_unit_interval(peak)

  # mode of lognormal = exp(mu - sigma ^ 2)
  mu <- log(mode) + sigma ^ 2

  peak * exp(dlnorm(age, mu, sigma, TRUE) - dlnorm(mode, mu, sigma, TRUE))
}

##' @title age_spline_gamma
##' @param age a vector of ages
##' @param pars a vector of length three giving the gamma coefficients
##' mode, shape, peak (i.e. height of peak at the mode), mode should be > 0,
##' shape > 1 and peak should be in unit interval.
##' @return gamma-based spline at input ages
##' @export
##' @importFrom stats dgamma
age_spline_gamma <- function(age, pars) {
  assert_length(pars, 3)
  assert_nonnegative(age)

  mode  <- pars[1]
  shape <- pars[2]
  peak  <- pars[3]

  assert_positive(mode)
  stopifnot(shape > 1)
  assert_unit_interval(peak)

  ## mode of gamma = (shape - 1) / rate
  rate <- (shape - 1) / mode
  peak * exp(dgamma(age, shape, rate, log = TRUE) -
               dgamma(mode, shape, rate, log = TRUE))
}

##' @title mean_age_spline
##' @param age_start a vector denoting the start of each age bracket
##' @param age_end a vector denoting the end of each age bracket
##' @param pars a vector of length three giving the coefficients for `spline`
##' @param spline a spline function, options include `age_spline_polynomial`,
##' `age_spline_gamma` or `age_spline_lognormal`
##' @return vector giving mean of spline output for ages within each bracket
##' @export
mean_age_spline <- function(age_start, age_end, pars, spline) {

  stopifnot(length(age_start) == length(age_end))
  stopifnot(all(age_start < age_end))
  if (any(age_start[-1] == age_end[-length(age_end)])) {
    stop("age brackets must not overlap")
  }

  f <- function(age_start, age_end) {
    age <- seq(age_start, age_end)
    mean(spline(age, pars))
  }
  mapply(f, age_start, age_end)
}
