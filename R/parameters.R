##' @name check_gas_parameters
##' @title Check model parameters are within permitted ranges
##' @description Check model parameters are within permitted ranges
##' @param pars a named list of model parameters
##' @param n_group the number of model age-groups, defaults to 1
##' @return NULL
##' @export

check_gas_parameters <- function(pars, n_group = 1) {
  with(pars, {
    assert_nonnegative(beta, 1)
    assert_nonnegative(sigma, 1)
    assert_positive(t_s, 1)
    assert_unit_interval(p_S, 1)
    assert_unit_interval(p_R, 1)
    assert_unit_interval(p_I, 1)
    assert_unit_interval(p_F, 1)
    assert_positive(delta_A, 1)
    assert_positive(delta_E, 1)
    assert_positive(delta_S, 1)
    assert_positive(delta_P, 1)
    assert_positive(delta_F, 1)
    assert_positive(delta_R, 1)
    assert_positive_integer(k_A, 1)
    assert_positive_integer(k_E, 1)
    assert_positive_integer(k_S, 1)
    assert_positive_integer(k_P, 1)
    assert_positive_integer(k_F, 1)
    assert_positive_integer(k_R, 1)
    assert_unit_interval(theta_A, 1)
    assert_nonnegative_integer(alpha)
    assert_nonnegative_integer(U0, n_group)
    assert_nonnegative_integer(A0, n_group * k_A)
    assert_nonnegative_integer(E0, n_group * k_E)
    assert_nonnegative_integer(S0, n_group * k_S)
    assert_nonnegative_integer(P0, n_group * k_P)
    assert_nonnegative_integer(F0, n_group * k_F)
    assert_nonnegative_integer(R0, n_group * k_R)
    assert_dim(m, c(n_group, n_group))
    assert_nonnegative(m)
    assert_nonnegative(r_age, 1)
    assert_length(omega, n_group)
  })
}

example_gas_parameters <- function(n_group = 1) {
  list(prev_A = rep(0.1, n_group),
       prev_R = rep(0.5, n_group),
       n_group = n_group,
       beta = 2,
       sigma = 0.6,
       t_s = 100,
       p_S = 0.6,
       p_R = 0.3,
       p_I = 0.0001,
       p_F = 0.001,
       p_T = 1,
       delta_A = 30,
       delta_R = 365 * 5,
       k_gp = 1,
       k_hpr = 1,
       theta_A = 1,
       phi_S = rep(0.25, n_group),
       q_F = 0.3)
}

##' @name example_parameters
##' @title Example of transformed model parameters for use in testing
##' @inheritParams check_gas_parameters
##' @return A list of named model parameters
##' @export
example_parameters <- function(n_group = 1) {
  pars <- example_gas_parameters(n_group)
  model_parameters(pars)
}

##' @name no_gas_parameters
##' @title Model parameters with no gas for use in testing
##' @description Model parameters with no gas for use in testing
##' @inheritParams example_parameters
##' @return A list of named model parameters
no_gas_parameters <- function(n_group = 1) {
  pars <- list(prev_A = rep(0, n_group),
               prev_R = rep(0, n_group),
               n_group = n_group,
               beta = 0,
               sigma = 0,
               t_s = 1,
               p_S = 0,
               p_R = 0,
               p_I = 0,
               p_F = 0,
               p_T = 0,
               delta_A = 1,
               delta_R = 1,
               delta_E = 1,
               delta_I = 1,
               delta_S = 1,
               delta_F = 1,
               theta_A = 0,
               phi_S = rep(1, n_group),
               q_F = 0)
}

##' @name model_parameters
##' @title Demographic model parameters
##' @description Demographic model parameters
##' @param gas_pars A list of named fitted model parameters
##' @param initial_pars default NULL uses `N0`, `prev_A` and `prev_R` to
##' calculate the starting position
##' @param demographic_pars default NULL uses `demographic_parameters()`)
##' @inheritParams check_gas_parameters
##' @return A list of all model parameters
##' @export
model_parameters <- function(gas_pars, initial_pars = NULL,
                             demographic_pars = NULL, n_group = NULL) {
  n_group <- n_group %||% (gas_pars$n_group  %||% 1)
  demographic_pars <- demographic_pars %||% demographic_parameters(n_group)
  pars <- c(demographic_pars, gas_pars)

  pars$dt <- 1 / 7 # Data is weekly, model steps are daily
  pars$exp_noise <- 1e6 # exponential noise parameter for observation dist

  # add fixed model parameters (i.e. not fitted)
  pars$delta_E <- 2.1 # mean days in incubation period
  pars$delta_S <- 4.3 # mean days with pharyngitis symptoms
  pars$delta_P <- 1.6 # mean days pharyngitis -> scarlet fever
  pars$delta_F <- 6   # mean days with scarlet fever

  pars$k_A <- 1
  pars$k_E <- 1
  pars$k_S <- 2
  pars$k_P <- 1
  pars$k_F <- 2
  pars$k_R <- 1

  pars$k_gp <- 1
  pars$k_hpr <- 1

  # convert duration in days to duration in weeks
  for (i in grep("^delta_", names(pars))) {
    pars[[i]] <- pars[[i]] * pars$dt
  }

  initial_pars <- initial_pars %||% initial_parameters(pars)
  pars <- c(pars, initial_pars)
  check_gas_parameters(pars, n_group)
  pars
}

##' @name demographic_parameters
##' @title Demographic model parameters
##' @description Demographic model parameters
##' @inheritParams check_gas_parameters
##' @return A list of demographic parameters
##' @export
demographic_parameters <- function(n_group = 1) {
  N0 <- 56e6 # England population
  x <- 11000 # births and deaths per week - i.e. 572000 per year
  # convert annual mortality to weekly
  list(N0 = round(rep(N0 / n_group, n_group)),
       alpha = round(x),
       omega = rep(x / N0, n_group),
       r_age = 0,
       m = matrix(1 / n_group, n_group, n_group)) # uniform mixing matrix
}

##' @name initial_parameters
##' @title Initial conditions for the model
##' @description Initial conditions for the model
##' @param pars A parameter list containing `N0`, `prev_A` and `prev_R`.
##' @return A named list of initial model states
##' @export
initial_parameters <- function(pars) {

  # check parameters
  assert_nonnegative_integer(pars$N0)
  n_group <- length(pars$N0)
  assert_unit_interval(pars$prev_A, n_group)
  assert_unit_interval(pars$prev_R, n_group)

  ret <- list()
  # set initial asymptomatic prevalence and prevalence of immunity
  ret$A0 <- matrix(round(pars$N0 * pars$prev_A / pars$k_A), ncol = pars$k_A)
  ret$E0 <- matrix(0, n_group, pars$k_E)
  ret$S0 <- matrix(0, n_group, pars$k_S)
  ret$P0 <- matrix(0, n_group, pars$k_P)
  ret$F0 <- matrix(0, n_group, pars$k_F)
  ret$R0 <- matrix(round((pars$N0 - rowSums(ret$A0)) * pars$prev_R / pars$k_R),
                   ncol = pars$k_R)

  # set initial uninfecteds
  ret$U0 <- pars$N0 - rowSums(ret$A0) - rowSums(ret$R0)

  ret
}
