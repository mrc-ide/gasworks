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
    assert_positive(delta_I, 1)
    assert_positive(delta_S, 1)
    assert_positive(delta_F, 1)
    assert_positive(delta_R, 1)
    assert_unit_interval(theta_A, 1)
    assert_nonnegative_integer(alpha, 1)
    assert_unit_interval(phi_S, n_group)
    assert_nonnegative_integer(U0, n_group)
    assert_nonnegative_integer(A0, n_group)
    assert_nonnegative_integer(E0, n_group)
    assert_nonnegative_integer(I0, n_group)
    assert_nonnegative_integer(S10, n_group)
    assert_nonnegative_integer(S20, n_group)
    assert_nonnegative_integer(F0, n_group)
    assert_nonnegative_integer(R0, n_group)
  })
}


##' @name example_gas_parameters
##' @title Example of fitted gas parameters for use in testing
##' @description Example of fitted gas parameters for use in testing
##' @inheritParams check_gas_parameters
##' @return A list of named model parameters
example_gas_parameters <- function(n_group = 1) {
  pars <- list(prev_A = rep(0.1, n_group),
               prev_R = rep(0.5, n_group),
               n_group = n_group,
               beta = 2,
               sigma = 0.6,
               t_s = 100,
               p_S = 0.6,
               p_R = 0.3,
               p_I = 0.0001,
               p_F = 0.001,
               delta_A = 30,
               delta_R = 365 * 5,
               k_gp = 1,
               k_hpr = 1,
               theta_A = 1,
               phi_S = rep(0.25, n_group))
  transform(pars)
}

##' @name no_gas_parameters
##' @title Model parameters with no gas for use in testing
##' @description Model parameters with no gas for use in testing
##' @inheritParams example_gas_parameters
##' @return A list of named model parameters
no_gas_parameters <- function(n_group = 1) {
  pars <- list(m = diag(0, n_group),
               prev_A = rep(0, n_group),
               prev_R = rep(0, n_group),
               n_group = n_group,
               beta = 0,
               sigma = 0,
               t_s = 1,
               p_S = 0,
               p_R = 0,
               p_I = 0,
               p_F = 0,
               delta_A = 1,
               delta_R = 1,
               delta_E = 1,
               delta_I = 1,
               delta_S = 1,
               delta_F = 1,
               theta_A = 0,
               phi_S = rep(0, n_group))
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
  n_group <- gas_pars$n_group %||% 1
  demographic_pars <- demographic_pars %||% demographic_parameters(n_group)
  pars <- c(demographic_pars, gas_pars)

  pars$dt <- 1 / 7 # Data is weekly, model steps are daily
  pars$exp_noise <- 1e6 # exponential noise parameter for observation dist

  # add fixed model parameters (i.e. not fitted)
  pars$delta_E <- 2   # mean days in incubation period
  pars$delta_I <- 14  # mean days iGAS
  pars$delta_S <- 2.3 # mean days with pharyngitis symptoms (x 2)
  pars$delta_F <- 7   # mean days with scarlet fever

  # convert duration in days to duration in weeks
  for (i in grep("^delta_", names(pars))) {
    pars[[i]] <- pars[[i]] * pars$dt
  }

  initial_pars <- initial_pars %||% initial_parameters(pars)
  pars <- c(pars, initial_pars)
  check_gas_parameters(pars)
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
  nms <- paste0(model_compartments(), "0")
  ret <- named_list(nms, list(rep(0, n_group)))
  # set initial asymptomatic prevalence and prevalence of immunity
  ret$A0 <- round(pars$N0 * pars$prev_A)
  ret$R0 <- round((pars$N0 - ret$A0) * pars$prev_R)

  # set initial uninfecteds
  ret$U0 <- pars$N0 - ret$A0 - ret$R0

  ret
}

##' Transform fitted parameters into gas params
##' @name transform
##' @title Transform fitted parameters into gas params
##' @description Transform fitted parameters into gas params
##' @param pars list of fitted parameters
##' @return A list of parameters for use in the model
##' @export
transform <- function(pars) {
  pars <- as.list(pars)
  if (!is.null(pars$dt)) {
    stop("Parameters have already been transformed")
  }
  # add initial conditions and demographic parameters
  model_parameters(pars)
}
