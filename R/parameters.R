##' @name check_gas_parameters
##' @title Check model parameters are within permitted ranges
##' @description Check model parameters are within permitted ranges
##' @param pars a named list of model parameters
##' @param n_group the number of model age-groups, defaults to 1
##' @return NULL
##' @export

check_gas_parameters <- function(pars, n_group = 1) {
  with(pars, {
    assert_positive(beta, 1)
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
    assert_nonnegative_integer(alpha, n_group)
    assert_unit_interval(omega, n_group)
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
##' @return A list of named model parameters
example_gas_parameters <- function() {
  pars <- list(prev_A = 0.1,
               prev_R = 0.5,
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
               phi_S = 0.25)
  transform(pars)
}


##' @name model_parameters
##' @title Demographic model parameters
##' @description Demographic model parameters
##' @param gas_pars A list of named fitted model parameters
##' @param initial_pars default NULL uses `N0`, `prev_A` and `prev_R` to
##' calculate the starting position
##' @param demographic_pars default NULL uses `demographic_parameters()`)
##' @return A list of all model parameters
##' @export
model_parameters <- function(gas_pars, initial_pars = NULL,
                             demographic_pars = NULL) {
  demographic_pars <- demographic_pars %||% demographic_parameters()
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
##' @return A list of demographic parameters
##' @export
demographic_parameters <- function() {
  list(N0 = 56e6,
       alpha = 7e5,
       omega = 0.0125,
       m = matrix(1, 1, 1)) # mixing matrix
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
