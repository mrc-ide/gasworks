##' @name check_gas_parameters
##' @title Check model parameters are within permitted ranges
##' @description Check model parameters are within permitted ranges
##' @param pars a named list of model parameters
##' @return NULL
##' @export

check_gas_parameters <- function(pars) {
  with(pars, {
    assert_scalar_positive(beta)
    assert_scalar_nonnegative(sigma)
    assert_scalar_positive(t_s)
    assert_scalar_unit_interval(p_S)
    assert_scalar_unit_interval(p_R)
    assert_scalar_unit_interval(p_I)
    assert_scalar_unit_interval(p_F)
    assert_scalar_positive(delta_A)
    assert_scalar_positive(delta_E)
    assert_scalar_positive(delta_I)
    assert_scalar_positive(delta_S)
    assert_scalar_positive(delta_F)
    assert_scalar_positive(delta_R)
    assert_scalar_nonnegative_integer(alpha)
    assert_scalar_unit_interval(omega)
    assert_scalar_nonnegative_integer(U0)
    assert_scalar_nonnegative_integer(A0)
    assert_scalar_nonnegative_integer(E0)
    assert_scalar_nonnegative_integer(I0)
    assert_scalar_nonnegative_integer(S0)
    assert_scalar_nonnegative_integer(F0)
    assert_scalar_nonnegative_integer(R0)
  })
}


##' @name example_gas_parameters
##' @title Example of fitted gas parameters for use in testing
##' @description Example of fitted gas parameters for use in testing
##' @return A list of named model parameters
example_gas_parameters <- function() {
  list(prev_A = 0.1,
       prev_R = 0.5,
       beta = 0.01,
       sigma = 0.1,
       t_s = 100,
       p_S = 0.6,
       p_R = 0.3,
       p_I = 0.0001,
       p_F = 0.001,
       delta_A = 30,
       delta_R = 365 * 5)
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
       omega = 0.0125)
}

##' @name initial_parameters
##' @title Initial conditions for the model
##' @description Initial conditions for the model
##' @param pars A parameter list containing `N0`, `prev_A` and `prev_R` elements.
##' @return A named list of initial model states
##' @export
initial_parameters <- function(pars) {
  assert_scalar_unit_interval(pars$prev_A)
  assert_scalar_unit_interval(pars$prev_R)
  assert_scalar_nonnegative_integer(pars$N0)

  # set initial asymptomatic prevalence and prevalence of immunity
  A0 <- round(pars$N0 * pars$prev_A)
  R0 <- round((pars$N0 - A0) * pars$prev_R)

  # set initial uninfecteds
  U0 <- pars$N0 - A0 - R0

  list(U0 = U0, A0 = A0, E0 = 0, I0 = 0, S10 = 0, S20 = 0, F0 = 0, R0 = R0)
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

  # add fixed model parameters (i.e. not fitted)
  pars$delta_E <- 2 # mean days in incubation period
  pars$delta_I <- 14 # mean days iGAS
  pars$delta_S <- 2.3 # mean days with pharyngitis symptoms (x 2)
  pars$delta_F <- 7 # mean days with scarlet fever
  pars$dt <- 1 / 7 # fit to weekly data

  # add initial conditions and demographic parameters
  model_parameters(pars)
}
