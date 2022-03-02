
test_that("model_day functions work", {

  expect_equal(model_day("2014-01-05"), 7)
  expect_equal(model_day("2013-12-30"), 1)
  expect_equal(model_day(as.Date("2013-12-28")), -1)

  expect_equal(model_date(7), as.Date("2014-01-05"))
  expect_equal(model_date(1), as.Date("2013-12-30"))
  expect_equal(model_date(-1), as.Date("2013-12-28"))
})

test_that("model_week functions work", {

  expect_equal(model_week("2014-01-01"), 1)
  expect_equal(model_week("2014-01-05"), 1)
  expect_equal(model_week("2014-01-06"), 2)
  expect_equal(model_week(as.Date("2013-12-29")), 0)

  expect_equal(model_week_date(1), as.Date("2014-01-05"))
  expect_equal(model_week_date(2), as.Date("2014-01-12"))
  expect_equal(model_week_date(0), as.Date("2013-12-29"))

  # Check for 2015 which has 53 isoweeks
  date <- "2015-12-31"
  expect_equal(lubridate::isoweek(date) + 52, model_week(date))
})


test_that("ll_nbinom", {
  f <- function(model) {
    set.seed(1)
    dnbinom(10, 10, mu = model + rexp(length(model), rate = 1e6), log = TRUE)
  }

  set.seed(1)
  expect_equal(ll_nbinom(10, 10, 0.1, 1e6), f(10))

  x <- 1:10
  set.seed(1)
  expect_equal(ll_nbinom(10, x, 0.1, 1e6), f(x))

  x <- rep(10, 10)
  expect_lt(diff(range(ll_nbinom(10, x, 0.1, 1e6))),
            diff(range(ll_nbinom(10, x, 0.1, 1e2))))
})

test_that("ll_norm", {
  f <- function(model) {
    set.seed(1)
    dnorm(10, model, sqrt(model * 1.1) + rexp(length(model), 1e6),
          log = TRUE)
  }

  set.seed(1)
  expect_equal(ll_norm(10, 10, 0.1, 1e6), f(10))

  x <- 1:10
  set.seed(1)
  expect_equal(ll_norm(10, x, 0.1, 1e6), f(x))

  x <- rep(10, 10)
  expect_lt(diff(range(ll_norm(10, x, 0.1, 1e6))),
            diff(range(ll_norm(10, x, 0.1, 1e2))))
})

test_that("ll_nbinom returns a vector of zeros if data missing", {
  expect_equal(ll_nbinom(NA, 10, 0.1, 1e6), 0)
  expect_equal(ll_nbinom(NA, rep(10, 5), 0.1, 1e6), rep(0, 5))
})

test_that("ll_norm returns a vector of zeros if data missing", {
  expect_equal(ll_norm(NA, 10, 0.1, 1e6), 0)
  expect_equal(ll_norm(NA, rep(10, 5), 0.1, 1e6), rep(0, 5))
})

test_that("ll_multinom", {
  data <- c(1, 5, 10)
  f <- function(state) {
    dmultinom(data, prob = state + 1e-6, log = TRUE)
  }
  ## rows: observations (eg by age)
  ## cols: particles
  state <- matrix(c(1, 2,
                    5, 6,
                    10, 11), byrow = TRUE, nrow = 3)
  expect_equal(apply(state, 2, f), ll_multinom(data, state, noise = 1e-6))

  ## check that zero probabilities are treated as equal
  zero_state <- array(0, dim(state))
  expect_true(all(diff(ll_multinom(data, zero_state, noise = 1e-6)) == 0))

  ## check can deal with missing data
  expect_equal(ll_multinom(c(NA, 1, 1), state, noise = 1e-6), rep(0, 2))
})

test_that("age_spline_polynomial", {
  b0 <- -0.4
  b1 <- 0.4
  b2 <- -0.2
  pars <- c(b0, b1, b2)
  age <- seq(0, 10)
  x <- log(age + 1)

  tmp <- age_spline_polynomial(age, pars)
  expect_equal(tmp, exp(b0 + b1 * x + b2 * x ^ 2))
  expect_equal(tmp[1], exp(b0))

  ## check can only input three coefs
  expect_error(age_spline_polynomial(age, c(b0, b1)))
  expect_error(age_spline_polynomial(age, c(pars, 1)))

  ## check cannot use negative ages
  expect_error(age_spline_polynomial(-1, pars))
})

test_that("age_spline_lognormal", {

  mode <- 10
  sd <- 1
  peak <- 0.2

  pars <- c(mode, sd, peak)
  age <- seq_len(10)

  f <- function(x, mode, sd, peak) {
    v <- sd ^ 2
    mu <- log(mode) + v
    ## algebraic derivation of peak * p(x) / p(mode) where p is lognormal pdf
    peak * exp(mu - v / 2 - (log(x) - mu) ^ 2 / (2 * v)) / x
  }

  expect_equal(age_spline_lognormal(age, pars),  f(age, mode, sd, peak))
  expect_equal(age_spline_lognormal(0, pars), 0)

  ## check can only input three coefs
  expect_error(age_spline_lognormal(age, pars[-1]))
  expect_error(age_spline_lognormal(age, c(pars, 1)))

  ## check coefs are valid
  expect_error(age_spline_lognormal(age, c(-1, sd, peak)))
  expect_error(age_spline_lognormal(age, c(mode, 0, peak)))
  expect_error(age_spline_lognormal(age, c(mode, sd, 2)))

  ## check cannot use negative ages
  expect_error(age_spline_lognormal(-1, pars))
})
