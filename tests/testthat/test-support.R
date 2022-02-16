
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
