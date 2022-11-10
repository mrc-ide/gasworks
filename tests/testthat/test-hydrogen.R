
test_that("hydrogen_fitted_states", {
  expect_equal(hydrogen_fitted_states(),
               c("daily_pharyngitis_rate", "daily_scarlet_fever_rate",
                 "scarlet_fever_cases", "igas_inc"))
})

test_that("hydrogen_index", {
  p <- example_parameters(1)
  mod <- model$new(p, 0, 10)
  idx <- hydrogen_index(mod$info())

  expect_equal(names(idx), c("run", "state"))
  expect_equal(names(idx$run),
               c("daily_gas_pharyngitis_rate", "daily_scarlet_fever_rate",
                 "scarlet_fever_cases", "igas_inc"))

  # check can only use on hydrogen model
  mod <- model$new(example_parameters(2), 0, 10)
  expect_error(hydrogen_index(mod$info()))
})

test_that("hydrogen_compare", {
  pars <- example_parameters(1)
  state <- rbind(daily_gas_pharyngitis_rate = 10:15,
                 daily_scarlet_fever_rate = 1:6,
                 scarlet_fever_cases = 100:105,
                 igas_inc = 90:95)
  observed <- list(daily_pharyngitis_rate = 13 / pars$phi_S,
                   daily_scarlet_fever_rate = 4,
                   scarlet_fever_cases = 103,
                   igas_inc = 93)

  y <- hydrogen_compare(state, observed, pars)

  expect_equal(length(y), ncol(state))
  expect_true(all(y > hydrogen_compare(state * 5, observed, pars)))

  # check loglikelihood is maximised at data point in univariate sensitivity
  # analysis
  x <- seq(0.5, 1.5, length.out = 101)
  for (i in seq_along(nrow(state))) {
    tmp <- matrix(state[, 4], nrow = nrow(state), ncol = 101,
                  dimnames = list(rownames(state), NULL))
    tmp[i, ] <- tmp[i, ] * x
    y <- hydrogen_compare(tmp, observed, pars)
    expect_equal(which.max(y), which(x == 1))
  }

  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  idx <- hydrogen_index(mod$info())$run
  state <- drop(mod$simulate(7))[idx, ]
  rownames(state) <- names(idx)
  observed <- list(igas_inc = 100, scarlet_fever_cases = 700,
                   daily_pharyngitis_rate = 1400,
                   daily_scarlet_fever_rate = 0.17)
  ll <- hydrogen_compare(state, observed, pars)
  expect_equal(ll, c(-17.4284917966389, -17.4969689904182, -17.6260807979267,
                     -17.5544631699755, -17.6074746321239))

  # NA data returns 0
  expect_equal(hydrogen_compare(state, replace(observed, 1:4, NA), pars),
               rep(0, ncol(state)))
  expect_equal(hydrogen_compare(state, replace(observed, 1:4, -1), pars),
               rep(-Inf, ncol(state)))

  expect_error(hydrogen_compare(state, unname(observed), pars),
               "missing or misnamed data")
  expect_error(hydrogen_compare(state, observed, example_parameters(2)))
})

test_that("hydrogen_filter", {
  data <- data.frame(model_week = 1:6,
                     daily_pharyngitis_rate = 10:15,
                     daily_scarlet_fever_rate = 1:6,
                     scarlet_fever_cases = 100:105,
                     igas_inc = 90:95)
  x <- hydrogen_prepare_data(data)
  expect_true(all(diff(x$time_start) == 7))

  filter <- hydrogen_filter(data, 3)
  expect_equal(filter$n_particles, 3)
  expect_equal(filter$inputs()$compare, hydrogen_compare)
  expect_equal(filter$inputs()$index, hydrogen_index)
  expect_equal(filter$inputs()$data, x)

  # can run filter
  pars <- example_parameters()
  filter$run(pars)
})

test_that("hydrogen_create_transform", {
  dem_pars <- list(N0 = 10,
                   alpha = 1,
                   omega = 0,
                   r_age = 2,
                   m = matrix(0))
  transform <- hydrogen_create_transform(dem_pars)
  gas_pars <- no_gas_parameters()
  gas_pars$q_F_2014 <- 0.3
  gas_pars$q_F_2020 <- 0.6
  pars <- transform(gas_pars)

  expect_equal(names(pars), unique(names(pars)))
  expect_equal(pars[names(dem_pars)], dem_pars)
})
