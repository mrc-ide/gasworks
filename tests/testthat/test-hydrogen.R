
test_that("hydrogen_fitted_states", {
  expect_equal(hydrogen_fitted_states(),
               c("daily_pharyngitis_rate", "daily_scarlet_fever_rate",
                 "scarlet_fever_cases", "igas_inc"))
})

test_that("hydrogen_saved_states", {
  nms <- hydrogen_saved_states()
  expect_equal(nms, unique(nms))

  ## all are states of the model
  sys <- dust2::dust_system_create(model, example_parameters(1))
  expect_true(all(nms %in% names(dust2::dust_unpack_index(sys))))

  ## includes the modelled equivalent of everything we fit to
  fitted <- sub("pharyngitis", "gas_pharyngitis", hydrogen_fitted_states())
  expect_true(all(fitted %in% nms))
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

  ## The R implementation agrees with that compiled into the model, up to
  ## the (tiny) exponential noise, which uses different random number streams
  sys <- system_at_week(pars)
  state <- system_state(sys)[rownames(state), ]
  observed <- list(igas_inc = 100, scarlet_fever_cases = 700,
                   daily_pharyngitis_rate = 1400,
                   daily_scarlet_fever_rate = 0.17)
  ll <- hydrogen_compare(state, observed, pars)
  expect_true(all(is.finite(ll)))
  expect_equal(compare_compiled(sys, observed, hydrogen_prepare_data), ll,
               tolerance = 1e-5)

  observed_empty <- replace(observed, 1:4, NA)
  expect_equal(compare_compiled(sys, observed_empty, hydrogen_prepare_data),
               rep(0, ncol(state)))
  observed_some <- replace(observed, c(1, 3), NA)
  expect_equal(compare_compiled(sys, observed_some, hydrogen_prepare_data),
               hydrogen_compare(state, observed_some, pars),
               tolerance = 1e-5)

  # NA data returns 0
  expect_equal(hydrogen_compare(state, replace(observed, 1:4, NA), pars),
               rep(0, ncol(state)))
  expect_equal(hydrogen_compare(state, replace(observed, 1:4, -1), pars),
               rep(-Inf, ncol(state)))

  expect_error(hydrogen_compare(state, unname(observed), pars),
               "missing or misnamed data")
  expect_error(hydrogen_compare(state, observed, example_parameters(2)))
})

test_that("hydrogen_prepare_data", {
  data <- data.frame(model_week = 1:6,
                     daily_pharyngitis_rate = 10:15,
                     daily_scarlet_fever_rate = 1:6,
                     scarlet_fever_cases = 100:105,
                     igas_inc = 90:95,
                     daily_scarlet_fever_rate_04 = 6:11)
  x <- hydrogen_prepare_data(data)
  expect_s3_class(x, "dust_filter_data")
  expect_equal(x$time, data$model_week * 7)
  expect_equal(names(x), c("time", paste0("obs_", model_data_names())))
  for (nm in hydrogen_fitted_states()) {
    expect_equal(x[[paste0("obs_", nm)]], data[[nm]])
  }
  ## rates by age are not fitted in the hydrogen model
  expect_true(all(is.na(x[grep("_[0-9]+$", names(x))])))

  expect_error(hydrogen_prepare_data(data[-2]),
               "missing or misnamed data: daily_pharyngitis_rate")
})

test_that("hydrogen_filter", {
  data <- data.frame(model_week = 1:6,
                     daily_pharyngitis_rate = 10:15,
                     daily_scarlet_fever_rate = 1:6,
                     scarlet_fever_cases = 100:105,
                     igas_inc = 90:95)
  pars <- example_parameters()

  filter <- hydrogen_filter(data, 3, seed = 1L)
  expect_s3_class(filter, "dust_likelihood")
  expect_false(filter$deterministic)
  expect_equal(filter$n_particles, 3)
  ll <- dust2::dust_likelihood_run(filter, pars)
  expect_true(is.finite(ll))
  expect_identical(
    dust2::dust_likelihood_run(hydrogen_filter(data, 3, seed = 1L), pars), ll)

  unfilter <- hydrogen_filter(data, 3, deterministic = TRUE)
  expect_true(unfilter$deterministic)
  ll_det <- dust2::dust_likelihood_run(unfilter, pars)
  expect_identical(dust2::dust_likelihood_run(unfilter, pars), ll_det)

  ## saved trajectories line up with the data
  dust2::dust_likelihood_run(unfilter, pars, save_trajectories = TRUE)
  h <- dust2::dust_likelihood_last_trajectories(unfilter)
  y <- model_simulate(pars, data$model_week * 7, deterministic = TRUE)
  expect_equal(h, unname(y[, 1, ]))
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

test_that("hydrogen_packer", {
  dem_pars <- list(N0 = 1000,
                   alpha = 1,
                   omega = 0,
                   r_age = 2,
                   m = matrix(0))
  gas_pars <- example_gas_parameters()
  gas_pars$q_F <- NULL
  nms <- c("beta", "q_F_2014", "q_F_2020", "k_gp", "delta_R")
  fixed <- gas_pars[setdiff(names(gas_pars), nms)]

  packer <- hydrogen_packer(nms, dem_pars, fixed)
  expect_equal(packer$names(), nms)

  x <- c(1, 0.3, 0.6, 0.7, 100)
  pars <- packer$unpack(x)
  transform <- hydrogen_create_transform(dem_pars)
  expect_mapequal(pars, transform(c(setNames(as.list(x), nms), fixed)))
  expect_equal(unlist(pars[nms]), setNames(x, nms))
  expect_equal(range(pars$q_F), c(0.3, 0.6))

  ## can be used to create a likelihood for fitting
  data <- data.frame(model_week = 1:3,
                     daily_pharyngitis_rate = 10:12,
                     daily_scarlet_fever_rate = 1:3,
                     scarlet_fever_cases = 100:102,
                     igas_inc = 90:92)
  filter <- hydrogen_filter(data, 3, deterministic = TRUE)
  likelihood <- dust2::dust_likelihood_monty(
    filter, packer, save_trajectories = hydrogen_saved_states())
  expect_equal(monty::monty_model_density(likelihood, x),
               dust2::dust_likelihood_run(filter, pars))
})
