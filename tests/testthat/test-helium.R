test_that("helium age groups are parametrised correctly", {
  groups <- helium_age_groups()
  expect_equal(16L, groups$n_group)
  expect_equal(length(groups$age_start), groups$n_group)
  expect_equal(length(groups$age_end), groups$n_group)
  expect_true(all(groups$age_end > groups$age_start))
  expect_equivalent(unlist(groups$idx_ukhsa), seq_len(groups$n_group))
})

test_that("helium_index", {
  p <- example_gas_parameters(16)
  mod <- model$new(p, 0, 10)
  idx <- helium_index(mod$info())

  expect_equal(names(idx), c("run", "state"))
  state_nms <- c("scarlet_fever_inc", "igas_inc", "daily_pharyngitis_rate",
                 paste0("N", seq_len(16)),
                 "pharyngitis_prop_04", "pharyngitis_prop_05_14",
                 "pharyngitis_prop_15_44", "pharyngitis_prop_45_64",
                 "pharyngitis_prop_65_74", "pharyngitis_prop_75",
                 "scarlet_fever_prop_04", "scarlet_fever_prop_05_14",
                 "scarlet_fever_prop_15_44", "scarlet_fever_prop_45_64",
                 "scarlet_fever_prop_65_74", "scarlet_fever_prop_75")
  expect_equal(names(idx$run), state_nms)
  expect_true(all(state_nms %in% names(idx$state)))

  # check can only use on helium model
  mod <- model$new(example_gas_parameters(2), 0, 10)
  expect_error(helium_index(mod$info()))
})

test_that("helium_compare", {
  pars <- example_gas_parameters(16)
  mod <- model$new(pars, 0, 5, seed = 1L)
  info <- mod$info()

  full_state <- mod$run(7)
  rownames(full_state) <- names(unlist(info$index))
  state <- full_state[helium_index(info)$run, ]
  observed <- as.list(full_state[helium_fitted_states(), 3])

  ll <- helium_compare(state, observed, pars)

  expect_equal(length(ll), ncol(state))
  expect_true(all(ll > helium_compare(state * 5, observed, pars)))
  expect_equal(ll, c(-103.298262175227, -84.421679833637, -68.4963345542519,
                     -107.834192364211, -83.0004051778352))

  # check matches hydrogen_compare when age-based fitting is not used
  idx <- hydrogen_fitted_states()
  observed_na <- replace(observed, setdiff(names(observed), idx), NA)
  expect_equal(hydrogen_compare(state[idx, ], observed[idx],
                                replace(pars, "n_group", 1)),
               helium_compare(state, observed_na, pars))

  # check loglikelihood is maximised at data point in univariate sensitivity
  # analysis
  idx <- grep("^N", rownames(state), invert = TRUE)
  par(mfrow = c(3, 5), bty = "n", mar = c(3, 3, 1, 1))
  for (i in idx) {
    tmp <- matrix(state[, 3], nrow = nrow(state), ncol = 101,
                  dimnames = list(rownames(state), NULL))
    tmp[i, ] <- tmp[i, ] * seq(0.5, 1.5, length.out = 101)
    y <- helium_compare(tmp, observed, pars)
    expect_equal(which.max(y), ceiling(ncol(tmp) / 2))
  }


  # NA data returns 0
  expect_equal(helium_compare(state, replace(observed, seq_along(observed), NA),
                              pars), rep(0, ncol(state)))

  expect_error(helium_compare(state, unname(observed), pars),
               "missing or misnamed data")
  expect_error(helium_compare(state, observed, example_gas_parameters(1)))
})

test_that("create_constant_log_likelihood", {
  df <- data.frame(age_start = seq(0, 20, 5),
                   age_end = seq(4, 24, 5),
                   N = 1e3)
  data <- list(etiologic_fraction = df,
               asymptomatic_carriage = df)
  constant_ll <- create_constant_log_likelihood(data)
  pars <- list(b0_phi_S = 5, b1_phi_S = 2, b2_phi_S = 0.2,
               b0_prev_A = 10, b1_prev_A = 3, b2_prev_A = 0.3)

  coefs <- helium_get_spline_coefficients(c("phi_S", "prev_A"), pars)
  phi_S <- mean_age_spline(df$age_start, df$age_end, coefs$phi_S,
                           age_spline_gamma)
  prev_A <- mean_age_spline(df$age_start, df$age_end, coefs$prev_A,
                           age_spline_gamma)
  data$etiologic_fraction$n <- round(df$N * phi_S)
  data$asymptomatic_carriage$n <- round(df$N * prev_A)

  # check likelihood is maximised at true parameters
  x <- seq(0.6, 1.4, length.out = 101)
  for (i in seq_along(pars)) {
    y <- sapply(x, function(x) {
      pars[[i]] <- pars[[i]] * x
      constant_ll(pars)
    })
    expect_equal(which.max(y), which(x == 1))
  }
})

test_that("helium_filter", {
  data <- data.frame(model_week = 1:6,
                     daily_pharyngitis_rate = 10:15,
                     scarlet_fever_inc = 100:105,
                     igas_inc = 90:95,
                     daily_pharyngitis_rate_04 = 0:5,
                     daily_pharyngitis_rate_05_14 = 1:6,
                     daily_pharyngitis_rate_15_44 = 2:7,
                     daily_pharyngitis_rate_45_64 = 3:8,
                     daily_pharyngitis_rate_65_74 = 4:9,
                     daily_pharyngitis_rate_75 = 5:10,
                     daily_scarlet_fever_rate_04 = 6:11,
                     daily_scarlet_fever_rate_05_14 = 7:12,
                     daily_scarlet_fever_rate_15_44 = 8:13,
                     daily_scarlet_fever_rate_45_64 = 9:14,
                     daily_scarlet_fever_rate_65_74 = 10:15,
                     daily_scarlet_fever_rate_75 = 11:16)
  x <- helium_prepare_data(data)
  expect_true(all(diff(x$step_start) == 7))

  df <- data.frame(age_start = 0, age_end = 5, N = 1e3, n = 1e2)
  constant_data <- list(etiologic_fraction = df, asymptomatic_carriage = df)
  constant_ll <- create_constant_log_likelihood(constant_data)

  filter <- helium_filter(data, constant_data, 3)
  expect_equal(filter$n_particles, 3)
  expect_equal(filter$inputs()$compare, helium_compare)
  expect_equal(filter$inputs()$index, helium_index)
  expect_equal(filter$inputs()$data, x)
  expect_equal(filter$inputs()$constant_log_likelihood, constant_ll)

  ## check filter runs and we can reclaim the same value
  pars <- example_gas_parameters(16)
  spline_pars <- list(b0_phi_S = 5, b1_phi_S = 2, b2_phi_S = 0.2,
                      b0_prev_A = 10, b1_prev_A = 3, b2_prev_A = 0.3)
  pars <- c(pars, spline_pars)
  set.seed(1)
  expect_equal(filter$run(pars), -21506639)
})
