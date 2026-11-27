test_that("helium age groups are parametrised correctly", {
  groups <- helium_age_groups()
  expect_equal(16L, groups$n_group)
  expect_equal(length(groups$age_start), groups$n_group)
  expect_equal(length(groups$age_end), groups$n_group)
  expect_true(all(groups$age_end > groups$age_start))
  expect_equivalent(unlist(groups$idx_ukhsa), seq_len(groups$n_group))
})

test_that("helium_saved_states", {
  nms <- helium_saved_states()
  expect_equal(nms, unique(nms))

  ## all are states of the model
  sys <- dust2::dust_system_create(model, example_parameters(16))
  expect_true(all(nms %in% names(dust2::dust_unpack_index(sys))))

  ## includes the modelled equivalent of everything we fit to
  fitted <- sub("pharyngitis", "gas_pharyngitis", helium_fitted_states())
  expect_true(all(fitted %in% nms))
})

test_that("helium_compare", {
  transform <- helium_create_transform(NULL)
  pars <- transform(example_helium_parameters())
  sys <- system_at_week(pars)
  full_state <- system_state(sys)
  groups <- names(helium_age_groups()$idx_ukhsa)
  nms <- c("scarlet_fever_cases", "igas_inc", "daily_gas_pharyngitis_rate",
           "etiologic_fraction",
           paste0("etiologic_fraction_", groups),
           paste0("daily_gas_pharyngitis_rate_", groups),
           paste0("daily_scarlet_fever_rate_", groups))
  state <- full_state[nms, ]

  sf_rate <- state[grep("daily_sca", rownames(state)), 3]
  phi_S <- state[grep("^etio", rownames(state)), 3]
  pharyngitis_rate <- state[grep("^daily_gas", rownames(state)), 3] / phi_S
  names(pharyngitis_rate) <- gsub("gas_", "", names(pharyngitis_rate))

  observed <- as.list(c(state[c("scarlet_fever_cases", "igas_inc"), 3],
                        sf_rate, pharyngitis_rate))

  ## The R and compiled versions draw their (tiny) exponential noise from
  ## different random number streams
  tolerance <- 1e-5

  ll <- helium_compare(state, observed, pars)

  expect_equal(length(ll), ncol(state))
  expect_true(all(ll > helium_compare(state * 5, observed, pars)))
  expect_equal(compare_compiled(sys, observed, helium_prepare_data), ll,
               tolerance = tolerance)
  expect_equal(
    compare_compiled(sys, observed, helium_prepare_data, state * 5),
    helium_compare(state * 5, observed, pars), tolerance = tolerance)
  compare_compiled(sys, observed, helium_prepare_data, state)

  # check loglikelihood is maximised at data point in univariate sensitivity
  # analysis
  idx <- grep("^N", rownames(state), invert = TRUE, value = TRUE)
  x <- seq(0.5, 1.5, length.out = 101)

  for (i in idx) {
    tmp <- matrix(state[, 3], nrow = nrow(state), ncol = 101,
                  dimnames = list(rownames(state), NULL))
    tmp[i, ] <- tmp[i, ] * x

    y <- helium_compare(tmp, observed, pars)
    expect_equal(max(y), y[x == 1], tolerance = 1 / sqrt(pars$exp_noise))
  }

  # fits use aggregate data when pharyngitis by age is all NA
  observed_agg <-
    replace(observed, grep("daily_pharyngitis_rate_", names(observed)), NA)
  observed_none <-
    replace(observed, grep("daily_pharyngitis_rate", names(observed)), NA)
  expect_equal(
    helium_compare(state, observed_agg, pars),
    helium_compare(state, observed_none, pars) +
    ll_norm(observed$daily_pharyngitis_rate * state["etiologic_fraction", ],
            state["daily_gas_pharyngitis_rate", ], pars$k_gp),
    tolerance = tolerance)
  expect_equal(
    compare_compiled(sys, observed_agg, helium_prepare_data),
    helium_compare(state, observed_agg, pars), tolerance = tolerance)
  expect_equal(
    compare_compiled(sys, observed_none, helium_prepare_data),
    helium_compare(state, observed_none, pars), tolerance = tolerance)

  # NA data returns 0
  observed_empty <- replace(observed, seq_along(observed), NA)
  expect_equal(helium_compare(state, observed_empty, pars),
               rep(0, ncol(state)))
  expect_equal(compare_compiled(sys, observed_empty, helium_prepare_data),
               rep(0, ncol(state)))
  expect_error(helium_compare(state, unname(observed), pars),
               "missing or misnamed data")

  expect_error(helium_compare(state, observed, example_parameters(1)))
})

test_that("create_constant_log_likelihood", {
  df <- data.frame(age_start = seq(0, 20, 5),
                   age_end = seq(4, 24, 5),
                   N = 1e3)
  data <- list(etiologic_fraction = df,
               asymptomatic_carriage = df)
  pars <- list(b0_phi_S = 5, b1_phi_S = 2, b2_phi_S = 0.2,
               b0_prev_A = 10, b1_prev_A = 3, b2_prev_A = 0.3)

  coefs <- helium_get_spline_coefficients(c("phi_S", "prev_A"), pars)
  phi_S <- mean_age_spline(df$age_start, df$age_end, coefs$phi_S,
                           age_spline_gamma)
  prev_A <- mean_age_spline(df$age_start, df$age_end, coefs$prev_A,
                           age_spline_gamma)
  data$etiologic_fraction$n <- round(df$N * phi_S)
  data$asymptomatic_carriage$n <- round(df$N * prev_A)
  constant_ll <- create_constant_log_likelihood(data)

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

example_helium_data <- function() {
  data.frame(model_week = 1:6,
             daily_pharyngitis_rate = 10:15,
             daily_scarlet_fever_rate = 1:6,
             scarlet_fever_cases = 100:105,
             igas_inc = 90:95,
             daily_pharyngitis_rate_04 = c(0:4, NA),
             daily_pharyngitis_rate_05_14 = c(1:5, NA),
             daily_pharyngitis_rate_15_44 = c(2:6, NA),
             daily_pharyngitis_rate_45_64 = c(3:7, NA),
             daily_pharyngitis_rate_65_74 = c(4:8, NA),
             daily_pharyngitis_rate_75 = c(5:9, NA),
             daily_scarlet_fever_rate_04 = 6:11,
             daily_scarlet_fever_rate_05_14 = 7:12,
             daily_scarlet_fever_rate_15_44 = 8:13,
             daily_scarlet_fever_rate_45_64 = 9:14,
             daily_scarlet_fever_rate_65_74 = 10:15,
             daily_scarlet_fever_rate_75 = 11:16)
}

test_that("helium_prepare_data", {
  data <- example_helium_data()
  x <- helium_prepare_data(data)
  expect_s3_class(x, "dust_filter_data")
  expect_equal(x$time, data$model_week * 7)
  expect_equal(names(x), c("time", paste0("obs_", model_data_names())))

  ## the all-age scarlet fever rate is not fitted in the helium model
  expect_true(all(is.na(x$obs_daily_scarlet_fever_rate)))
  ## the all-age pharyngitis rate is only used where there is none by age
  expect_equal(x$obs_daily_pharyngitis_rate, c(rep(NA, 5), 15))
  expect_equal(x$obs_daily_pharyngitis_rate_75, c(5:9, NA))
  expect_equal(x$obs_igas_inc, data$igas_inc)

  expect_error(helium_prepare_data(data[names(data) != "igas_inc"]),
               "missing or misnamed data: igas_inc")
  expect_error(helium_prepare_data(data[names(data) != "model_week"]),
               "missing or misnamed data: model_week")
})

test_that("helium_filter", {
  data <- example_helium_data()
  transform <- helium_create_transform(NULL)
  pars <- transform(example_helium_parameters())

  filter <- helium_filter(data, 3, seed = 1L)
  expect_s3_class(filter, "dust_likelihood")
  expect_false(filter$deterministic)
  expect_equal(filter$n_particles, 3)
  ll <- dust2::dust_likelihood_run(filter, pars)
  expect_true(is.finite(ll))

  ## same seed gives same result, different seed does not
  filter2 <- helium_filter(data, 3, seed = 1L)
  expect_identical(dust2::dust_likelihood_run(filter2, pars), ll)
  filter3 <- helium_filter(data, 3, seed = 2L)
  expect_false(dust2::dust_likelihood_run(filter3, pars) == ll)

  ## deterministic likelihood is the sum of the comparison over the data
  unfilter <- helium_filter(data, 3, deterministic = TRUE)
  expect_true(unfilter$deterministic)
  ll_det <- dust2::dust_likelihood_run(unfilter, pars)
  expect_identical(dust2::dust_likelihood_run(unfilter, pars), ll_det)

  x <- helium_prepare_data(data)
  sys <- dust2::dust_system_create(model, pars, deterministic = TRUE)
  dust2::dust_system_set_state_initial(sys)
  expected <- 0
  for (i in seq_len(nrow(x))) {
    dust2::dust_system_run_to_time(sys, x$time[[i]])
    expected <- expected + dust2::dust_system_compare_data(
      sys, as.list(x[i, names(x) != "time"]))
  }
  expect_equal(ll_det, expected)
})

test_that("helium_constant_likelihood", {
  df <- data.frame(age_start = c(0, 5), age_end = c(4, 9), N = 1e3, n = 1e2)
  constant_data <- list(etiologic_fraction = df, asymptomatic_carriage = df)
  pars <- example_helium_parameters()
  nms <- grep("^b[0-9]_|^q_F_", names(pars), value = TRUE)
  packer <- helium_packer(nms, NULL,
                          fixed = pars[setdiff(names(pars), c(nms, "q_F"))])

  likelihood <- helium_constant_likelihood(constant_data, packer)
  expect_s3_class(likelihood, "monty_model")
  expect_equal(likelihood$parameters, nms)

  x <- packer$pack(pars)
  expect_equal(monty::monty_model_density(likelihood, x),
               create_constant_log_likelihood(constant_data)(pars))

  ## can be combined with the filter and a prior to give a posterior
  filter <- helium_filter(example_helium_data(), 3, deterministic = TRUE)
  ll <- dust2::dust_likelihood_monty(filter, packer)
  posterior <- ll + likelihood
  expect_equal(
    monty::monty_model_density(posterior, x),
    dust2::dust_likelihood_run(filter, packer$unpack(x)) +
      monty::monty_model_density(likelihood, x))
})

test_that("helium_create_transform", {
  n_group <- helium_age_groups()$n_group
  dem_pars <- list(N0 = seq_len(n_group),
                   alpha = 1,
                   omega = seq(0.1, 0.3, length.out = n_group),
                   r_age = 2,
                   m = matrix(0, n_group, n_group))
  transform <- helium_create_transform(dem_pars)
  gas_pars <- example_helium_parameters()
  pars <- transform(gas_pars)

  expect_equal(names(pars), unique(names(pars)))
  expect_equal(pars[names(dem_pars)], dem_pars)
})

test_that("helium_packer", {
  n_group <- helium_age_groups()$n_group
  dem_pars <- list(N0 = seq_len(n_group) * 100,
                   alpha = 1,
                   omega = seq(0.1, 0.3, length.out = n_group),
                   r_age = 2,
                   m = matrix(0, n_group, n_group))
  transform <- helium_create_transform(dem_pars)
  gas_pars <- example_helium_parameters()
  nms <- c("beta", "sigma", "b0_phi_S", "k_gp", "k_hpr", "delta_A")
  fixed <- gas_pars[setdiff(names(gas_pars), c(nms, "q_F"))]

  packer <- helium_packer(nms, dem_pars, fixed)
  expect_equal(packer$names(), nms)

  x <- c(1, 0.5, 6, 0.7, 0.2, 20)
  pars <- packer$unpack(x)
  expected <- transform(c(setNames(as.list(x), nms), fixed))
  expect_mapequal(pars, expected)

  ## fitted values reach the model unchanged
  expect_equal(unlist(pars[nms]), setNames(x, nms))
  expect_equal(packer$pack(pars), x)

  ## and the result can be used to run the model
  y <- model_simulate(pars, 7)
  expect_true(all(y >= 0))
})

test_that("packer errors if the transform changes a fitted parameter", {
  transform <- function(pars) {
    pars$a <- pars$a * 2
    pars$c <- 1
    pars
  }
  packer <- create_packer(c("a", "b"), transform)
  expect_error(packer$unpack(c(1, 2)),
               "transform must not change the value of: a")

  packer <- create_packer(c("a", "b"), function(pars) c(pars, list(c = 1)),
                          fixed = list(d = 4))
  expect_equal(packer$unpack(c(1, 2)), list(a = 1, b = 2, d = 4, c = 1))
})
