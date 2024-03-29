test_that("helium age groups are parametrised correctly", {
  groups <- helium_age_groups()
  expect_equal(16L, groups$n_group)
  expect_equal(length(groups$age_start), groups$n_group)
  expect_equal(length(groups$age_end), groups$n_group)
  expect_true(all(groups$age_end > groups$age_start))
  expect_equivalent(unlist(groups$idx_ukhsa), seq_len(groups$n_group))
})

test_that("helium_index", {
  p <- example_parameters(16)
  mod <- model$new(p, 0, 10)
  idx <- helium_index(mod$info())

  expect_equal(names(idx), c("run", "state"))
  state_nms <- c("scarlet_fever_cases", "igas_inc",
                 "daily_gas_pharyngitis_rate", "etiologic_fraction",
                 "etiologic_fraction_04",
                 "etiologic_fraction_05_14",
                 "etiologic_fraction_15_44",
                 "etiologic_fraction_45_64",
                 "etiologic_fraction_65_74",
                 "etiologic_fraction_75",
                 "daily_gas_pharyngitis_rate_04",
                 "daily_gas_pharyngitis_rate_05_14",
                 "daily_gas_pharyngitis_rate_15_44",
                 "daily_gas_pharyngitis_rate_45_64",
                 "daily_gas_pharyngitis_rate_65_74",
                 "daily_gas_pharyngitis_rate_75",
                 "daily_scarlet_fever_rate_04",
                 "daily_scarlet_fever_rate_05_14",
                 "daily_scarlet_fever_rate_15_44",
                 "daily_scarlet_fever_rate_45_64",
                 "daily_scarlet_fever_rate_65_74",
                 "daily_scarlet_fever_rate_75")
  expect_equal(names(idx$run), state_nms)
  expect_true(all(state_nms %in% names(idx$state)))

  # check can only use on helium model
  mod <- model$new(example_parameters(2), 0, 10)
  expect_error(helium_index(mod$info()))
})

test_that("helium_compare", {
  transform <- helium_create_transform(NULL)
  pars <- transform(example_helium_parameters())
  mod <- model$new(pars, 0, 5, seed = 1L)
  info <- mod$info()

  full_state <- mod$run(7)
  rownames(full_state) <- names(unlist(info$index))
  state <- full_state[helium_index(info)$run, ]
  groups <- helium_age_groups()

  sf_rate <- full_state[grep("daily_sca", rownames(full_state)), 3]
  phi_S <- state[grep("^etio", rownames(state)), 3]
  pharyngitis_rate <- state[grep("^daily_gas", rownames(state)), 3] / phi_S
  names(pharyngitis_rate) <- gsub("gas_", "", names(pharyngitis_rate))

  observed <- as.list(c(state[c("scarlet_fever_cases", "igas_inc"), 3],
                        sf_rate, pharyngitis_rate))

  compare_compiled <- function(state, observed, pars, mod) {
    idx <- helium_index(mod$info())$run
    state_full <- mod$state()
    state_full[idx, ] <- state
    mod$update_state(pars = pars, state = state_full)
    mod$set_data(list(list(mod$time(), observed)))
    mod$compare_data()
  }

  set.seed(1)
  ll <- helium_compare(state, observed, pars)

  expect_equal(length(ll), ncol(state))
  expect_true(all(ll > helium_compare(state * 5, observed, pars)))
  expect_equal(ll, c(-24.813554694153, -33.7484013304405, -24.0115119493748,
                     -25.1792899314225, -24.2923960810031))
  expect_equal(compare_compiled(state, observed, pars, mod), ll)

  # check loglikelihood is maximised at data point in univariate sensitivity
  # analysis
  idx <- grep("^N", rownames(state), invert = TRUE, value = TRUE)
  x <- seq(0.5, 1.5, length.out = 101)

  for (i in idx) {
    tmp <- matrix(state[, 3], nrow = nrow(state), ncol = 101,
                  dimnames = list(rownames(state), NULL))
    tmp[i, ] <- tmp[i, ] * x
    # rebase so that scarlet fever cases sum to same number
    nms <- grep("scarlet_fever_inc_", rownames(tmp), value = TRUE)
    tmp[nms, ] <- round(t(t(tmp[nms, ]) / colSums(tmp[nms, ])) *
                          sum(state[nms, 3]))

    y <- helium_compare(tmp, observed, pars)
    expect_equivalent(max(y), y[x == 1], tolerance = 1 / sqrt(pars$exp_noise))
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
            state["daily_gas_pharyngitis_rate", ], pars$k_gp))
  expect_equal(
    helium_compare(state, observed_agg, pars),
    compare_compiled(state, observed_agg, pars, mod))

  # NA data returns 0
  observed_empty <- replace(observed, seq_along(observed), NA)
  expect_equal(helium_compare(state, observed_empty, pars),
               rep(0, ncol(state)))
  expect_equal(compare_compiled(state, observed_empty, pars, mod),
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

test_that("helium_filter", {
  data <- data.frame(model_week = 1:6,
                     daily_pharyngitis_rate = 10:15,
                     daily_scarlet_fever_rate = 1:6,
                     scarlet_fever_cases = 100:105,
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
  expect_true(all(diff(x$time_start) == 7))

  df <- data.frame(age_start = 0, age_end = 5, N = 1e3, n = 1e2)
  constant_data <- list(etiologic_fraction = df, asymptomatic_carriage = df)
  constant_ll <- create_constant_log_likelihood(constant_data)

  filter <- helium_filter(data, constant_data, 3)
  expect_equal(filter$n_particles, 3)
  expect_equal(filter$inputs()$compare, helium_compare)
  expect_equal(filter$inputs()$index, helium_index)
  expect_equal(filter$inputs()$data, x)
  expect_equal(filter$inputs()$constant_log_likelihood, constant_ll)


  ## check filter runs
  transform <- helium_create_transform(NULL)
  pars <- example_helium_parameters()
  set.seed(1)
  filter$run(transform(pars))

  ## Check compiled filter runs - hard to test values here because the
  ## variance is large with this set of parameters.
  filter_compiled <- helium_filter(data, constant_data, 3, TRUE)
  expect_null(filter_compiled$inputs()$compare)
  filter_compiled$run(transform(pars))

  expect_true(is.function(filter$.__enclos_env__$private$compare))
  expect_null(filter_compiled$.__enclos_env__$private$compare)
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
