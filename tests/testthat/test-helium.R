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
