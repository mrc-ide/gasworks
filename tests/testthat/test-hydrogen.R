
test_that("hydrogen_fitted_states", {
  expect_equal(hydrogen_fitted_states(),
               c("daily_pharyngitis_rate", "scarlet_fever_inc", "igas_inc"))
})

test_that("hydrogen_index", {
  p <- example_gas_parameters(1)
  mod <- model$new(p, 0, 10)
  idx <- hydrogen_index(mod$info())

  expect_equal(names(idx), c("run", "state"))
  expect_equal(names(idx$run), hydrogen_fitted_states())
  expect_true(all(hydrogen_fitted_states() %in% names(idx$state)))

  # check can only use on hydrogen model
  mod <- model$new(example_gas_parameters(2), 0, 10)
  expect_error(hydrogen_index(mod$info()))
})

test_that("hydrogen_compare", {
  state <- rbind(daily_pharyngitis_rate = 10:15,
                 scarlet_fever_inc = 100:105,
                 igas_inc = 90:95)
  observed <- list(daily_pharyngitis_rate = 13, scarlet_fever_inc = 103,
                   igas_inc = 93)
  pars <- example_gas_parameters(1)
  y <- hydrogen_compare(state, observed, pars)

  expect_equal(max(y), y[3])
  expect_true(all(y > hydrogen_compare(state * 5, observed, pars)))
})
