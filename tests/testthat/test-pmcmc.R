
test_that("compare function works", {
  pars <- example_gas_parameters()
  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  idx <- index(mod$info())$run
  state <- drop(mod$simulate(7))[idx, ]
  rownames(state) <- names(idx)
  observed <- list(igas_inc = 25, scarlet_fever_inc = 53,
                   pharyngitis_scarlet_fever_rate = 50)
  ll <- compare(state, observed, pars)
  expect_equal(ll, c(-2446.69560810172, -2445.72002766639, -2445.66270072129,
                     -2446.24777780703, -2449.83251103714))

  observed2 <- as.list(state[, 1])

  ll_max <- compare(state, observed2, pars)
  expect_true(mean(ll_max) > mean(ll))

  # check data names
  names(observed2) <- gsub("_inc", "", names(observed2))
  expect_error(compare(state, observed2, pars), "missing or misnamed data")

  # check can deal with zero model trajectories
  state[, ] <- 0
  ll <- compare(state, observed, pars)
  expect_true(all(ll > -Inf))
})
