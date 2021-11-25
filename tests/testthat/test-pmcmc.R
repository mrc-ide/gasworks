
test_that("compare function works", {
  pars <- example_gas_parameters()
  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  idx <- index(mod$info())$run
  state <- drop(mod$simulate(7))[idx, ]
  rownames(state) <- names(idx)
  observed <- list(igas_inc = 25, scarlet_fever_inc = 53, pharyngitis_rate = 50)
  ll <- compare(state, observed, pars)
  expect_equal(ll, c(-1045.82878036054, -1043.67263594512, -1043.36126271083,
                     -1045.77522055165, -1045.71156016202))

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
