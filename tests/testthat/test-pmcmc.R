
test_that("compare function works", {
  pars <- example_gas_parameters()
  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  state <- drop(mod$simulate(7))
  rownames(state) <- names(model_index())
  observed <- list(igas = 25, scarlet_fever = 53, pharyngitis = 50)
  ll <- compare(state, observed, pars)
  expect_equal(ll, c(-2074.18293142027, -2069.99939812045, -2069.33021034521,
                     -2074.03514771718, -2074.01605760047))

  nms <- c("igas_inc", "scarlet_fever_inc", "pharyngitis_inc")
  observed2 <- state[nms, ]
  rownames(observed2) <- gsub("_inc", "", nms)
  observed2["pharyngitis", ] <- calculate_pharyngitis_incidence(state, pars)

  ll_max <- compare(state, as.list(observed2[, 1]), pars)
  expect_true(mean(ll_max) > mean(ll))

  # check can deal with zero model trajectories
  state[nms, ] <- 0
  ll <- compare(state, observed, pars)
  expect_true(all(ll > -Inf))

})
