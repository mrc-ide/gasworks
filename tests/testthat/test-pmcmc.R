
test_that("compare function works", {
  pars <- example_gas_parameters()
  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  state <- drop(mod$simulate(7))
  rownames(state) <- names(model_index())
  observed <- list(igas_inc = 25, scarlet_fever_inc = 53, pharyngitis_rate = 50)
  ll <- compare(state, observed, pars)
  expect_equal(ll, c(-1045.82877080168, -1043.67265912834, -1043.36126286318,
                     -1045.77521389557, -1045.71150508569))

  nms <- c("igas_inc", "scarlet_fever_inc", "pharyngitis_rate")
  observed2 <- as.list(state[nms, 1])

  ll_max <- compare(state, observed2, pars)
  expect_true(mean(ll_max) > mean(ll))

  # check data names
  names(observed2) <- gsub("_inc", "", names(observed2))
  expect_error(compare(state, observed2, pars), "missing or misnamed data")

  # check can deal with zero model trajectories
  state[nms, ] <- 0
  ll <- compare(state, observed, pars)
  expect_true(all(ll > -Inf))
})
