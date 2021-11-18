
test_that("compare function works", {
  pars <- transform(example_gas_parameters())
  set.seed(1L)
  mod <- model$new(pars, 0, 5, seed = 1L)
  state <- mod$simulate(7)
  rownames(state) <- names(model_index())
  observed <- list(igas = 25, scarlet_fever = 53, pharyngitis = 50)
  ll <- compare(state, observed, pars)
  expect_equal(ll,
               c(-22.3664096220503, -22.3187335697625, -270.393415382723,
                 -22.2261152715995, -21.9178967447101))

  observed <- state[c("igas_inc", "scarlet_fever_inc", "pharyngitis_inc"), , ]
  rownames(observed) <- gsub("_inc", "", rownames(observed))
  observed["pharyngitis", ] <-
    calculate_pharyngitis_incidence(state, model_index(), pars)

  ll_max <- compare(state, as.list(observed[, 1]), pars)
  expect_true(mean(ll_max) > mean(ll))
})
