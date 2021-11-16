pars <- transform(example_gas_parameters())

test_that("model runs", {
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  idx <- names(mod$info()$index) %in% model_compartments()

  # check that states sum to N
  expect_equal(colSums(y[idx, , ]), y[idx["N"], , ])
  expect_true(all(y >= 0))

  tmp <- c(1, 25700346, 5359260, 12089, 0, 2345, 229, 0, 24929013, 56003282,
           1, 25700753, 5359326, 11754, 0, 2401, 252, 0, 24927187, 56001673,
           1, 25699514, 5358338, 11795, 1, 2354, 247, 0, 24928223, 56000472,
           1, 25700244, 5359021, 11898, 0, 2489, 257, 0, 24926944, 56000853,
           1, 25700111, 5358832, 11903, 2, 2444, 245, 0, 24927715, 56001252)
  expect_equal(y, array(tmp, dim = c(10L, 5L, 1L)))
})

test_that("there are no infections when beta is 0", {
  pars$beta <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there are no infections when A0 = 0", {
  pars$A0 <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})


test_that("there are no symptomatic infections when p_S = 0", {
  pars$p_S <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there are no asymptomatic infections when p_S = 1", {
  pars$p_S <- 1
  pars$E0 <- pars$A0
  pars$A0 <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A == 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there is no iGAS when p_I = 0", {
  pars$p_I <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(all(y$I == 0))
  expect_true(any(y$F > 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there is no pharyngitis when p_I = 1", {
  pars$p_I <- 1
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there is no scarlet fever when p_F = 0", {
  pars$p_F <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(any(y$I > 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there is no pharyngitis when p_F = 1", {
  pars$p_F <- 1
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 == 0))
  expect_true(any(y$I > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("there is no immunity when p_S = 0 and p_R = 0", {
  pars$p_S <- 0
  pars$p_R <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})


test_that("the foi is calculated correctly", {

})


test_that("incidence time series output correctly", {

})
