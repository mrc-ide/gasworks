
pars <- list(U0 = 5e7,
             A0 = 6e6,
             beta = 0.01,
             sigma = 0.1,
             t_s = 100,
             p_S = 0.6,
             p_R = 0.3,
             p_I = 0.0001,
             p_F = 0.001,
             delta_A = 30,
             delta_E = 2,
             delta_I = 14,
             delta_S = 2.3,
             delta_F = 7,
             delta_R = 365 * 5,
             alpha = 7e5 / 365,
             omega = 0.01)

test_that("model runs", {
  mod <- model$new(pars, 0, 5)
  y <- mod$simulate(7)
  idx <- unlist(mod$info()$index)

  # check that states sum to N
  expect_equal(colSums(y[seq(2, 9), , ]), y[idx["N"], , ])
  expect_true(all(y > 0))
})

test_that("there are no infections when beta is 0", {
  pars$beta <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there are no infections when A0 = 0", {
  pars$A0 <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R == 0))

  expect_true(all(unlist(y) >= 0))
})


test_that("there are no symptomatic infections when p_S = 0", {
  pars$p_S <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there are no asymptomatic infections when p_S = 1", {
  pars$p_S <- 1
  pars$E0 <- pars$A0
  pars$A0 <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A == 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there is no iGAS when p_I = 0", {
  pars$p_I <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there is no pharyngitis when p_I = 1", {
  pars$p_I <- 1
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there is no scarlet fever when p_F = 0", {
  pars$p_F <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 > 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there is no pharyngitis when p_F = 1", {
  pars$p_F <- 1
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S1 > 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_true(all(unlist(y) >= 0))
})

test_that("there is no immunity when p_S = 0 and p_R = 0", {
  pars$p_S <- 0
  pars$p_R <- 0
  mod <- model$new(pars, 0, 5)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E == 0))
  expect_true(all(y$S1 == 0))
  expect_true(all(y$S2 == 0))
  expect_true(all(y$I == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R == 0))

  expect_true(all(unlist(y) >= 0))
})


test_that("the foi is calculated correctly", {

})


test_that("incidence time series output correctly", {

})

