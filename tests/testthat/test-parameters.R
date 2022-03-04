
test_that("demographic_parameters works as expected", {
  pars <- demographic_parameters()
  expect_equal(names(pars), c("N0", "alpha", "omega", "r_age", "m"))
  expect_equal(pars$N0, 56000000)
  expect_equal(pars$alpha / pars$N0, pars$omega)
  expect_equivalent(pars$m, 1)
  expect_equal(pars$r_age, 0)
})

test_that("initial_parameters works as expected", {
  pars <- c(example_parameters(), demographic_parameters())
  init_pars <- initial_parameters(pars)
  expect_equal(names(init_pars),
               paste0(c("A", "E", "S", "P", "F", "R", "U"), 0))
  expect_equivalent(pars$N0, sum(unlist(init_pars)))
  expect_equivalent(init_pars$A0, pars$N0 * pars$prev_A)
  expect_equivalent(init_pars$R0, (pars$N0 - init_pars$A0) * pars$prev_R)
  expect_equivalent(init_pars$E0, rep(0, pars$k_E))
  expect_equivalent(init_pars$S0, rep(0, pars$k_S))
  expect_equivalent(init_pars$P0, rep(0, pars$k_P))
  expect_equivalent(init_pars$F0, rep(0, pars$k_F))
})
