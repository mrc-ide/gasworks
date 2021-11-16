
test_that("demographic_parameters works as expected", {
  pars <- demographic_parameters()
  expect_equal(names(pars), c("N0", "alpha", "omega"))
  expect_equal(pars$N0, 56000000)
  expect_equal(pars$alpha / pars$N0, pars$omega)
})

test_that("initial_parameters works as expected", {
  pars <- c(example_gas_parameters(), demographic_parameters())
  init_pars <- initial_parameters(pars)
  expect_equal(names(init_pars),
               sprintf("%s0", c("U", "A", "E", "I", "S", "F", "R")))
  expect_equal(pars$N0, sum(unlist(init_pars)))
  expect_equal(init_pars$A0, pars$N0 * pars$prev_A)
  expect_equal(init_pars$R0, (pars$N0 - init_pars$A0) * pars$prev_R)
  expect_equal(init_pars$E0, 0)
  expect_equal(init_pars$I0, 0)
  expect_equal(init_pars$S0, 0)
  expect_equal(init_pars$F0, 0)
})

test_that("model_parameters works as expected", {
  gas_pars <- transform(example_gas_parameters())
})


