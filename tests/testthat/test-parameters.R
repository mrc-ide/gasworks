
test_that("demographic_parameters works as expected", {
  pars <- demographic_parameters()
  expect_equal(names(pars), c("N0", "alpha", "omega", "m"))
  expect_equal(pars$N0, 56000000)
  expect_equal(pars$alpha / pars$N0, pars$omega)
  expect_equivalent(pars$m, 1)
})

test_that("initial_parameters works as expected", {
  pars <- c(example_gas_parameters(), demographic_parameters())
  init_pars <- initial_parameters(pars)
  expect_equal(names(init_pars),
               sprintf("%s0", model_compartments()))
  expect_equal(pars$N0, sum(unlist(init_pars)))
  expect_equal(init_pars$A0, pars$N0 * pars$prev_A)
  expect_equal(init_pars$R0, (pars$N0 - init_pars$A0) * pars$prev_R)
  expect_equal(init_pars$E0, 0)
  expect_equal(init_pars$S0, 0)
  expect_equal(init_pars$P0, 0)
  expect_equal(init_pars$F0, 0)
})

test_that("model_parameters works as expected", {
  expect_error(transform(example_gas_parameters()),
                         "Parameters have already been transformed")
})
