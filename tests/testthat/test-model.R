

test_that("model runs", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(mod$info()$index)

  # check that states sum to N
  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))


  tmp <- c(1, 4071984, 589563, 573, 140, 7e+05, 699998, 1.94454194361955,
           4211.16413531557, 1.02321424917092, 22515830, 5879765, 1031524,
           130, 812108, 418016, 515, 25342114, 56000002, 1, 4075865, 589714,
           612, 136, 7e+05, 699998, 1.94454194361955, 4212.24278192424,
           1.09285712334184, 22511937, 5879812, 1032891, 116, 813412, 418599,
           531, 25342704, 56000002, 1, 4074471, 590417, 612, 145, 7e+05,
           699998, 1.94454194361955, 4217.26405978943, 1.09285708431123,
           22513302, 5879136, 1032043, 131, 813519, 418575, 546, 25342750,
           56000002, 1, 4073423, 591173, 570, 144, 7e+05, 7e+05,
           1.94454194361955, 4222.66428571429, 1.01785714285714, 22513673,
           5877721, 1031419, 136, 813451, 418864, 501, 25344235, 5.6e+07, 1,
           4074409, 590356, 573, 136, 7e+05, 699998, 1.94454194361955,
           4216.82849612806, 1.0232142674426, 22512511, 5880570, 1030853, 123,
           813357, 418644, 500, 25343444, 56000002)
  expect_equivalent(y, array(tmp, dim = c(19L, 5L, 1L)))
})

test_that("there are no infections when beta is 0", {
  pars <- example_gas_parameters()
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


test_that("there are no infections when theta_A is 0", {
  pars <- example_gas_parameters()
  pars$theta_A <- 0
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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
  pars <- example_gas_parameters()
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


test_that("the transmission rate is calculated correctly", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, 365 * 5, 7), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(model_index())
  time <- y["time", 1, -1]
  beta <- y["beta_t", , -1]

  max_date <- model_week_date(time[apply(beta, 1, which.max)])
  # when there is a seasonal effect all max foi dates should be the same and
  # in the same week as the input t_s
  expect_true(all(max_date == max_date[1]))
  expect_true(as.numeric(max_date[1] - model_date(pars$t_s)) %% 365 <= 7)
  # beta does not vary by particle
  expect_equal(colMeans(beta), beta[1, ])
  expect_equal(max(beta), pars$beta * (1 + pars$sigma), tol = 1e-4)
  expect_equal(min(beta), pars$beta * (1 - pars$sigma), tol = 1e-4)

  pars$sigma <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, 365 * 5, 7), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(model_index())

  # when there is no seasonal effect, all beta should be the same
  expect_true(all(y["beta_t", , -1] == pars$beta))
})


test_that("incidence time series output correctly", {
  pars <- example_gas_parameters()
  pars$omega <- 0
  pars$alpha <- 0
  pars$delta_I <- 1e6
  pars$delta_F <- 1e6
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(model_index())
  expect_equal(y["F", , 5], rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(y["F", , 2], y["scarlet_fever_inc", , 2])
  expect_equal(y["I", , 5], rowSums(y["igas_inc", , ]))
  expect_equal(y["I", , 2], y["igas_inc", , 2])
  expect_true(all(y["entrants_inc", , ] == 0))
  expect_true(all(y["leavers_inc", , ] == 0))
  expect_equal(y["pharyngitis_rate", , ] * pars$phi_S,
               y["pharyngitis_inc", , ] / y["N", , ] * 100000)
  expect_equal(y["scarlet_fever_rate", , ],
               y["scarlet_fever_inc", , ] / y["N", , ] * 100000)

  pars$delta_S <- 1e6
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- mcstate::array_bind(arrays = y, along = 3L)
  rownames(y) <- names(model_index())
  expect_equal(y["S2", , 5],
               rowSums(y["pharyngitis_inc", , ] - y["scarlet_fever_inc", , ]))
  expect_equal(y["S2", , 2],
               y["pharyngitis_inc", , 2] - y["scarlet_fever_inc", , 2])
})

test_that("aging works", {
  pars <- example_gas_parameters(10)

  check_compartment_aging <- function(nm, pars) {
    pars$alpha[] <- pars$omega[] <- 0 # no entrants or leavers
    pars$r_age <- 1 / pars$dt # everyone ages a group per week
    pars$beta <- 0 # no transmission
    pars[grep("delta", names(pars))] <- Inf # no movement between compartments
    pars$A0[] <- pars$U0[] <- pars$R0[] <- 0 # clear population
    pars[[paste0(nm, 0)]][1] <- N0 <- 1e4  # everyone initially in group 1

    # run model
    mod <- model$new(pars, 0, 2, seed = 1L)
    y <- lapply(seq(0, 9), mod$simulate)
    y <- mcstate::array_bind(arrays = y, along = 3L)
    rownames(y) <- nms <- names(model_index(10))

    # check no-one moves compartments
    expect_equivalent(y[grep(nm, nms), , ], y[grep("N_", nms), , ])
    # check all demographic changes are deterministic
    expect_equal(y[, 1, ], y[, 2, ])
    ## check pop size is constant
    expect_true(all(apply(y[grep("N_", nms), , ], c(2, 3), sum) == N0))
    ## check everyone advances one group per day
    expect_equivalent(y[grep(nm, nms), 1, ], diag(pars$n_group) * N0)
    ## nothing going on in other compartments
    expect_equal(sum(y[grep(paste0("N_|time|", nm), nms, invert = TRUE), , ]),
                 0)
  }

  for (i in model_compartments()) {
    check_compartment_aging(i, pars)
  }

})
