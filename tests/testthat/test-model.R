

test_that("model runs", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(mod$info()$index)

  # check that states sum to N
  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))

  tmp <- c(1, 4018819, 1369730, 707, 974, 10997, 10998, 1.94454194361955,
           9785.04838901872, 1.26250002254464, 22185442, 5921153, 1039293,
           989017, 638, 600, 25863856, 55999999, 1, 4019221, 1369672, 691,
           950, 10997, 10998, 1.94454194361955, 9784.60553186796,
           1.23392859346301,
           22185224, 5920925, 1039404, 988145, 646, 575, 25865080, 55999999,
           1, 4018070, 1368033, 653, 1004, 10997, 10998, 1.94454194361955,
           9772.83053165769, 1.16607144939413, 22185734, 5921853, 1039906,
           987374, 671, 564, 25863897, 55999999, 1, 4015364, 1369906, 669,
           968, 10997, 10998, 1.94454194361955, 9786.23767475424,
           1.19464287847577,
           22188926, 5917633, 1037129, 987837, 704, 562, 25867208, 55999999,
           1, 4017723, 1371228, 742, 934, 10997, 10998, 1.94454194361955,
           9795.81088921091, 1.32500002366071, 22185562, 5918936, 1038619,
           989248, 661, 621, 25866352, 55999999)
  expect_equivalent(y, array(tmp, dim = c(18L, 5L, 1L)))
})

test_that("there are no infections when beta is 0", {
  pars <- example_gas_parameters()
  pars$beta <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$P == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc == 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})


test_that("there are no infections when theta_A is 0", {
  pars <- example_gas_parameters()
  pars$theta_A <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$P == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc == 0))

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
  expect_true(all(y$S == 0))
  expect_true(all(y$P == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc == 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})


test_that("there are no symptomatic infections when p_S = 0", {
  pars <- example_gas_parameters()
  pars$p_S <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$E == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$P == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc > 0))

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
  expect_true(all(y$S > 0))
  expect_true(all(y$P > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc > 0))

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
  expect_true(all(y$S > 0))
  expect_true(all(y$P > 0))
  expect_true(any(y$F > 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc == 0))

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
  expect_true(all(y$S > 0))
  expect_true(all(y$P == 0))
  expect_true(all(y$F == 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$igas_inc > 0))

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
  expect_true(all(y$S == 0))
  expect_true(all(y$P > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))

  expect_equal(y$N, Reduce(`+`, y[model_compartments()]))
  expect_true(all(unlist(y) >= 0))
})

test_that("no pharyngitis cases are reported when p_T = 0", {
  pars <- example_gas_parameters()
  pars$p_T <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$transform_variables(mod$simulate(7))

  expect_true(all(y$A > 0))
  expect_true(all(y$E > 0))
  expect_true(all(y$S > 0))
  expect_true(all(y$P > 0))
  expect_true(all(y$F > 0))
  expect_true(all(y$R > 0))
  expect_true(all(y$pharyngitis_inc > 0))
  expect_equal(y$pharyngitis_scarlet_fever_rate,
               y$scarlet_fever_inc / y$N * 1e5, tolerance = 1e-6)
  expect_equal(y$pharyngitis_scarlet_fever_rate,
               y$scarlet_fever_rate)

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
  expect_true(all(y$S == 0))
  expect_true(all(y$P == 0))
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
  pars$delta_F <- Inf
  pars$delta_S <- Inf
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(model_index())
  expect_equal(y["F", , 5], rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(y["F", , 2], y["scarlet_fever_inc", , 2])
  expect_equal(y["S", , 5], rowSums(y["pharyngitis_inc", , ]))
  expect_equal(y["S", , 2], y["pharyngitis_inc", , 2])


  expect_true(all(y["births_inc", , ] == 0))
  expect_true(all(y["net_leavers_inc", , ] == 0))
  expect_true(all(y["N", , ] == pars$N0))
  expect_equal(y["pharyngitis_scarlet_fever_rate", , ] * y["N", , ] / 100000,
               y["pharyngitis_inc", , ] / pars$phi_S * pars$p_T +
                 y["scarlet_fever_inc", , ])
  expect_equal(y["scarlet_fever_rate", , ],
               y["scarlet_fever_inc", , ] / y["N", , ] * 100000)

  pars$delta_S <- Inf
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- mcstate::array_bind(arrays = y, along = 3L)
  rownames(y) <- names(model_index())

})

test_that("aging works", {
  pars <- model_parameters(no_gas_parameters(10))
  pars$alpha <- pars$omega[] <- 0
  pars$r_age <- 1 / pars$dt # everyone ages a group per week
  pars$U0[] <- 0 # clear population
  pars[grep("delta", names(pars))] <- Inf # no movement between compartments

  check_compartment_aging <- function(nm, pars) {
    pars[[paste0(nm, 0)]][1] <- N0 <- 1e4 # everyone initially in group 1
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

test_that("aging does not affect model dynamics", {
  pars <- example_gas_parameters()
  pars_a <- example_gas_parameters(2)
  pars$alpha <- pars_a$alpha <- 0
  pars$omega[] <- pars_a$omega[] <- 0

  f <- function(pars) {
    mod <- model$new(pars, 0, 5, seed = 1L)
    y <- lapply(seq(0, 9), mod$simulate)
    y <- mcstate::array_bind(arrays = y, along = 3L)
    mod$transform_variables(y)
  }

  y <- f(pars)
  y_a <- f(pars_a)

  absdiff <- function(state) {
    x <- drop(y[[state]])
    y <- colSums(y_a[[state]])
    max(abs(x - y) / x, na.rm = TRUE)
  }

  # check larger compartments are withing 1% (others too stochastic)
  expect_true(absdiff("U") < 0.01)
  expect_true(absdiff("A") < 0.01)
  expect_true(absdiff("R") < 0.01)
  expect_true(absdiff("infections_inc") < 0.01)
  expect_true(absdiff("pharyngitis_inc") < 0.01)

  expect_equal(Reduce("+", y[model_compartments()]), y$N)
  expect_equal(Reduce("+", y_a[model_compartments()]), y_a$N)
  expect_equal(y$beta_t, y_a$beta_t)

  # check all compartments equal
  expect_equal(sum(y$net_leavers_inc), 0)
  expect_equal(sum(y$births_inc), 0)
  expect_equal(sum(y_a$net_leavers_inc), 0)
  expect_equal(sum(y_a$births_inc), 0)
})

test_that("time-varying births works", {
  pars <- model_parameters(no_gas_parameters(1))
  pars$alpha <- rep(seq_len(52), each = 7) * 7

  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, length(pars$alpha), 7), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(model_index())

  expect_equal(y["births_inc", 1, ], seq(0, 364, 7))

  # check that states sum to N
  expect_equivalent(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})
