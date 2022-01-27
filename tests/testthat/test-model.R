

test_that("model runs", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  # check that states sum to N
  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))

  tmp <- c(1, 4045133, 1373194, 727, 938, 10997, 10996, 1.94454194361955,
           9809.82678571429, 1.29821428571429, 22159247, 56000001, 5928661,
           1052289, 775606, 416796, 704, 542, 166, 25665990, 1, 4044289,
           1372467, 692, 939, 10997, 10997, 1.94454194361955, 9804.57142857143,
           1.23571428571429, 22160502, 5.6e+07, 5930566, 1051231, 775591,
           415890, 736, 513, 147, 25664824, 1, 4047132, 1371527, 678, 1016,
           10997, 10997, 1.94454194361955, 9797.83214285714, 1.21071428571429,
           22156462, 5.6e+07, 5933201, 1054026, 775298, 416302, 704, 519,
           139, 25663349, 1, 4047245, 1373647, 698, 926, 10997, 10996,
           1.94454194361955,
           9813.01071428571, 1.24642857142857, 22157315, 56000001, 5930844,
           1052602, 776369, 416513, 710, 513, 157, 25664978, 1, 4045950,
           1373480, 680, 939, 10997, 10996, 1.94454194361955, 9811.78571428571,
           1.21428571428571, 22157734, 56000001, 5928868, 1053003, 775962,
           416807, 677, 501, 151, 25666298)
  expect_equivalent(y, array(tmp, dim = c(20L, 5L, 1L)))
})

test_that("there are no infections when beta is 0", {
  pars <- example_gas_parameters()
  pars$beta <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())


  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})


test_that("there are no infections when theta_A is 0", {
  pars <- example_gas_parameters()
  pars$theta_A <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there are no infections when A0 = 0", {
  pars <- example_gas_parameters()
  pars$A0[] <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})


test_that("there are no symptomatic infections when p_S = 0", {
  pars <- example_gas_parameters()
  pars$p_S <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there are no asymptomatic infections when p_S = 1", {
  pars <- example_gas_parameters()
  pars$p_S <- 1
  pars$E0[] <- pars$A0
  pars$A0[] <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] == 0))
  expect_true(all(y["S1", , ] > 0))
  expect_true(all(y["S2", , ] > 0))
  expect_true(all(y["P", , ] > 0))
  expect_true(all(y["F1", , ] > 0))
  expect_true(all(y["F2", , ] > 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there is no iGAS when p_I = 0", {
  pars <- example_gas_parameters()
  pars$p_I <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] > 0))
  expect_true(all(y["S2", , ] > 0))
  expect_true(all(y["P", , ] > 0))
  expect_true(all(y["F1", , ] > 0))
  expect_true(all(y["F2", , ] > 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there is no scarlet fever when p_F = 0", {
  pars <- example_gas_parameters()
  pars$p_F <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] > 0))
  expect_true(all(y["S2", , ] > 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_rate", , ] == 0))
  expect_true(all(y["pharyngitis_scarlet_fever_rate", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there is no pharyngitis when p_F = 1", {
  pars <- example_gas_parameters()
  pars$p_F <- 1
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] > 0))
  expect_true(all(y["F1", , ] > 0))
  expect_true(all(y["F2", , ] > 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_rate", , ] > 0))
  expect_true(all(y["pharyngitis_inc", , ] == 0))
  expect_true(all(y["pharyngitis_scarlet_fever_rate", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("no pharyngitis cases are reported when p_T = 0", {
  pars <- example_gas_parameters()
  pars$p_T <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] > 0))
  expect_true(all(y["S2", , ] > 0))
  expect_true(all(y["P", , ] > 0))
  expect_true(all(y["F1", , ] > 0))
  expect_true(all(y["F2", , ] > 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_rate", , ] > 0))
  expect_true(all(y["pharyngitis_inc", , ] > 0))
  expect_equal(y["pharyngitis_scarlet_fever_rate", , ],
               y["scarlet_fever_rate", , ])
  expect_equal(y["pharyngitis_scarlet_fever_rate", , ],
               y["scarlet_fever_inc", , ] / y["N", , ] * 1e5, tolerance = 1e-6)

})

test_that("there is no immunity when p_S = 0 and p_R = 0", {
  pars <- example_gas_parameters()
  pars$p_S <- 0
  pars$p_R <- 0
  pars$R0[] <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(model_index())

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] == 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_rate", , ] == 0))
  expect_true(all(y["pharyngitis_inc", , ] == 0))
  expect_true(all(y["pharyngitis_scarlet_fever_rate", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
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
  expect_equal(y["F1", , 5] + y["F2", , 5], rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(y["F1", , 2] + y["F2", , 2], y["scarlet_fever_inc", , 2])
  expect_equal(y["S1", , 5] + y["S2", , 5], rowSums(y["pharyngitis_inc", , ]))
  expect_equal(y["S1", , 2] + y["S2", , 5], y["pharyngitis_inc", , 2])


  expect_true(all(y["births_inc", , ] == 0))
  expect_true(all(y["net_leavers_inc", , ] == 0))
  expect_true(all(y["N", , ] == pars$N0))
  expect_equal(y["pharyngitis_scarlet_fever_rate", , ] * y["N", , ] / 100000,
               y["pharyngitis_inc", , ] / pars$phi_S * pars$p_T +
                 y["scarlet_fever_inc", , ])
  expect_equal(y["scarlet_fever_rate", , ],
               y["scarlet_fever_inc", , ] / y["N", , ] * 100000)
})

test_that("aging works", {
  pars <- model_parameters(no_gas_parameters(10))
  pars$alpha <- pars$omega[] <- 0
  pars$r_age <- 1 / pars$dt # everyone ages a group per week
  pars[grep("delta", names(pars))] <- Inf # no movement between compartments
  pars[grep("^k_", names(pars))] <- 1 # all exponential
  init <- initial_parameters(pars)
  pars[names(init)] <- init
  pars$U0[] <- 0 # clear population

  check_compartment_aging <- function(nm, pars) {
    pars[[paste0(nm, 0)]][1] <- N0 <- 1e4 # everyone initially in group 1
    # run model
    mod <- model$new(pars, 0, 2, seed = 1L)
    y <- lapply(seq(0, 9), mod$simulate)
    y <- mcstate::array_bind(arrays = y, along = 3L)
    rownames(y) <- nms <- names(unlist(mod$info()$index))

    # check no-one moves compartments
    expect_equivalent(y[grep(nm, nms), , ], y[grep("^N", nms), , ])
    # check all demographic changes are deterministic
    expect_equal(y[, 1, ], y[, 2, ])
    ## check pop size is constant

    expect_true(all(apply(y[grep("^N", nms), , ], c(2, 3), sum) == N0))
    ## check everyone advances one group per day
    expect_equivalent(y[grep(nm, nms), 1, ], diag(pars$n_group) * N0)
    ## nothing going on in other compartments

    expect_equal(sum(y[grep(paste0("^N|time|", nm), nms, invert = TRUE), , ]),
                 0)
  }

  for (i in unique(substr(model_compartments(), 1, 1))) {
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
    x <- y[[state]]
    z <- apply(y_a[[state]], seq(2, length(dim(y_a[[state]]))), sum)
    dim(z) <- c(1, dim(z))
    max(abs(x - z) / x, na.rm = TRUE)
  }

  # check larger compartments are within 1% (others too stochastic)
  expect_true(absdiff("U") < 0.01)
  expect_true(absdiff("A") < 0.01)
  expect_true(absdiff("R") < 0.01)
  expect_true(absdiff("infections_inc") < 0.01)
  expect_true(absdiff("pharyngitis_inc") < 0.01)

  expect_true(all(y$N == pars$N0))
  expect_true(all(y_a$N == pars_a$N0))
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

test_that("can vary erlang compartments", {
  pars <- model_parameters(example_gas_parameters(1))
  pars$alpha[] <- 0
  pars$omega[] <- 0
  pars$k_A <- 2
  pars$k_E <- 3
  pars$k_F <- 4
  pars$k_P <- 5
  pars$k_R <- 6
  pars$k_S <- 7
  init_pars <- initial_parameters(pars)
  pars[names(init_pars)] <- init_pars

  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, 52, 7), mod$simulate)
  y <- mcstate::array_bind(arrays = y)
  rownames(y) <- names(unlist(mod$info()$index))

  expect_equal(length(grep("S", rownames(y))), pars$k_S)
  expect_equal(length(grep("A", rownames(y))), pars$k_A)
  expect_equal(length(grep("E", rownames(y))), pars$k_E)
  expect_equal(length(grep("P", rownames(y))), pars$k_P)
  expect_equal(length(grep("F", rownames(y))), pars$k_F)
  expect_equal(length(grep("R", rownames(y))), pars$k_R)

  # check that states sum to N
  expect_true(all(y["N", , ] == pars$N0))
  expect_true(all(y >= 0))
})
