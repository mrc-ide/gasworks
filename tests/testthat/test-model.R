
test_that("model runs", {
  pars <- example_parameters()
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  # check that states sum to N
  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))

  ## one row per state, one column per particle
  expect_equal(dim(y), c(length(model_index()), 5, 1))
  expect_equal(rownames(y), names(model_index()))

  ## particles differ, but are reproducible given the seed
  expect_false(all(y[, 1, ] == y[, 2, ]))
  expect_identical(model_simulate(pars, 7, n_particles = 5, seed = 1L), y)
  expect_false(identical(model_simulate(pars, 7, n_particles = 5, seed = 2L),
                         y))

  ## age-specific outputs are only produced by the helium model
  expect_true(all(y[grep("_[0-9]+(_[0-9]+)?$", rownames(y)), , ] == 0))
})

test_that("deterministic model is close to mean of stochastic model", {
  pars <- example_parameters()
  times <- seq(0, 28, by = 7)
  y_det <- model_simulate(pars, times, deterministic = TRUE)
  y <- model_simulate(pars, times, n_particles = 200, seed = 1L, n_threads = 2)

  expect_identical(model_simulate(pars, times, deterministic = TRUE), y_det)
  nms <- c("U", "A", "R", "N", "infections_inc", "gas_pharyngitis_inc",
           "igas_inc")
  expect_equal(apply(y[nms, , ], c(1, 3), mean), y_det[nms, 1, ],
               tolerance = 0.01)
})

test_that("there are no infections when beta is 0", {
  pars <- example_parameters()
  pars$beta <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)


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
  pars <- example_parameters()
  pars$theta_A <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))
  expect_true(all(y["prev_A", , ] > 0))
  expect_true(all(y["prev_R", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there are no infections when A0 = 0", {
  pars <- example_parameters()
  pars$A0[] <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y["E", , ] == 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] == 0))
  expect_true(all(y["prev_A", , ] == 0))
  expect_true(all(y["prev_R", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})


test_that("there are no symptomatic infections when p_S = 0", {
  pars <- example_parameters()
  pars$p_S <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

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
  pars <- example_parameters()
  pars$p_S <- 1
  pars$E0[] <- pars$A0
  pars$A0[] <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

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
  pars <- example_parameters()
  pars$p_I <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

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
  pars <- example_parameters()
  pars$p_F <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] > 0))
  expect_true(all(y["S2", , ] > 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_inc", , ] == 0))
  expect_true(all(y["daily_gas_pharyngitis_rate", , ] > 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there is no pharyngitis when p_F = 1", {
  pars <- example_parameters()
  pars$p_F <- 1
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] > 0))
  expect_true(all(y["F1", , ] > 0))
  expect_true(all(y["F2", , ] > 0))
  expect_true(all(y["R", , ] > 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_inc", , ] > 0))
  expect_true(all(y["daily_gas_pharyngitis_rate", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("there is no immunity when p_S = 0 and p_R = 0", {
  pars <- example_parameters()
  pars$p_S <- 0
  pars$p_R <- 0
  pars$R0[] <- 0
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y["A", , ] > 0))
  expect_true(all(y["S1", , ] == 0))
  expect_true(all(y["S2", , ] == 0))
  expect_true(all(y["P", , ] == 0))
  expect_true(all(y["F1", , ] == 0))
  expect_true(all(y["F2", , ] == 0))
  expect_true(all(y["R", , ] == 0))
  expect_true(all(y["igas_inc", , ] > 0))
  expect_true(all(y["scarlet_fever_inc", , ] == 0))
  expect_true(all(y["gas_pharyngitis_inc", , ] == 0))
  expect_true(all(y["daily_gas_pharyngitis_rate", , ] == 0))

  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})


test_that("the transmission rate is calculated correctly", {
  pars <- example_parameters()
  y <- model_simulate(pars, seq(0, 365 * 5, 7), n_particles = 5, seed = 1L)
  time <- seq(0, 365 * 5, 7)[-1] / 7
  beta <- y["beta_t", , -1]

  max_date <- model_week_date(time[apply(beta, 1, which.max)])
  # when there is a seasonal effect all max foi dates should be the same and
  # in the same week as the input t_s
  expect_true(all(max_date == max_date[1]))
  expect_true(as.numeric(max_date[1] - model_date(pars$t_s)) %% 365 <= 7)
  # beta does not vary by particle
  expect_equal(colMeans(beta), beta[1, ])
  expect_equal(max(beta), pars$beta * (1 + pars$sigma), tolerance = 1e-4)
  expect_equal(min(beta), pars$beta * (1 - pars$sigma), tolerance = 1e-4)

  pars$sigma <- 0
  y <- model_simulate(pars, seq(0, 365 * 5, 7), n_particles = 5, seed = 1L)

  # when there is no seasonal effect, all beta should be the same
  expect_true(all(y["beta_t", , -1] == pars$beta))
})


test_that("incidence time series output correctly", {
  pars <- example_parameters()
  pars$omega <- 0
  pars$alpha <- 0
  pars$delta_F <- Inf
  pars$delta_S <- Inf
  pars$p_T <- 0.1
  y <- model_simulate(pars, seq(0, by = 7, length.out = 5),
                      n_particles = 5, seed = 1L)
  expect_equal(y["F1", , 5] + y["F2", , 5], rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(y["F1", , 2] + y["F2", , 2], y["scarlet_fever_inc", , 2])
  expect_equal(y["S1", , 5] + y["S2", , 5],
               rowSums(y["gas_pharyngitis_inc", , ]))
  expect_equal(y["S1", , 2] + y["S2", , 5], y["gas_pharyngitis_inc", , 2])

  ## check scarlet_fever_cases is calculated correctly
  ## (the model rounds half away from zero, R rounds half to even)
  expect_equal(y["scarlet_fever_cases", , ],
               floor(pars$q_F * y["scarlet_fever_inc", , ] + 0.5))

  expect_true(all(y["births_inc", , ] == 0))
  expect_true(all(y["net_leavers_inc", , ] == 0))
  expect_true(all(y["N", , ] == pars$N0))
  ## output rates are daily so multiply by 7 for weekly
  ## divide by 1e5 for per person
  ## multiply by population size

  expect_equal(y["daily_gas_pharyngitis_rate", , ] * 7 / 1e5 * y["N", , ],
    y["gas_pharyngitis_inc", , ] * pars$p_T)
  expect_equal(y["daily_scarlet_fever_rate", , ] * 7 / 1e5 * y["N", , ],
               y["scarlet_fever_inc", , ])
})

test_that("etiologic_fraction outputs_correctly", {
  pars <- example_parameters()
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)
  expect_equal(y["etiologic_fraction", , ], rep(pars$phi_S, 5))

  pars <- example_parameters(16)
  y <- model_simulate(pars, 7, n_particles = 5, seed = 1L)

  expect_true(all(y[grep("^etio", rownames(y)), , ] == pars$phi_S[1]))
})

test_that("can vary probabililty of reporting over time", {
  pars <- example_parameters()
  pars$q_F <- seq(0.3, 0.6, length.out = 7)

  y <- model_simulate(pars, seq(0, 7), n_particles = 5, seed = 1L)

  ## check scarlet_fever_cases is calculated correctly
  ## include adjustment for rounding conventions (0.5 -> up vs 0.5 -> even)
  expect_equal(t(y["scarlet_fever_cases", , ]),
               floor(c(0, pars$q_F) * t(y["scarlet_fever_inc", , ]) + 0.5))
})


test_that("rates are calculated correctly when n_group == 16", {
  pars <- example_parameters(16)
  pars$omega[] <- 0
  pars$alpha <- 0
  pars$delta_F <- Inf
  pars$delta_S <- Inf
  pars$p_T <- 0.1
  pars$prev_A <- pars$prev_A <- seq(0.1, 0.2, length.out = 16)
  ip <- initial_parameters(pars)
  pars[names(ip)] <- ip
  y <- model_simulate(pars, seq(0, by = 7, length.out = 5),
                      n_particles = 5, seed = 1L)
  nms <- rownames(y)

  expect_equal(colSums(y[grep("^F", nms), , 5]),
               rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(colSums(y[grep("^F", nms), , 2]), y["scarlet_fever_inc", , 2])
  expect_equal(colSums(y[grep("^S", nms), , 5]),
               rowSums(y["gas_pharyngitis_inc", , ]))
  expect_equal(colSums(y[grep("^S", nms), , 2]), y["gas_pharyngitis_inc", , 2])
  expect_equivalent(pars$prev_A,
               y[grep("^prev_A", nms), 1, 1], tolerance = 1e-6)
  expect_equivalent((1 - pars$prev_A) * pars$prev_R,
                    y[grep("^prev_R", nms), 1, 1], tolerance = 1e-6)

  expect_true(all(y["births_inc", , ] == 0))
  expect_true(all(y["net_leavers_inc", , ] == 0))
  N <- y[grep("^N", nms), , ]
  expect_equivalent(N, array(pars$N0, dim(N)))

  ## output rates are daily so multiply by 7 for weekly
  ## divide by 1e5 for per person
  ## multiply by population size

  expect_equivalent(y["daily_gas_pharyngitis_rate", , ] * 7 / 1e5 * colSums(N),
    y["gas_pharyngitis_inc", , ] * pars$p_T)

  age <- helium_age_groups()$age_start
  N_04   <- colSums(N[age < 5, , , drop = FALSE])
  N_5_14  <- colSums(N[age >= 5 & age < 15, , , drop = FALSE])
  N_15_44 <- colSums(N[age >= 15 & age < 45, , , drop = FALSE])
  N_45_64 <- colSums(N[age >= 45 & age < 65, , , drop = FALSE])
  N_65_74 <- colSums(N[age >= 65 & age < 75, , , drop = FALSE])
  N_75    <- colSums(N[age >= 75, , , drop = FALSE])


  expect_equivalent(N_04 + N_5_14 + N_15_44 + N_45_64 + N_65_74 + N_75,
                    colSums(N))

  expect_equivalent(y["gas_pharyngitis_inc", , ] * pars$p_T,
    (y["daily_gas_pharyngitis_rate_04", , ] * N_04 +
       y["daily_gas_pharyngitis_rate_05_14", , ] * N_5_14 +
       y["daily_gas_pharyngitis_rate_15_44", , ] * N_15_44 +
       y["daily_gas_pharyngitis_rate_45_64", , ] * N_45_64 +
       y["daily_gas_pharyngitis_rate_65_74", , ] * N_65_74 +
       y["daily_gas_pharyngitis_rate_75", , ] * N_75) * 7 / 1e5)

  expect_equivalent(y["scarlet_fever_inc", , ],
    (y["daily_scarlet_fever_rate_04", , ] * N_04 +
       y["daily_scarlet_fever_rate_05_14", , ] * N_5_14 +
       y["daily_scarlet_fever_rate_15_44", , ] * N_15_44 +
       y["daily_scarlet_fever_rate_45_64", , ] * N_45_64 +
       y["daily_scarlet_fever_rate_65_74", , ] * N_65_74 +
       y["daily_scarlet_fever_rate_75", , ] * N_75) * 7 / 1e5)

  expect_equivalent(y["scarlet_fever_inc", , ],
       y["scarlet_fever_inc_04", , ] +
       y["scarlet_fever_inc_05_14", , ] +
       y["scarlet_fever_inc_15_44", , ] +
       y["scarlet_fever_inc_45_64", , ] +
       y["scarlet_fever_inc_65_74", , ] +
       y["scarlet_fever_inc_75", , ])

})

test_that("aging works", {
  pars <- model_parameters(no_gas_parameters(10))
  pars$alpha <- pars$omega[] <- 0
  pars$r_age <- 7 # everyone ages a group per day
  pars[grep("delta", names(pars))] <- Inf # no movement between compartments
  pars[grep("^k_", names(pars))] <- 1 # all exponential
  init <- initial_parameters(pars)
  pars[names(init)] <- init
  pars$U0[] <- 0 # clear population

  check_compartment_aging <- function(nm, pars) {

    pars[[paste0(nm, 0)]][1] <- N0 <- 1e4 # everyone initially in group 1
    # run model
    y <- model_simulate(pars, seq(0, 9), n_particles = 2, seed = 1L)
    nms <- rownames(y)

    # check no-one moves compartments
    expect_equivalent(y[startsWith(nms, nm), , ], y[grep("^N", nms), , ])
    # check all demographic changes are deterministic
    expect_equal(y[, 1, ], y[, 2, ])
    ## check pop size is constant

    expect_true(all(apply(y[grep("^N", nms), , ], c(2, 3), sum) == N0))
    ## check everyone advances one group per day
    expect_equivalent(y[startsWith(nms, nm), 1, ], diag(pars$n_group) * N0)
    ## nothing going on in other compartments

    expect_equal(sum(y[grep(paste0("^N|time|prev|etio|", nm), nms,
                            invert = TRUE), , ]), 0)
  }

  for (i in unique(substr(model_compartments(), 1, 1))) {
    check_compartment_aging(i, pars)
  }
})

test_that("aging does not affect model dynamics", {
  pars <- example_parameters()
  pars_a <- example_parameters(2)
  pars$alpha <- pars_a$alpha <- 0
  pars$omega[] <- pars_a$omega[] <- 0

  f <- function(pars) {
    simulate_unpacked(pars, seq(0, 9))
  }

  y <- f(pars)
  y_a <- f(pars_a)

  ## total over age groups and sub-compartments, by particle and time
  total <- function(x) {
    n <- length(dim(x))
    if (n == 2) x else apply(x, c(n - 1, n), sum)
  }

  absdiff <- function(state) {
    x <- total(y[[state]])
    z <- total(y_a[[state]])
    max(abs(x - z) / x, na.rm = TRUE)
  }

  # check larger compartments are within 1% (others too stochastic)
  expect_true(absdiff("U") < 0.01)
  expect_true(absdiff("A") < 0.01)
  expect_true(absdiff("R") < 0.01)
  expect_true(absdiff("infections_inc") < 0.01)
  expect_true(absdiff("gas_pharyngitis_inc") < 0.01)

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

  y <- model_simulate(pars, seq(0, length(pars$alpha), 7),
                      n_particles = 5, seed = 1L)

  expect_equal(y["births_inc", 1, ], seq(0, 364, 7))

  # check that states sum to N
  expect_equivalent(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))
})

test_that("can vary erlang compartments", {
  pars <- example_parameters(1)
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

  y <- model_simulate(pars, seq(0, 52, 7), n_particles = 5, seed = 1L)

  expect_equal(length(grep("^S", rownames(y))), pars$k_S)
  expect_equal(length(grep("^A", rownames(y))), pars$k_A)
  expect_equal(length(grep("^E", rownames(y))), pars$k_E)
  expect_equal(length(grep("^P", rownames(y))), pars$k_P)
  expect_equal(length(grep("^F", rownames(y))), pars$k_F)
  expect_equal(length(grep("^R", rownames(y))), pars$k_R)

  # check that states sum to N
  expect_true(all(y["N", , ] == pars$N0))
  expect_true(all(y >= 0))
})
