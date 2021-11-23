

test_that("model runs", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- mod$simulate(7)
  rownames(y) <- names(mod$info()$index)

  # check that states sum to N
  expect_equal(colSums(y[model_compartments(), , ]), y["N", , ])
  expect_true(all(y >= 0))


  y <- c(1, 21687448, 6174233, 1112981, 159, 950850, 569569,
           749, 25514822, 56010811, 4950895, 897364, 892, 179, 700590, 689779,
           0.254421609239827, 1, 21685182, 6176629, 1112823, 169, 949518,
           569192, 776, 25515215, 56009504, 4951768, 896877, 910, 182, 699526,
           690022, 0.254437293682947, 1, 21683343, 6178466, 1112581, 150,
           951060, 570309, 741, 25513542, 56010192, 4954365, 898339, 899,
           177, 699812, 689620, 0.25453578045828, 1, 21688454, 6173169,
           1112842, 170, 950277, 569084, 720, 25515656, 56010372, 4949707,
           897336, 857, 208, 699963, 689591, 0.254319583239358, 1, 21680889,
           6178399, 1113810, 160, 951629, 570990, 747, 25513261, 56009885,
           4956400, 897754, 892, 180, 700225, 690340, 0.254589465871049)
  expect_equivalent(y, array(y, dim = c(17L, 5L, 1L)))
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


test_that("the foi is calculated correctly", {
  pars <- example_gas_parameters()
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, 365 * 5, 7), mod$simulate)
  y <- abind::abind(y, along = 3L)
  rownames(y) <- names(model_index())
  t <- y["time", 1, ]
  par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
  N <- y["N", , ]
  infectious <- y["A", , ] + y["S1", , ] + y["S2", , ]
  I <- y["E", , ] + infectious + y["F", , ] + y["I", , ]
  matplot(t, t(y["U", , ] / N), type = "l", ylim = c(0, 1),
          col = 2, ylab = "%")
  matlines(t, t(I / N), type = "l", col = 3)
  matlines(t, t(y["R", , ] / N), type = "l", col = 4)
  legend("topright", legend = c("S", "I", "R"), fill = 2:4, bty = "n")

  foi <- t(y["foi", , ] / infectious * N)
  max_date <- model_week_date(t[apply(foi, 2, which.max)])
  # when there is a seasonal effect all max foi dates should be the same and
  # in the same week as the input t_s
  expect_true(all(max_date == max_date[1]))
  expect_true(as.numeric(max_date[1] - model_date(pars$t_s)) %% 365 <= 7)

  matplot(t, foi, type = "l", col = 1)
  matplot(t, t(y["igas_inc", , ]), type = "l", col = 1, ylab = "iGAS / 100k")
  matplot(t, t(y["scarlet_fever_inc", , ]), type = "l", col = 1,
          ylab = "SF / 100k")

  pars$sigma <- 0
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, 365 * 5, 7), mod$simulate)
  y <- abind::abind(y, along = 3L)
  rownames(y) <- names(model_index())
  infectious <- y["A", , ] + y["S1", , ] + y["S2", , ]
  foi <- t(y["foi", , ] / infectious * y["N", , ])
  max_date <- model_week_date(t[apply(foi, 2, which.max)])

  # when there is no seasonal effect, date of max foi should differ by particle
  expect_true(any(max_date != max_date[1]))
  expect_true(as.numeric(max_date[1] - model_date(pars$t_s)) %% 365 > 7)
})


test_that("incidence time series output correctly", {
  pars <- example_gas_parameters()
  pars$omega <- 0
  pars$alpha <- 0
  pars$delta_I <- 1e6
  pars$delta_F <- 1e6
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- abind::abind(y, along = 3L)
  rownames(y) <- names(model_index())
  expect_equal(y["F", , 5], rowSums(y["scarlet_fever_inc", , ]))
  expect_equal(y["F", , 2], y["scarlet_fever_inc", , 2])
  expect_equal(y["I", , 5], rowSums(y["igas_inc", , ]))
  expect_equal(y["I", , 2], y["igas_inc", , 2])
  expect_true(all(y["entrants_inc", ,] == 0))
  expect_true(all(y["leavers_inc", ,] == 0))

  pars$delta_S <- 1e6
  mod <- model$new(pars, 0, 5, seed = 1L)
  y <- lapply(seq(0, by = 7, length.out = 5), mod$simulate)
  y <- abind::abind(y, along = 3L)
  rownames(y) <- names(model_index())
  expect_equal(y["S2", , 5],
               rowSums(y["pharyngitis_inc", , ] - y["scarlet_fever_inc", , ]))
  expect_equal(y["S2", , 2],
               y["pharyngitis_inc", , 2] - y["scarlet_fever_inc", , 2])
})
