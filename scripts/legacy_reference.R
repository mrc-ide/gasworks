#!/usr/bin/env Rscript
## Generate reference output from the last version of gasworks built on
## odin.dust, dust and mcstate (0.3.6), used in tests/testthat/test-legacy.R to
## check that the odin2 port reproduces it.
##
## This must be run from the package root, with gasworks 0.3.6 (commit 74be5c2)
## and its dependencies installed. Deterministic runs are used, as the order of
## random number draws differs between dust and dust2.
stopifnot(packageVersion("gasworks") == "0.3.6")
library(gasworks)

convert_pars <- function(pars) {
  ## Durations were converted to weeks in model_parameters(), and are now
  ## passed to the model in days
  i <- grep("^delta_", names(pars))
  pars[i] <- lapply(pars[i], function(x) x / pars$dt)
  pars$dt <- NULL
  pars
}

reference <- function(pars, compare, fitted, index, observed_missing = NULL) {
  times <- seq(0, 70, by = 7)
  mod <- model$new(pars, 0, 1, deterministic = TRUE)
  y <- drop(mod$simulate(times))
  rownames(y) <- names(unlist(mod$info()$index))

  ## Make up some data close to the modelled values in week 4
  ## (two identical columns, as helium_compare required a matrix of particles)
  state <- y[, c(5, 5)]
  nms <- sub("pharyngitis", "gas_pharyngitis", fitted)
  observed <- setNames(as.list(state[nms, 1] * 1.1), fitted)
  i <- grep("pharyngitis", fitted)
  observed[i] <- lapply(observed[i], function(x) x / 0.3)
  observed$scarlet_fever_cases <- round(observed$scarlet_fever_cases)
  observed$igas_inc <- round(observed$igas_inc)
  observed[observed_missing] <- NA

  set.seed(1)
  ll <- compare(state[names(index(mod$info())$run), ], observed, pars)[[1]]

  list(pars = convert_pars(pars), times = times,
       y = y[rownames(y) != "time", ], observed = observed, ll = ll)
}

pars_hydrogen <- example_parameters(1)
pars_hydrogen$alpha <- round(seq(9000, 12000, length.out = 25))
pars_hydrogen$q_F <- seq(0.2, 0.6, length.out = 40)
pars_hydrogen$k_gp <- 2

transform <- helium_create_transform(NULL)
pars_helium <- transform(gasworks:::example_helium_parameters())
pars_helium$alpha <- round(seq(9000, 12000, length.out = 25))
pars_helium$r_age <- 1 / (5 * 52)
pars_helium$omega <- seq(-1e-4, 5e-4, length.out = 16)
pars_helium$k_gp <- 2

by_age <- grep("daily_pharyngitis_rate_", helium_fitted_states(), value = TRUE)

ret <- list(
  hydrogen = reference(pars_hydrogen, hydrogen_compare,
                       hydrogen_fitted_states(), hydrogen_index),
  helium = reference(pars_helium, helium_compare,
                     helium_fitted_states(), helium_index),
  helium_aggregated = reference(pars_helium, helium_compare,
                                helium_fitted_states(), helium_index, by_age))
saveRDS(ret, "tests/testthat/legacy-reference.rds", version = 2)
