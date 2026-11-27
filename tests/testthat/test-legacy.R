## Reference output created with gasworks 0.3.6, the last version built on
## odin.dust, dust and mcstate; see scripts/legacy_reference.R

test_that("deterministic trajectories match gasworks 0.3.6", {
  reference <- readRDS("legacy-reference.rds")
  for (ref in reference[c("hydrogen", "helium")]) {
    y <- model_simulate(ref$pars, ref$times, deterministic = TRUE)
    expect_setequal(setdiff(rownames(ref$y), rownames(y)), character())
    expect_equal(y[rownames(ref$y), 1, ], ref$y, tolerance = 1e-12)
  }
})

test_that("compiled comparison matches gasworks 0.3.6", {
  reference <- readRDS("legacy-reference.rds")
  prepare <- list(hydrogen = hydrogen_prepare_data,
                  helium = helium_prepare_data,
                  helium_aggregated = helium_prepare_data)
  for (nm in names(reference)) {
    ref <- reference[[nm]]
    sys <- system_at_week(ref$pars, 4, n_particles = 1, deterministic = TRUE)
    ll <- compare_compiled(sys, ref$observed, prepare[[nm]])
    ## The legacy value includes exponential noise with mean 1e-6
    expect_equal(ll, ref$ll, tolerance = 1e-6)
  }
})
