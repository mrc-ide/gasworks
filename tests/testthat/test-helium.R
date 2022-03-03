test_that("helium age groups are parametrised correctly", {
  groups <- helium_age_groups()
  expect_equal(16L, groups$n_group)
  expect_equal(length(groups$age_start), groups$n_group)
  expect_equal(length(groups$age_end), groups$n_group)
  expect_true(all(groups$age_end > groups$age_start))
})

test_that("helium_index", {
  p <- example_gas_parameters(19)
  mod <- model$new(p, 0, 10)
  idx <- helium_index(mod$info())

  expect_equal(names(idx), c("run", "state"))
  state_nms <- c("scarlet_fever_inc", "igas_inc", "daily_pharyngitis_rate",
                 paste0("N", seq_len(19)),
                 "pharyngitis_prop_04", "pharyngitis_prop_05_14",
                 "pharyngitis_prop_15_44", "pharyngitis_prop_45_64",
                 "pharyngitis_prop_65_74", "pharyngitis_prop_75",
                 "scarlet_fever_prop_04", "scarlet_fever_prop_05_14",
                 "scarlet_fever_prop_15_44", "scarlet_fever_prop_45_64",
                 "scarlet_fever_prop_65_74", "scarlet_fever_prop_75")
  expect_equal(names(idx$run), state_nms)
  expect_true(all(state_nms %in% names(idx$state)))

  # check can only use on helium model
  mod <- model$new(example_gas_parameters(2), 0, 10)
  expect_error(helium_index(mod$info()))
})
