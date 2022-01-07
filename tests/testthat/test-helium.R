test_that("helium age groups are parametrised correctly", {
  groups <- helium_age_groups()
  expect_equal(length(groups$age_start), groups$n_group)
  expect_equal(length(groups$age_end), groups$n_group)
  expect_true(all(groups$age_end > groups$age_start))
})
