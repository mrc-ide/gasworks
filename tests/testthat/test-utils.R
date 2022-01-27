context("utils")

test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})

test_that("named_list", {
  expect_equal(named_list(letters[1:3], 1:3), list(a = 1, b = 2, c = 3))
  expect_equal(named_list(letters[1:3]), list(a = NULL, b = NULL, c = NULL))
})
