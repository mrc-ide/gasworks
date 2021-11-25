context("utils (assert)")

test_that("assert_is", {
  expect_error(assert_is("x", "foo"), "must be a foo")
  expect_silent(assert_is(structure("x", class = "foo"), "foo"))
})


test_that("assert_named", {
  expect_error(assert_named(1), "must be named")
  expect_error(assert_named(setNames(1:2, c("a", "a")), TRUE),
               "must have unique names")
  expect_silent(assert_named(setNames(1:2, c("a", "a")), FALSE))
})


test_that("assert_character", {
  expect_silent(assert_character("string"))
  expect_error(assert_character(1), "must be a character")
  expect_error(assert_character(TRUE), "must be a character")
})


test_that("assert_integer", {
  expect_error(assert_integer(pi), "'pi' must be an integer")
  expect_error(assert_integer(1L, 2), "'1L' must be of length 2")
  expect_identical(assert_integer(1L), 1L)
  expect_identical(assert_integer(1.0), 1L)
  expect_identical(assert_integer(1 + 1e-15), 1L)
})

test_that("assert_length", {
  x <- c(1, 1)
  expect_error(assert_length(x, 1), "'x' must be of length 1")
  expect_identical(assert_length(x), x)
  expect_identical(assert_length(x, 2), x)
})

test_that("assert_strictly_increasing", {
  x <- c(0, 1, 2)
  expect_silent(assert_strictly_increasing(x))
  expect_error(assert_strictly_increasing(x, 2), "'x' must be of length 2")
  expect_error(assert_strictly_increasing(c(0, 0, 1)),
               "must be strictly increasing")
  expect_error(assert_strictly_increasing(c(0, -1, -2)),
               "must be strictly increasing")
})

test_that("assert_unit_interval", {
  expect_identical(assert_unit_interval(0), 0)
  expect_identical(assert_unit_interval(0.5), 0.5)
  expect_identical(assert_unit_interval(1), 1)
  expect_error(assert_unit_interval(-0.1), "'-0.1' must be between 0 and 1")
  expect_error(assert_unit_interval(1.1), "'1.1' must be between 0 and 1")
  expect_error(assert_unit_interval(1L, 2), "'1L' must be of length 2")
})

test_that("assert_positive", {
  expect_identical(assert_positive(pi), pi)
  expect_identical(assert_positive(0.1), 0.1)
  expect_error(assert_positive(0), "'0' must be greater than 0")
  expect_error(assert_positive(pi, 2), "'pi' must be of length 2")
})

test_that("assert_nonnegative", {
  expect_identical(assert_nonnegative(pi), pi)
  expect_identical(assert_nonnegative(0), 0)
  expect_error(assert_nonnegative(-1),
               "'-1' must be greater than or equal to 0")
  expect_error(assert_nonnegative(pi, 2), "'pi' must be of length 2")
})

test_that("assert_positive_integer", {
  expect_identical(assert_positive_integer(1L), 1L)
  expect_error(assert_positive_integer(pi), "'pi' must be an integer")
  expect_error(assert_positive_integer(0L), "'0L' must be greater than 0")
  expect_error(assert_positive_integer(1L, 2), "'1L' must be of length 2")
})

test_that("assert_nonnegative_integer", {
  expect_identical(assert_nonnegative_integer(1L), 1L)
  expect_identical(assert_nonnegative_integer(0L), 0L)
  expect_error(assert_nonnegative_integer(pi), "'pi' must be an integer")
  expect_error(assert_nonnegative_integer(-1L),
               "'-1L' must be greater than or equal to 0")
  expect_error(assert_nonnegative_integer(1L, 2), "'1L' must be of length 2")
})
