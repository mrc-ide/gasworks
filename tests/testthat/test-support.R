
test_that("model_day functions work", {

  expect_equal(model_day("2014-01-01"), 1)
  expect_equal(model_day("2014-01-07"), 7)
  expect_equal(model_day("2013-12-30"), -1)
  expect_equal(model_day(as.Date("2013-12-30")), -1)

  expect_equal(model_date(1), as.Date("2014-01-01"))
  expect_equal(model_date(7), as.Date("2014-01-07"))
  expect_equal(model_date(-1), as.Date("2013-12-30"))
})

test_that("model_week functions work", {

  expect_equal(model_week("2014-01-01"), 1)
  expect_equal(model_week("2014-01-07"), 1)
  expect_equal(model_week("2014-01-08"), 2)
  expect_equal(model_week(as.Date("2013-12-30")), 0)

  expect_equal(model_week_date(1), as.Date("2014-01-07"))
  expect_equal(model_week_date(2), as.Date("2014-01-14"))
  expect_equal(model_week_date(0), as.Date("2013-12-31"))
})
