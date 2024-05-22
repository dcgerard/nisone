test_that("1F1 works", {
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 1, method = "integral"), exp(1) - 1)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "integral"), (exp(2) - 1) / 2)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "integral"), (exp(2) - 1) / 2)
  expect_equal(hypergeo1f1(a = 2, b = 3, z = 1, method = "integral"), 2)
  expect_equal(hypergeo1f1(a = 3, b = 5, z = 3, method = "integral"), (exp(3) - 4) * 4 / 9)

  expect_equal(hypergeo1f1(a = 1, b = 2, z = 1, method = "series"), exp(1) - 1, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "series"), (exp(2) - 1) / 2, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 1, b = 2, z = 2, method = "series"), (exp(2) - 1) / 2, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 3, b = 5, z = 3, method = "series"), (exp(3) - 4) * 4 / 9, tolerance = 1e-3)
  expect_equal(hypergeo1f1(a = 4, b = 3, z = 1, method = "series"), exp(1) * 4 / 3, tolerance = 1e-3)
})
