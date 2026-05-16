test_that("idist agrees with other implementations", {
  expect_equal(
    didist(x = 1, center = 2, scale = 3, fam = "normal"),
    dinvnorm(x = 1, imean = 2, isd = 3)
  )

  expect_equal(
    didist(x = 1, center = 2, scale = 3, fam = "cauchy"),
    didist(x = 1, center = 2, scale = 3, ddist = dt, df = 1)
  )

  expect_equal(
    pidist(q = 1, center = 2, scale = 3, fam = "normal"),
    pinvnorm(q = 1, imean = 2, isd = 3)
  )

  expect_equal(
    pidist(q = 1, center = 2, scale = 3, fam = "cauchy"),
    pidist(q = 1, center = 2, scale = 3, pdist = pt, df = 1)
  )

  expect_equal(
    pidist(q = -1, center = 2, scale = 3, fam = "normal"),
    pinvnorm(q = -1, imean = 2, isd = 3)
  )

  expect_equal(
    pidist(q = -1, center = 2, scale = 3, fam = "cauchy"),
    pidist(q = -1, center = 2, scale = 3, pdist = pt, df = 1)
  )
})
