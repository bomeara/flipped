test_that("flat linear works", {
  nheads <- 2
  nflips <- 3
  expect_equal(dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=0.5, slope=0), stats::dbinom(x=nheads,size=nflips,prob=0.5, log=FALSE), tolerance=0.001)
  expect_equal(dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=0.5, slope=0, log=TRUE), stats::dbinom(x=nheads,size=nflips,prob=0.5, log=TRUE), tolerance=0.001)
})


test_that("exponential works", {
  nheads <- 3
  nflips <- 3
  expect_equal(dcoin_exponential(nheads=nheads, nflips=nflips, halflife=Inf), stats::dbinom(x=nheads,size=nflips,prob=1, log=FALSE), tolerance=0.001)
  expect_equal(dcoin_exponential(nheads=nheads, nflips=nflips, halflife=Inf, log=TRUE), stats::dbinom(x=nheads,size=nflips,prob=1, log=TRUE), tolerance=0.001)
})

test_that("profiling works", {
  nheads <- 7
  nflips <- 10
  result <- profile_linear_model(nheads, nflips, slope=0.1)
  expect_s3_class(result, "data.frame")
  
  result <- profile_exponential_model(nheads, nflips)
  expect_s3_class(result, "data.frame")
})