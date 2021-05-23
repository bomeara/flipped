test_that("flat linear works", {
  nheads <- 2
  nflips <- 3
  nflips <- 3
  expect_equal(dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=0.5, slope=0), stats::dbinom(x=nheads,size=nflips,prob=0.5, log=FALSE), tolerance=0.001)
  nflips <- 3
  expect_equal(dcoin_linear(nheads=nheads, nflips=nflips, preflip_prob=0.5, slope=0, log=TRUE), stats::dbinom(x=nheads,size=nflips,prob=0.5, log=TRUE), tolerance=0.001)
})
