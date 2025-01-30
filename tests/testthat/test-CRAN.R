library(testthat)
library(ggplot2)
library(data.table)
test_that("param is TODO", {
  degrees.vec <- c(10,30,50)
  (result <- pfpop::pfpop(degrees.vec, Inf))
  gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
  map_dt <- melt(
    data.table(result$iterations)[, data_i := .I-1],
    measure.vars=measure(limit, value.name, pattern="(min|max)_(.*)")
  )[
  , N := data_i+1
  ][, let(
    Lmean=Linear/N,
    Cmean=Constant/N,
    cost_mean=cost/N
  )][]
  geodesichange::plot_model(gres$model)+
    geom_point(aes(
      param, cost_mean, color=limit),
      size=4,
      shape=21,
      fill=NA,
      data=map_dt)+
    geom_abline(aes(
      slope=Lmean, intercept=Cmean, color=limit),
      data=map_dt)
  expect_equal(result$segments$param, 1.5)
})
## test_that("param is 0.1", {
##   (result <- pfpop::pfpop(c(0.1,1,359), Inf))
##   expect_equal(result$segments$param, 0.1)
## })
## test_that("param is 359", {
##   (result <- pfpop::pfpop(c(358,1,359), Inf))
##   expect_equal(result$segments$param, 359)
## })
## test_that("param is same as data", {
##   angle.vec <- c(0.1, 6, 6.1, 6.2, 3, 3.1, 2.9)
##   (result <- pfpop::pfpop(angle.vec, 0))
##   expect_equal(result$segments$param, angle.vec)
## })
## test_that("params for reasonable penalties", {
##   angle.vec <- c(0.1, 359, 359.5, 3, 3.1, 2.9)
##   (result <- pfpop::pfpop(angle.vec, 1))
##   expect_equal(result$segments$param, c(359.5,3))
##   (result <- pfpop::pfpop(angle.vec, 0.001))
##   expect_equal(result$segments$param, angle.vec)
##   (result <- pfpop::pfpop(angle.vec, 1000))
##   expect_equal(nrow(result$segments),1)
## })
## test_that("pieces disappear", {
##   fit <- pfpop::pfpop(c(0,180,0,180),Inf)
##   expect_equal(fit$iterations$cost, c(0,90,60,90))
##   expect_equal(fit$iterations$pieces, c(2,1,2,1))
## })
