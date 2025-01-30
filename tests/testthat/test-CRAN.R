library(testthat)
library(ggplot2)
library(data.table)

test_that("argmin in 3rd step is 350", {
  degrees.vec <- c(10,310, 350)
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
  expect_equal(result$iterations$min_param, c(10,340,350))
})

test_that("argmin in 2nd step is 340", {
  degrees.vec <- c(310,10)
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
  expect_equal(result$iterations$min_param, c(310,340))
})

## works.
degrees.vec <- c(310,10, 300, 320, 20, 30, 40, 50)
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

## two of same data does not work for other solver?
degrees.vec <- c(310,10, 300, 310)
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

test_that("map size goes to 1 then 0", {
  degrees.vec <- c(180, 0, 180, 0)
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
  expect_equal(result$iterations$map_size, c(1, 0, 1, 0))
})
