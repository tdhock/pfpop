library(testthat)
library(ggplot2)
library(data.table)

mean_cost <- function(result)melt(
  data.table(result$iterations)[, data_i := .I-1],
  measure.vars=measure(limit, value.name, pattern="(min|max)_(.*)")
)[
, N := data_i+1
][, let(
  Lmean=Linear/N,
  Cmean=Constant/N,
  cost_mean=cost/N
)][]
plot_check <- function(gres, result){
  map_dt <- mean_cost(result)
  geodesichange::plot_model(gres$model)+
    theme_bw()+
    theme(panel.margin=grid::unit(0,"lines"))+
    facet_grid(data_i ~ ., scales="free")+
    geom_point(aes(
      param, cost_mean, color=limit),
      size=4,
      shape=21,
      fill=NA,
      data=map_dt)+
    geom_abline(aes(
      slope=Lmean, intercept=Cmean, color=limit),
      data=map_dt)
}

test_that("argmin in 3rd step is 350", {
  degrees.vec <- c(10,310, 350)
  (result <- pfpop::pfpop_map(degrees.vec, Inf))
  gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
  plot_check(gres, result)
  expect_equal(result$iterations$min_param, c(10,340,350))
})

test_that("argmin in 2nd step is 340", {
  degrees.vec <- c(310,10)
  (result <- pfpop::pfpop_map(degrees.vec, Inf))
  gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
  plot_check(gres, result)
  expect_equal(result$iterations$min_param, c(310,340))
})

## works.
degrees.vec <- c(310,10, 300, 320, 20, 30, 40, 50)
(result <- pfpop::pfpop_map(degrees.vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
plot_check(gres, result)

## two of same data does not work for other solver?
degrees.vec <- c(310,10, 300, 310)
(result <- pfpop::pfpop_map(degrees.vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
plot_check(gres, result)

test_that("map size goes to 1 then 0", {
  degrees.vec <- c(180, 0, 180, 0)
  (result <- pfpop::pfpop_map(degrees.vec, Inf))
  gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
  plot_check(gres, result)
  expect_equal(result$iterations$map_size, c(1, 0, 1, 0))
})

test_that("map size goes to 3 then 0", {
  degrees.vec <- c(180, 90, 0, 270)
  (result <- pfpop::pfpop_map(degrees.vec, Inf))
  gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
  plot_check(gres, result)
  expect_equal(result$iterations$map_size, c(1, 3, 2, 0))
  expect_equal(result$iterations$min_cost, c(0, 90, 180, 360))
})

set.seed(1)
data_vec <- runif(8,0,360)
(result <- pfpop::pfpop_map(data_vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
plot_check(gres, result)

## Multiple global min/max.
N <- 9
data_vec <- seq(1, 90, l=N)+rep(c(0,180),l=N)
(result <- pfpop::pfpop_map(data_vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
plot_check(gres, result)

## This seems like a worst case number of pointer moves.
N <- 5
data_vec <- (seq(0, 1, l=N)^2+1)*90/2+rep(c(0,180),l=N)
(result <- pfpop::pfpop_map(data_vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
plot_check(gres, result)

data_vec <- c(10, 20+180,40,50+180,80,90+180)
(result <- pfpop::pfpop_map(data_vec, Inf))
gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
plot_check(gres, result)

