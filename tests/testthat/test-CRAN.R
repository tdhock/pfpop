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

map_list_atime <- atime::atime(
  setup={
    set.seed(1)
    data_vec <- runif(N,0,360)
  },
  map={
    mres <- pfpop::pfpop_map(data_vec, Inf)
    with(mres$iterations, data.frame(
      cost=min_cost[N]/N,
      space=max(map_size),
      moves=sum(pointer_moves)))
  },
  list={
    lres <- pfpop::pfpop_list(data_vec, Inf)
    data.frame(
      cost=lres$iterations$best_cost[N],
      space=max(lres$iterations$num_pieces),
      moves=NA)
  },
  seconds.limit=0.1,
  result=TRUE
)
plot(map_list_atime)

library(data.table)
dcast(
  map_list_atime$measurements,
  N ~ expr.name,
  value.var="cost"
)[list<map]

map_list_refs <- atime::references_best(map_list_atime)
plot(map_list_refs)
