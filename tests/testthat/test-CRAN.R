library(testthat)
library(ggplot2)
library(data.table)

pfpop_map_verbose <- function(degrees_vec, penalty=Inf, weight_vec = rep(1, length(degrees_vec)), verbose_file=tempfile()){
  result <- pfpop::pfpop_map(degrees_vec, penalty, weight_vec, verbose_file)
  result$clusters <- fread(verbose_file)
  result
}

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
cldt <- function(data_i, opt, start, end){
  sedt <- data.table(
    data_i,
    start=as.numeric(start),
    end=as.numeric(end)
  )[
    ##start != end
  ][
  , noInf := start <= end
  ]
  data.table(opt, rbind(
    sedt[noInf==TRUE],
    sedt[noInf==FALSE][, let(start=-Inf)],
    sedt[noInf==FALSE][, let(end=Inf)]))
}
plot_check <- function(gres, result){
  map_dt <- mean_cost(result)
  cluster.dt <- result$clusters[, rbind(
    cldt(data_i, "before", first_param, opt_param),
    cldt(data_i, "after", opt_param, last_param))]
  lab.dt <- melt(
    result$clusters,
    measure.vars=measure(
      pointer, value.name, pattern="(first|opt|last)_(param|diff)"))
  geodesichange::plot_model(gres$model)+
    theme_bw()+
    theme(panel.margin=grid::unit(0,"lines"))+
    facet_grid(data_i ~ ., scales="free")+
    geom_rect(aes(
      xmin=start, xmax=end,
      fill=opt,
      ymin=-Inf, ymax=Inf),
      data=cluster.dt,
      alpha=0.5,
      color="black")+
    geom_label(aes(
      param, ifelse(diff<0, -Inf, Inf),
      vjust=ifelse(diff<0, 0, 1),
      label=diff),
      data=lab.dt,
      alpha=0.5)+
    scale_fill_manual(values=c(
      before="grey50",
      after="violet"))+
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

## test_that("argmin in 3rd step is 350", {
##   degrees.vec <- c(10, 310, 350)
##   (result <- pfpop_map_verbose(degrees.vec))
##   gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
##   plot_check(gres, result)
##   expect_equal(result$iterations$min_param, c(10,340,350))
## })

## test_that("argmin in 2nd step is 340", {
##   degrees.vec <- c(310,10)
##   (result <- pfpop::pfpop_map(degrees.vec, Inf))
##   gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
##   plot_check(gres, result)
##   expect_equal(result$iterations$min_param, c(310,340))
## })

## ## works.
## degrees.vec <- c(310,10, 300, 320, 20, 30, 40, 50)
## (result <- pfpop::pfpop_map(degrees.vec, Inf))
## gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
## plot_check(gres, result)

## ## two of same data does not work for other solver?
## degrees.vec <- c(310,10, 300, 310)
## (result <- pfpop::pfpop_map(degrees.vec, Inf))
## gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
## plot_check(gres, result)

## test_that("map size goes to 1 then 0", {
##   degrees.vec <- c(180, 0, 180, 0)
##   (result <- pfpop::pfpop_map(degrees.vec, Inf))
##   gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
##   plot_check(gres, result)
##   expect_equal(result$iterations$map_size, c(1, 0, 1, 0))
## })

## test_that("map size goes to 3 then 0", {
##   degrees.vec <- c(180, 90, 0, 270)
##   (result <- pfpop::pfpop_map(degrees.vec, Inf))
##   gres <- geodesichange::geodesicFPOP_vec(degrees.vec, Inf, verbose=1)
##   plot_check(gres, result)
##   expect_equal(result$iterations$map_size, c(1, 3, 2, 0))
##   expect_equal(result$iterations$min_cost, c(0, 90, 180, 360))
## })

## set.seed(1)
## data_vec <- runif(8,0,360)
## (result <- pfpop::pfpop_map(data_vec, Inf))
## gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
## plot_check(gres, result)

## ## Multiple global min/max.
## N <- 9
## data_vec <- seq(1, 90, l=N)+rep(c(0,180),l=N)
## (result <- pfpop::pfpop_map(data_vec, Inf))
## gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
## plot_check(gres, result)

## ## This seems like a worst case number of pointer moves.
## N <- 5
## data_vec <- (seq(0, 1, l=N)^2+1)*90/2+rep(c(0,180),l=N)
## (result <- pfpop::pfpop_map(data_vec, Inf))
## gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
## plot_check(gres, result)

## data_vec <- c(10, 20+180,40,50+180,80,90+180)
## (result <- pfpop_map_verbose(data_vec))
## gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
## plot_check(gres, result)

data_vec <- c(30, 20, 335, 10, 325, 340, 330, 320, 310, 300)
data_vec <- c(310, 300, 50, 60, 70, 185, 190, 200)
(result <- pfpop_map_verbose(data_vec))
gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
plot_check(gres, result)

test_that("clusters split", {
  data_vec <- c(10, 200, 40)
  (result <- pfpop_map_verbose(data_vec))
  if(interactive()){
    gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
    plot_check(gres, result)
  }
  sign_data_counts <- result$clusters[, .(clusters=.N), keyby=.(data_i,sign)]
  expect_equal(sign_data_counts$clusters, c(1,1,1,1,3,3))
  cl2 <- result$clusters[data_i==2][order(first_param)]
  expected_param <- c(10,20,40,190,200,220)
  expect_equal(cl2$first_param, expected_param)
  expect_equal(cl2$opt_param, expected_param)
  expect_equal(cl2$last_param, expected_param)
  expect_equal(cl2$sign, c(1,-1,1,-1,1,-1))
})

gdist <- function(point, data_vec){
  abs_vec <- abs(point-data_vec)
  d_vec <- pmin(abs_vec, 360-abs_vec)
  sum(d_vec)
}

test_that("move max from 370 to 10 has correct cost", {
  data_vec <- c(170, 200, 190)
  (result <- pfpop_map_verbose(data_vec))
  if(interactive()){
    gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
    plot_check(gres, result)
  }
  expect_equal(result$iterations$max_cost, c(180, 330, 510))
})

test_that("move min from 355 to 25 has correct cost", {
  data_vec <- c(355, 25, 35)
  (result <- pfpop_map_verbose(data_vec))
  if(interactive()){
    gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
    plot_check(gres, result)
  }
  expect_equal(result$iterations$min_cost, c(0, 30, 40))
})
