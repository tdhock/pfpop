library(testthat)
library(pfpop)
test_that("default is not verbose", {
  data_vec <- c(40,50,60)
  penalty <- 30
  result <- pfpop_map_l1(data_vec, penalty)
  expect_identical(result$list_result, NULL)
  expect_identical(result$clusters, NULL)
  expect_identical(result$breaks, NULL)
  expect_equal(nrow(result$iterations), length(data_vec))
  result <- pfpop_list_l1(data_vec, penalty)
  expect_identical(result$model, NULL)
  expect_equal(nrow(result$iterations), length(data_vec))
})

result <- pfpop_map_l1_verbose(data_vec, penalty)
plot(result)

test_that("case 1 and 2: breaks ok", {
  data_vec <- c(40,50,60)
  penalty <- 30
  result <- pfpop_map_verbose(data_vec, penalty)
  gres <- geodesichange::geodesicFPOP_vec(data_vec, penalty, verbose=1)
  plot_check(gres, result)
  ## first push constant from 70 to 10, then grow constant end from 10
  ## to 25. Sign is not correct at 25, why? we grow cluster which
  ## replaces last
  cl1.0 <- result$clusters[data_i==1 & step_i==0]
  expect_equal(cl1.0[, last_param*last_Linear+last_Constant], c(penalty, 0))
  computed.dt <- result$breaks[data_i==1 & step_i==0][order(param)]
  expected.dt <- data.table(data_i=1L, step_i=0L, param=c(10,40,70), Linear_diff=c(-1,2,-1))
  expect_equal(computed.dt, expected.dt)
  computed.dt <- result$breaks[data_i==2 & step_i==0][order(param)]
  expected.dt <- data.table(data_i=2L, step_i=0L, param=c(25,40,50,65), Linear_diff=c(-2,2,2,-2))
  expect_equal(computed.dt, expected.dt)
})

