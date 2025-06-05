library(testthat)
devtools::load_all("~/R/pfpop")
library(pfpop)
test_that("default is not verbose", {
  data_vec <- c(355,365,375)
  penalty <- Inf
  result <- pfpop_list_l1(data_vec, penalty)
  expect_identical(result$iterations$best_param[3], 365)
  result <- pfpop_map_l1(data_vec, penalty)
  expect_identical(result$list_result, NULL)
  expect_identical(result$clusters, NULL)
  expect_identical(result$breaks, NULL)
  expect_equal(nrow(result$iterations), length(data_vec))
  result <- pfpop_list_l1(data_vec, penalty)
  expect_identical(result$model, NULL)
  expect_equal(nrow(result$iterations), length(data_vec))
})
test_that("verbose ok", {
  data_vec <- c(40,50,60)
  penalty <- Inf
  result <- pfpop_map_l1_verbose(data_vec, penalty)
  plot(result)
  model0 <- result$list_result$model[data_i==0]
  expect_equal(model0$min_param, 40)
  expect_equal(model0$max_param, 60)
})
