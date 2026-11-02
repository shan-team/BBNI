# Helper to create small valid network and data
make_small_data <- function(n = 3, m = 6, timeseries = TRUE) {
  set.seed(123)
  net <- GenerateNetwork(num.node = n)
  data <- GenerateSample(
    trans_matrix = net,
    SampleSize = m,
    para = rep(0.5, n),
    error = matrix(0, nrow = n, ncol = m),
    timeseries = timeseries
  )
  list(net = net, data = data, n = n, m = m)
}

test_that("run_bbni validates every argument", {
  d <- make_small_data(3, 6)
  expect_error(run_bbni(as.data.frame(d$data)),                       "numeric matrix")
  expect_error(run_bbni(matrix("a", 3, 3)),                            "numeric matrix")
  expect_error(run_bbni({g <- d$data; g[1,1] <- NA; g}),               "missing values")
  expect_error(run_bbni(d$data * 2),                                   "binary")
  expect_error(run_bbni(d$data, num.node = 0),                         "num.node")
  expect_error(run_bbni(d$data, num.node = 2.5),                       "num.node")
  expect_error(run_bbni(d$data, SampleSize = 0),                       "SampleSize")
  expect_error(run_bbni(d$data, num.node = 4),                         "must equal nrow")
  expect_error(run_bbni(d$data, SampleSize = 4),                       "must equal ncol")
  expect_error(run_bbni(d$data, prior_para = matrix(1, 2, 2)),         "prior_para")
  expect_error(run_bbni(d$data, num_update = 0),                       "num_update")
  expect_error(run_bbni(d$data, verbose = "yes"),                      "verbose")
  expect_error(run_bbni(d$data, timeseries = "yes"),                   "timeseries")
  expect_error(run_bbni(d$data, burn_in = 1),                          "burn_in")
  expect_error(run_bbni(d$data, burn_in = -0.1),                       "burn_in")
})

test_that("run_bbni uses flat prior when prior_para is NULL", {
  d <- make_small_data(3, 6)
  res <- run_bbni(GeneData = d$data, num_update = 2, prior_para = NULL)
  expect_true(is.matrix(res$post_edge_prob))
})

test_that("GenerateNetwork validates num.node", {
  expect_error(GenerateNetwork(1),  "at least 2")
  expect_error(GenerateNetwork(2.5), "integer")
  expect_error(GenerateNetwork("a"), "integer")
})

test_that("GenerateSample validates arguments", {
  net <- GenerateNetwork(4)
  expect_error(GenerateSample("a"), "square matrix")
  expect_error(GenerateSample(net[, 1:3]), "square matrix")
  expect_error(GenerateSample(net, SampleSize = 0), "SampleSize")
  expect_error(GenerateSample(net, SampleSize = 1.5), "SampleSize")
  expect_error(GenerateSample(net, num.node = 5), "num.node")
  expect_error(GenerateSample(net, para = rep(2, 4)), "para")
  expect_error(GenerateSample(net, error = matrix(2, 4, 50)), "binary")
  expect_error(GenerateSample(net, timeseries = "x"), "timeseries")
})
