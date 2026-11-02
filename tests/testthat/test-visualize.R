# Quick tests for visualization functions functions (without checking actual plots)

# Helper to run a quick MCMC and return results
get_small_results <- function() {
  set.seed(1)
  net <- GenerateNetwork(num.node = 4)
  data <- GenerateSample(net, SampleSize = 8, para = rep(0.5,4),
                         error = matrix(0,4,8))
  run_bbni(GeneData = data, num_update = 2, verbose = FALSE)
}

test_that("plot_bbni runs without error and returns invisible igraph", {
  res <- get_small_results()
  # Use pdf device to suppress actual plot
  pdf(NULL)
  p <- plot_bbni(res, threshold = 0.3)
  dev.off()
  expect_true(inherits(p, "igraph"))
})

test_that("plot_bbni with true_network color coding", {
  res <- get_small_results()
  true_net <- GenerateNetwork(num.node = 4)
  pdf(NULL)
  p <- plot_bbni(res, true_network = true_net, threshold = 0.2)
  dev.off()
  expect_true(inherits(p, "igraph"))
})

test_that("plot_bbni with custom node names", {
  res <- get_small_results()
  pdf(NULL)
  p <- plot_bbni(res, node_names = c("A","B","C","D"))
  dev.off()
  expect_true(inherits(p, "igraph"))
})

test_that("plot_network runs without error", {
  net <- GenerateNetwork(num.node = 4)
  pdf(NULL)
  p <- plot_network(net)
  dev.off()
  expect_true(inherits(p, "igraph"))
})

test_that("plot_trace runs without error (all iterations)", {
  res <- get_small_results()
  pdf(NULL)
  plot_trace(res)
  dev.off()
  expect_null(NULL) # plot_trace returns NULL invisibly; just check no error
})

test_that("plot_trace with thinning runs without error", {
  res <- get_small_results()
  pdf(NULL)
  plot_trace(res, every = 2)
  dev.off()
  expect_null(NULL)
})
