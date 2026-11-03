# Helper to run a quick MCMC and return results
get_small_results <- function() {
  set.seed(1)
  net <- GenerateNetwork(num.node = 4)
  data <- GenerateSample(net, SampleSize = 8, para = rep(0.5,4),
                         error = matrix(0,4,8))
  run_bbni(GeneData = data, num_update = 2, verbose = FALSE)
}

test_that("summary.bbni returns expected structure and prints", {
  res <- get_small_results()
  out <- summary(res, threshold = 0.5, n_top = 3)
  expect_named(out, c("num_nodes","num_update","burn_in",
                      "final_log_posterior","n_strong_edges","top_edges"))
  expect_equal(out$num_nodes, 4)
  expect_equal(nrow(out$top_edges), 3)
})

test_that("summary.bbni handles rownames NULL", {
  res <- get_small_results()
  rownames(res$post_edge_prob) <- NULL
  expect_output(summary(res), "BBNI MCMC summary")
})

test_that("summary.bbni and print.bbni fall back when attr(num_update) is NULL", {
  res <- get_small_results()
  attr(res, "num_update") <- NULL
  expect_output(summary(res), "BBNI MCMC summary")
  expect_output(print(res), "<bbni>")
})

test_that("plot.bbni delegates to plot_trace", {
  res <- get_small_results()
  pdf(NULL)
  expect_invisible(plot(res))
  dev.off()
})

test_that("plot.bbni returns NULL when burn_in is absent", {
  res <- get_small_results()
  res$burn_in <- NULL
  pdf(NULL)
  expect_null(plot(res))
  dev.off()
})
