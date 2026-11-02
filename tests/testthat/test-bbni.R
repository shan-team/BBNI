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

test_that("run_bbni executes & returns expected structure (timeseries)", {
  d <- make_small_data(3, 8, timeseries = TRUE)
  res <- run_bbni(
    GeneData = d$data,
    num.node = d$n,
    SampleSize = d$m,
    num_update = 2,
    verbose = FALSE
  )
  expect_type(res, "list")
  expect_named(res, c("networks", "log_posterior", "post_edge_prob", "burn_in"))
  expect_true(is.list(res$networks))
  expect_true(is.numeric(res$log_posterior))
  expect_true(is.matrix(res$post_edge_prob))
  expect_equal(dim(res$post_edge_prob), c(d$n, d$n))
  expect_equal(res$burn_in, 0.7)
})

test_that("run_bbni works w/ independent (non-timeseries) data", {
  d <- make_small_data(3, 8, timeseries = FALSE)
  res <- run_bbni(
    GeneData = d$data,
    num_update = 2,
    timeseries = FALSE
  )
  expect_true(is.matrix(res$post_edge_prob))
  expect_equal(dim(res$post_edge_prob), c(3, 3))
})

test_that("run_bbni handles custom priors, penalty, and prop.ratio", {
  d <- make_small_data(4, 6, timeseries = TRUE)
  prior <- matrix(3, nrow = 5, ncol = 2)
  prior[5, ] <- c(2, 100)
  res <- run_bbni(
    GeneData = d$data,
    prior_para = prior,
    penalty = 1,
    prop.ratio = 0.5,
    burn_in = 0.3,
    num_update = 2
  )
  expect_equal(res$burn_in, 0.3)
  expect_true(is.matrix(res$post_edge_prob))
})

test_that("run_bbni prints verbose output w/o error", {
  d <- make_small_data(3, 5, timeseries = TRUE)
  expect_output(
    res <- run_bbni(
      GeneData = d$data,
      num_update = 1,
      verbose = TRUE
    ),
    "Running BBNI MCMC Sampling"
  )
  expect_type(res, "list")
})

test_that("run_bbni validates penalty and prop.ratio", {
  d <- make_small_data(3, 5)
  expect_error(
    run_bbni(GeneData = d$data, penalty = 2),
    "penalty"
  )
  expect_error(
    run_bbni(GeneData = d$data, prop.ratio = -0.1),
    "prop.ratio"
  )
})

test_that("run_bbni retains rownames in post_edge_prob", {
  d <- make_small_data(3, 6)
  rownames(d$data) <- c("G1", "G2", "G3")
  res <- run_bbni(GeneData = d$data, num_update = 1)
  expect_equal(rownames(res$post_edge_prob), c("G1", "G2", "G3"))
  expect_equal(colnames(res$post_edge_prob), c("G1", "G2", "G3"))
})

test_that("run_bbni explores multiple parent configurations", {
  # Create a 4-node network with known structure:
  # Node1: root (no parents)
  # Node2: one parent (node1) – uses rule 11 (identity)
  # Node3: two parents (node1, node2) – uses rule 1 (AND)
  # Node4: one parent (node3) – uses rule 12 (negation)
  net <- matrix(0, nrow = 4, ncol = 4)
  net[2, 1] <- 11
  net[3, 1] <- 1
  net[3, 2] <- 1
  net[4, 3] <- 12

  set.seed(2024)
  data <- GenerateSample(net, SampleSize = 10, para = rep(0.5, 4),
                         error = matrix(0, 4, 10), timeseries = TRUE)

  # Run enough iterations to allow proposals of add/remove/swap moves
  set.seed(99)  # deterministic chain
  res <- run_bbni(GeneData = data,
                  num.node = 4,
                  SampleSize = 10,
                  num_update = 20, # more iterations to explore branches
                  verbose = FALSE)

  # Basic checks
  expect_true(is.matrix(res$post_edge_prob))
  expect_equal(dim(res$post_edge_prob), c(4, 4))
  expect_true(all(res$post_edge_prob >= 0 & res$post_edge_prob <= 1))

  # The chain should have visited at least some states with 1 and 2 parents
  # (can check by looking at the number of edges in sampled networks)
  edge_counts <- sapply(res$networks, function(m) sum(m > 0))
  expect_true(any(edge_counts == 1))
  expect_true(any(edge_counts == 2))
  expect_true(any(edge_counts == 3))
})

test_that("run_bbni w/ prop.ratio = 0.9 covers empirical proposal branches", {
  # A small 3-node network; multiple runs with different seeds to explore many branches
  for (seed in c(1, 2, 3, 4, 5)) {
    set.seed(seed)
    n <- 3
    m <- 8
    net <- GenerateNetwork(num.node = n)
    data <- GenerateSample(net, SampleSize = m, para = rep(0.5, n),
                           error = matrix(0, n, m), timeseries = TRUE)
    # Use prop.ratio = 0.9 for frequent use of the empirical proposal distribution
    res <- run_bbni(GeneData = data,
                    num_update = 15, # enough iterations to move through many states
                    prop.ratio = 0.9,
                    verbose = FALSE)
    # Basic checks (no need to assert coverage; just ensure it runs)
    expect_true(is.matrix(res$post_edge_prob))
  }
})

test_that("run_bbni on 3-node chain hits limited-parent edge cases", {
  # Chain: 1 -> 2 -> 3. This creates states where when updating node 2,
  # swap_candi (legal parents) has length 1, and when updating node 3,
  # it has length 0, covering many previously uncovered branches.
  n <- 3
  m <- 8
  net <- matrix(0, n, n)
  net[2, 1] <- 11   # node2 = identity(node1)
  net[3, 2] <- 11   # node3 = identity(node2)

  for (seed in 10:14) {
    set.seed(seed)
    data <- GenerateSample(net, SampleSize = m, para = rep(0.5, n),
                           error = matrix(0, n, m), timeseries = TRUE)
    res <- run_bbni(GeneData = data,
                    num_update = 15, # enough to explore states
                    prop.ratio = 0.9, # use empirical proposals often
                    verbose = FALSE)
    expect_true(is.matrix(res$post_edge_prob))
  }
})

test_that("run_bbni w/ independent data and prop.ratio = 0.9 covers non-timeseries paths", {
  # Independent (cross-sectional) mode with empirical proposals
  for (seed in 20:22) {
    set.seed(seed)
    n <- 3
    m <- 8
    net <- GenerateNetwork(num.node = n)
    data <- GenerateSample(net, SampleSize = m, para = rep(0.5, n),
                           error = matrix(0, n, m), timeseries = FALSE)
    res <- run_bbni(GeneData = data,
                    timeseries = FALSE,
                    num_update = 12,
                    prop.ratio = 0.9,
                    verbose = FALSE)
    expect_true(is.matrix(res$post_edge_prob))
  }
})

test_that("run_bbni w/ a network of two-parent node executes swap moves", {
  # Force initial state with a node having two parents, covering many two-parent case branches
  set.seed(42)
  n <- 4
  m <- 10
  # Create a network: node1 root, node2 root, node3 = node1 AND node2, node4 = node3
  net <- matrix(0, n, n)
  net[3, 1] <- 1
  net[3, 2] <- 1
  net[4, 3] <- 11
  data <- GenerateSample(net, SampleSize = m, para = rep(0.5, n),
                         error = matrix(0, n, m), timeseries = TRUE)
  # Run with high proposal usage and multiple seeds to cover various swap/remove branches
  for (seed in c(100, 101, 102)) {
    set.seed(seed)
    res <- run_bbni(GeneData = data,
                    num_update = 15,
                    prop.ratio = 0.9,
                    verbose = FALSE)
    expect_true(is.matrix(res$post_edge_prob))
  }
})


