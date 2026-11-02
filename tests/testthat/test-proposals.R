test_that("ProposalConstruction returns valid candidate matrices", {
  set.seed(42)
  n <- 3
  m <- 6
  # Generate simple binary data
  GeneData <- matrix(rbinom(n * m, 1, 0.5), nrow = n, ncol = m)
  res <- BBNI:::ProposalConstruction(GeneData, SampleSize = m, timeseries = TRUE)
  expect_type(res, "list")
  expect_length(res, 2)
  # triplet matrix
  expect_true(is.matrix(res[[1]]))
  expect_gt(ncol(res[[1]]), 0)
  # pairwise matrix
  expect_true(is.matrix(res[[2]]))
  expect_gt(ncol(res[[2]]), 0)

  # check that triplet has at least 6 cols (i,j,k,code,weight,raw)
  if (nrow(res[[1]]) > 0) {
    expect_equal(ncol(res[[1]]), 6)
  }
  # check that pairwise has at least 5 cols (i,k,code,weight,raw)
  if (nrow(res[[2]]) > 0) {
    expect_equal(ncol(res[[2]]), 5)
  }
})

test_that("ProposalConstruction works in independent mode", {
  set.seed(7)
  n <- 3
  m <- 5
  GeneData <- matrix(sample(0:1, n * m, replace = TRUE), nrow = n)
  res <- BBNI:::ProposalConstruction(GeneData, SampleSize = m, timeseries = FALSE)
  expect_type(res, "list")
  expect_true(is.matrix(res[[1]]))
  expect_true(is.matrix(res[[2]]))
})

test_that("ConstructInitial returns a valid DAG transition matrix", {
  set.seed(123)
  n <- 4
  m <- 6
  GeneData <- matrix(sample(0:1, n * m, replace = TRUE), nrow = n)
  cand <- BBNI:::ProposalConstruction(GeneData, SampleSize = m, timeseries = TRUE)
  init <- BBNI:::ConstructInitial(cand, num.node = n)
  expect_true(is.matrix(init))
  expect_equal(dim(init), c(n, n))
  # all values must be 0 or positive integers (1..12)
  expect_true(all(init == 0 | (init >= 1 & init <= 12)))
  # convert to incidence matrix and check DAG
  incid <- (init > 0) * 1
  ances <- BBNI:::update_ancestor_matrix(incid)
  expect_equal(sum(diag(ances)), 0)
  # max in-degree <= 2
  expect_true(all(rowSums(incid) <= 2))
})
