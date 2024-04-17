mat <- readRDS("/tests/testdata/adj_mat.Rds")
seu <- readRDS("data/seu.Rds")
test_that("adjacency matrix is correct", {
  expect_equal(adj_matrix(seu), matrix(mat, nrow = 112))
})
