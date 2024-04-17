mat <- readRDS("/tests/testdata/adj_mat.Rds")
test_that("adjacency matrix is correct", {
  expect_equal(matrix(mat, nrow = 112))
})
