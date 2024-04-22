data(seu)
test_that("adjacency matrix is correct", {
  expect_snapshot(cat(adj_matrix(seu)))
})
