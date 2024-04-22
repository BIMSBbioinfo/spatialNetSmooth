library(Seurat)
data(seu)

test_that("smoothing works", {
  expect_snapshot_value(union_smooth(seu), style = "serialize")
})
test_that("different graph throws no error", {
  expect_no_error(union_smooth(seu, graph = 'snn'))
})

test_that("different alpha throws no error", {
  expect_no_error(union_smooth(seu, a = 0.6))
})

test_that("different assay throws no error", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(union_smooth(seu, assay = 'RNA'))
})
