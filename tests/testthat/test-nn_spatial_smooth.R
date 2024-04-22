library(Seurat)
data(seu)

test_that("smoothing works", {
  expect_snapshot_value(nn_spatial_smooth(seu), style = "serialize")
})

test_that("different graph throws no error", {
  expect_no_error(nn_spatial_smooth(seu, graph = 'snn'))
})

test_that("different alpha throws no error", {
  expect_no_error(nn_spatial_smooth(seu, a1 = 0.6 , a2 = 0.7))
})

test_that("different assay works", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(nn_spatial_smooth(seu, assay = 'RNA'))
})
