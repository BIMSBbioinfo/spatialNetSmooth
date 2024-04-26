library(Seurat)
data(seu)

test_that("smoothing works", {
  expect_snapshot_value(spatial_smooth(seu), style = "serialize")
})

test_that("different alpha throws no error", {
  expect_no_error(spatial_smooth(seu, a = 0.6))
})

test_that("different assay works", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(spatial_smooth(seu, assay = 'RNA'))
})
