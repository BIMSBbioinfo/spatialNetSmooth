library(Seurat)
seu <- readRDS(system.file("data", "seu.Rds", package = "spatialNetSmooth"))

test_that("smoothing works", {
  expect_snapshot_value(nn_smooth(seu), style = "serialize")
})

test_that("different graph throws no error", {
  expect_no_error(nn_smooth(seu, graph = 'snn'))
})

test_that("different alpha throws no error", {
  expect_no_error(nn_smooth(seu, a = 0.6))
})

test_that("different assay works", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(nn_smooth(seu, assay = 'RNA'))
})

