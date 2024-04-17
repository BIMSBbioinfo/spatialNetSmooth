seu <- readRDS(system.file("data", "seu.Rds", package = "spatialNetSmooth"))

test_that("smoothing works", {
  expect_snapshot(cat(alpha_nn_spatial_smooth(seu)))
})

test_that("different graph throws no error", {
  expect_no_error(alpha_nn_spatial_smooth(seu, graph = 'snn'))
})

test_that("different alpha throws no error", {
  expect_no_error(alpha_nn_spatial_smooth(seu, a = 0.6 , alpha = 0.7))
})

test_that("different assay works", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(alpha_nn_spatial_smooth(seu, assay = 'RNA'))
})
