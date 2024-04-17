seu <- readRDS("data/seu.Rds")
gsea <- readRDS("data/gsea_spatial.Rds")
test_that("smoothing works", {
  seu2 <- spatial_smooth(seu)
  expect_equal(seu2$gsea_rat_norm, gsea)
})

test_that("different graph throws no error", {
  expect_no_error(spatial_smooth(seu, graph = 'snn'))
})

test_that("different alpha throws no error", {
  expect_no_error(spatial_smooth(seu, a = 0.6))
})

test_that("different assay works", {
  seu <-RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(spatial_smooth(seu, assay = 'RNA'))
})
