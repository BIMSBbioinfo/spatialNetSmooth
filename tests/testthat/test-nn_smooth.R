seu <- readRDS("data/seu.Rds")
gsea <- readRDS("tests/testdata/gsea_nn.Rds")
test_that("smoothing works", {
  seu2 <- nn_smooth(seu)
  expect_equal(seu2$gsea_rat_norm, gsea)
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
