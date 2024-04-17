seu <- readRDS("data/seu.Rds")
gsea <- readRDS("data/gsea_raw.Rds")
test_that("score calculation works", {
  seu2 <- gseaCalc(seu)
  expect_equal(seu2$gsea_rat_norm, gsea)
})

test_that("different assay works", {
  seu <- RenameAssays(object = seu, Spatial = 'RNA')
  seu2 <- gseaCalc(seu, assay = 'RNA')
  expect_equal(seu2$gsea_rat_norm, gsea)
})

