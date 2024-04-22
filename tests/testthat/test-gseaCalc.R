library(Seurat)
data(seu)
test_that("score calculation works", {
  expect_snapshot_value(gseaCalc(seu), style = "serialize")
})

test_that("different assay throws no error", {
  seu <- RenameAssays(object = seu, Spatial = 'RNA')
  expect_no_error(gseaCalc(seu, assay = 'RNA'))
})

