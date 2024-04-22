test_that("enrichment works", {
  mat <- readRDS(test_path("testdata","mat.Rds"))
  gset <- readRDS(test_path("testdata", "gset.Rds"))
  expect_snapshot(enrich_it(obj = mat, gene.sets = gset, groups = 1000, cores = 2))
})