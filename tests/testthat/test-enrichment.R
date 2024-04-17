test_that("enrichment works", {
  mat <- readRDS("tests/testdata/mat.Rds")
  gset <- readRDS("tests/testdata/gset.Rds")
  enriched <- readRDS("tests/testdata/enriched.Rds")
  se.it <- enrich_it(obj = mat, gene.sets = gset, groups = 1000, cores = 2)
  expect_equal(se.it, enriched)
})