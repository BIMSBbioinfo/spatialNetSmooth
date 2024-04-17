test_that("enrichment works", {
  genes = read.csv(genes) %>%
    dplyr::rename(type = X)
  glist = with(genes, split(genes, type))
  ndims = 10
  
  VariableFeatures(se) = unlist(glist)
  
  
  se = se %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures(selection.method = 'mvp') %>%
    RunPCA(verbose = FALSE, features = unlist(glist)) %>%
    FindNeighbors(reduction = "pca", dims = 1:ndims) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:ndims)
  seu2 <- nn_spatial_smooth(seu)
  lnams = names(glist)
  gset = GeneSetCollection(lapply(setNames(lnams, lnams), function(x)GeneSet(glist[[x]], setName=x)))
  mat = as.matrix(GetAssayData(se, assay))
  se.it <- enrich_it(obj = mat, gene.sets = gset, groups = 1000, cores = 2)
  expect_equal(seu2$gsea_rat_norm, gsea)
})