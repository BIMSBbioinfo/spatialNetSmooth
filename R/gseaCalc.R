#' calculate gsea-score
#' @param se    Seurat object
#' @param genes csv-file indicating which genes to use (optional)
#' @param assay   which assay to use (optional)
#' @return Seurat object with calculated scores in meta.data
#' @export
#' @import Seurat
#' @import dplyr
#' @importFrom GSEABase GeneSetCollection GeneSet
#' @importFrom magrittr set_rownames


gseaCalc <- function(se, genes = "Datasets - Ikarus - Gene_lists.csv", assay= "Spatial"){
  if(genes == "Datasets - Ikarus - Gene_lists.csv"){
  genes <- system.file("data", "Datasets - Ikarus - Gene_lists.csv", package = "spatialNetSmooth")
  }else{
    genes <- read.csv(genes)
  }
  genes = read.csv(genes) %>%
    dplyr::rename(type = X)
  glist = with(genes, split(genes, type))
  ndims = 10
  VariableFeatures(se, assay = assay) = unlist(glist)
  
  se = se %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures(selection.method = 'mvp') %>%
    RunPCA(verbose = FALSE, features = unlist(glist)) %>%
    FindNeighbors(reduction = "pca", dims = 1:ndims) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:ndims)

  lnams = names(glist)
  gset = GeneSetCollection(lapply(setNames(lnams, lnams), function(x)GeneSet(glist[[x]], setName=x)))
  mat = as.matrix(GetAssayData(se, assay))
  se.it <- enrich_it(obj = mat, gene.sets = gset, groups = 1000, cores = 2)
  rm(mat)

  seu = se
  seu@meta.data = seu@meta.data %>%
    cbind(
      se.it %>%
        as.data.frame() %>%
        magrittr::set_colnames(paste0('gsea_', colnames(.)))
    ) %>%
    mutate(gsea_Normal_norm      = gsea_Normal) %>%
    mutate(gsea_Normal_norm      = gsea_Normal_norm - min(gsea_Normal_norm)) %>%
    mutate(gsea_Normal_norm      = gsea_Normal_norm / max(gsea_Normal_norm)) %>%

    mutate(gsea_Tumor_norm       = gsea_Tumor) %>%
    mutate(gsea_Tumor_norm       = gsea_Tumor_norm  - min(gsea_Tumor_norm)) %>%
    mutate(gsea_Tumor_norm       = gsea_Tumor_norm / max(gsea_Tumor_norm)) %>%

    mutate(gsea_rat              = log2((gsea_Tumor + 1)/ (gsea_Normal + 1))) %>%
    mutate(gsea_rat_norm = log2((gsea_Tumor_norm+1) / (gsea_Normal_norm + 1))) %>%
    magrittr::set_rownames(.$cell_id)
  return(seu)
}
