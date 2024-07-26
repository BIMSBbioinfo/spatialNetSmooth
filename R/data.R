#' List of human genes
#'
#' Genes where high expression is 
#' either associated with normal cells or cancer cells
#'obtained using Ikarus
#'https://github.com/BIMSBbioinfo/ikarus
#'
#' @format A list of genes with with Type "Tumor" or "Normal"
#' @source ikarus machine learning pipeline (Dohmen et al. 2022)
#' data(geneset_ikarus)
"geneset_ikarus"

#' Human breast-cancer-tissue
#'
#' An Seurat-Object containing spatial and gene-expression data of human breats-tissue with cancer
#' 
#'
#' @format A Seurat-Object
#' @source from 10xGenomics support page, Human Breast Cancer (Block A Section 1) https://www.10xgenomics.com/datasets/
#' data(seu)
"seu"

