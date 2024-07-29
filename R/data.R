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
#' A Seurat-Object containing spatial and gene-expression data of human breats-tissue with cancer
#' 
#'
#' @format A Seurat-Object
#' @source from 10xGenomics support page, Human Breast Cancer (Block A Section 1) https://www.10xgenomics.com/datasets/
#' data(seu)
"seu"

#' Human breast-cancer-tissue
#'
#' A VoltRon-Object containing spatial and gene-expression data of human breats-tissue with cancer
#' 
#'
#' @format A VoltRon-Object
#' @source a paper by Anderson et al. (2021). The processed count matrices, HE-images and metadata with spot-annotation can be found here: https://zenodo.org/records/4751624. The sample used is B1
#' data(volt)
"volt"
