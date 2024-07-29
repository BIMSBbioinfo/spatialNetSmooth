spatialNetSmooth
========

spatialNetSmooth is an R package for smoothing of GSEA-scores from spatial omics breast cancer data. The smoothing is based on network-propagation and uses the [netSmooth](https://github.com/BIMSBbioinfo/netSmooth) package.
It contains functions that work with different methods based on spatial and/or profile-neighborhood. The input can be either a [Seurat](https://github.com/satijalab/seurat) or a [VoltRon](https://github.com/BIMSBbioinfo/VoltRon) object.

## Current Features
* Calculating GSEA-scores using gene lists from [ikarus](https://doi.org/10.1186/s13059-022-02683-1) or custom list
* spatial smoothing
* NN/SNN smoothing
* spatial+NN/SNN smoothing
* union/intersection/linear combination of NN/SNN and spatial smoothing
* Calculating F1-scores of different thresholds based on quantiles
* Calculating best thresholds based on ROC-curve
* Spatial plotting of prediction compared to known ground truth (thresholded with ROC-curve)
  
## Installation

```r
library(devtools)
install_github("BIMSBbioinfo/spatialNetSmooth", build_vignettes=FALSE, 
  dependencies=TRUE)
```

# How to use
See the [Tutorial](https://github.com/BIMSBbioinfo/spatialNetSmooth/blob/main/vignettes/Tutorial.Rmd)

