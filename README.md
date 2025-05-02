# Spatialsmooth  
*A Spatially-Aware Convolutional Autoencoder Framework for Enhanced Deconvolution of Spatial Transcriptomics Data*  
![](https://github.com/njjyxl/Spatialsmooth/blob/devel/Spatialsmooth.jpg)  
CAE for spatial smoothing aims at smoothing the cell type matrix obtained by various types of deconvolution methods and is used to ensure that this spatial distribution is smooth and continuous, meaning that changes in cell type do not suddenly and abruptly change spatially, but rather gradually. This is important for both the fidelity and reliability of biological data, as biological tissues are usually structured and continuous. CAE requires a deconvolutionally inferred cell type composition for each spatial location as well as localization information for spatial transcriptomics data. With these two inputs, the positional encoding is integrated with the matrix for convolutional self-encoding via CAE. 
## Installation  
You can install Spatialsmooth on Github with the following code.  
## Dependencies  
+ R version >= 4.2.0.  
+ R packages: tidyverse,spacexr,Matrix,data.table,Seurat(>=4.4.0),SeuratDisk,CARD,SPOTlight,Giotto,Redeconve,AdRoit,parallel, doParallel,dplyr,future, spatialreg,sp,igraph,fields,keras,tensorflow,reticulate,foreach.  
```  
# install devtools if necessary  
install.packages('devtools')  

# install the Spatialsmooth package  
devtools::install_github('njjyxl/Spatialsmooth')  

# load package  
library(Spatialsmooth)  
```
The R package has been installed successfully on Operating systems:  
+ CentOS Linux release 7.5.1804 (Core)  
+ Windows 11  

# Issues  
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R!  
# How to use Spatialsmooth  
Detail in Example file


