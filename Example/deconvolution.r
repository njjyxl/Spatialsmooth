setwd("D:/Bioinfomation_Project/Create_R_Package")

options(warn = -1)
suppressMessages(library(CARD))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratDisk))
suppressMessages(library(spacexr))
suppressMessages(library(Giotto))
suppressMessages(library(SPOTlight))
suppressMessages(library(Redeconve))
suppressMessages(library(AdRoit))
suppressMessages(library(SpatialDecon))
suppressMessages(library(Spatialsmooth))

scrna_path = 'Spatialsmooth/Example/sc_rna.h5seurat'
spatial_path = 'Spatialsmooth/Example/st_data.h5seurat'
celltype_final = 'celltype'
python_path = "D:/Application/anaconda3/envs/python_3.10"
output_path = 'output/test_output/'

# Run all methods
run_deconvolution_methods(scrna_path, spatial_path,celltype_final, python_path, output_path)


