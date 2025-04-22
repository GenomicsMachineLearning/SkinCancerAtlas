# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
output_dir <- '/scratch/project/stseq/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/fig5'
clean_dir <- '/scratch/project/stseq/Feng/SkinCancerAtlas/clean-data'
setwd(output_dir)

# library ----------------
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(hoodscanR)
library(SpatialExperiment)
library(ggpubr)
library(ggrepel)
library(Seurat)
set.seed(111)
source('/scratch/project/stseq/Feng/projects/SkinCancerAtlas/Revision/CosMX/functions.R')

# default path and parameters -----------------------
unique_colour <-  c(
    "DC" = "#5f9d9e",
    "Endothelial Cell" = "#f8a41e",
    "Fibroblast" = "#458b41",
    "KC Basal" = "#f16b6b",
    "KC Cornified" = "#9a1f61",
    "KC Granular" = "#c72685",
    "KC Differentiating" = "#9583bd",
    "KC Hair" = "#eb2627",
    "LC" = "#37479b",
    "Macrophage" = "#eae71d",
    "Melanocytes" = "#8b471f",
    "NK" = "#99ca3e",
    "T Cell" = "#41baeb",
    "Mast Cell" = "#7f8133",
    "B Cell" = "#fed9b9",
    "Pericytes" = "#dca566",
    "Sweat gland related" = "#f2634b",
   'Plasma' = '#f1ea9d',
   'KC IFN' = '#f06ba8', 
   "Monocytes"="#9cc7a1",
   "KC Dysplastic"= "#d8c0dd"
)



# xenium ----------------
xenium_dir <- "/home/uqfzha11/Feng/SkinCancerAtlas/process-data/Xenium_data"
patient_folder <- 'output-XETG00114__0017002__6475-07FC__20231129__074525'
patient <- '6475-07FC'
xenium.obj <- LoadXenium(file.path(xenium_dir,patient_folder), fov = "fov")
all_metadata <- read.csv(file.path(clean_dir,'all_metadata.csv'))
metadata <- all_metadata[all_metadata$patient_cancer==paste0(patient,'_Mel'),]
row.names(metadata) <- metadata$cell_id
xenium.obj <- AddMetaData(xenium.obj,metadata=metadata)
xenium.obj <- xenium.obj[,!is.na(xenium.obj$cell_type)]
xenium.obj$cell_type <- as.factor(xenium.obj$cell_type)
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")

fov_cell_type_p <- plot_cell_boundary(xenium.obj,fov_name='fov',group_by='cell_type',cols=unique_colour,molecules=c("TYR"),
                            is_legend=T,mols.size = 2,mols.cols=c('red'))+
                        coord_flip()+ scale_x_reverse()
ggsave(filename = paste0('xenium_',patient,'_celltype.pdf'),plot=fov_cell_type_p,width = 14993,height = 4299,units='px',dpi=200,limitsize = FALSE)

