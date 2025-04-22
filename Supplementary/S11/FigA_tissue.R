# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
raw_dir <- '/home/uqfzha11/Feng/SkinCancerAtlas/process-data/flat_file_format_raw_data/'
output_dir <- '/scratch/project/stseq/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/figS10'
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


# cosmx ----------------
samples <- readRDS(file.path(clean_dir, "cosmx_spatial_all.rds"))

## B18_SCC_15 ----------------------
# cell type
patient = "B18_SCC"
print(patient)
fov_id=15
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples@meta.data[samples$orig.ident==patient,]
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

# cell type
fov_cell_type_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='predicted.id',
                            cols=unique_colour,is_legend=F)+
                        coord_flip()
ggsave(filename = paste0(patient,'_',fov_id,'_celltype.pdf'),plot=fov_cell_type_p,width = 5472,height = 3648,units='px',dpi=200)
