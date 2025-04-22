# library --------------
rm(list=ls());options(stringsAsFactors=FALSE)
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)

# set path ------------------
input_dir <- '/scratch/project/stseq/Feng/SkinCancerAtlas/clean-data'
output_dir <- '/home/uqfzha11/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/figS8'
raw_dir <- "/home/uqfzha11/Feng/SkinCancerAtlas/process-data/Xenium_data"
setwd(output_dir)
set.seed(111)
source('/scratch/project/stseq/Feng/projects/SkinCancerAtlas/Revision/CosMX/functions.R')

# cell chat analysis ----------------------
all_metadata <- read.csv(file.path(input_dir,'all_metadata.csv'), check.names = FALSE)
patient_info <- list('23346-10SP'=list(patient_folder='output-XETG00114__0017002__23346-10SP__20231129__074525',
                                       height=5308,width=17879),
                     '30037-07BR'=list(patient_folder='output-XETG00114__0017002__30037-07BR__20231129__074525',
                                       height=8079,width=19827),
                     '6475-07FC'=list(patient_folder='output-XETG00114__0017002__6475-07FC__20231129__074525',
                                      height = 5308,width = 17879))
community_colour <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#6A3D9A", "#FDBF6F", "#FF7F00", "#CAB2D6", "#E31A1C")
names(community_colour) <- as.character(0:9)

# 23346-10SP
for(patient in names(patient_info)){
  # patient <- '6475-07FC'
  metadata <- all_metadata[all_metadata$batch==patient,]
  row.names(metadata) <- metadata$cell_id
  
  # plot heatmap
  cell_prop <- table(metadata$cluster10, metadata$cell_type) %>% scale()
  p <- Heatmap(cell_prop,show_row_dend = FALSE,show_column_dend = FALSE,rect_gp=gpar(col='gray'),column_names_rot = -45,name = "Proportion")
  pdf(file=paste0(patient,'_cell_prop_heatmap.pdf'),width = 7,height = 7)
  draw(p)
  dev.off()
  
  # plot tissue
  patient_folder <- patient_info[[patient]]$patient_folder
  xenium.obj <- LoadXenium(file.path(raw_dir,patient_folder), fov = "fov")
  xenium.obj <- AddMetaData(xenium.obj,metadata=metadata)
  xenium.obj <- xenium.obj[,!is.na(xenium.obj$cell_type)]
  xenium.obj$cluster10 <- as.factor(xenium.obj$cluster10)
  
  # community
  fov_community_p <- plot_cell_boundary(xenium.obj,fov_name='fov',group_by='cluster10',cols=community_colour,is_legend=T) +
    coord_flip() + scale_x_reverse()+
    labs(fill = "Community")
  ggsave(filename = paste0(patient,'_community.pdf'),plot=fov_community_p,height = patient_info[[patient]]$height,width =patient_info[[patient]]$width ,units = 'px',limitsize = FALSE)
}

