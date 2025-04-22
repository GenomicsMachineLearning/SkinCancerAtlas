# library ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
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
input_dir <- '/scratch/project/stseq/Feng/SkinCancerAtlas/clean-data'
output_dir <- '/home/uqfzha11/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/fig6'
setwd(output_dir)
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

community_colour <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#6A3D9A", "#FDBF6F", "#FF7F00", "#CAB2D6", "#E31A1C")
names(community_colour) <- as.character(0:9)

#samples <- readRDS(file=file.path(input_dir,'cosmx_cutoff50_label_list.rds'))
all_metadata <- read.csv(file.path(input_dir,'all_metadata.csv'))

# cosmx ----------------
cosmx_seu <- readRDS(file.path(input_dir, 'cosmx_spatial_all.rds'))
community <- all_metadata[all_metadata$dataset!='xenium',]
row.names(community) <- community$cell_id
cosmx_seu <- AddMetaData(cosmx_seu,metadata=community[,c('cell_id','cluster10')])
cosmx_seu <- cosmx_seu[,!is.na(cosmx_seu$cluster10)]
samples <- SplitObject(cosmx_seu,split.by='orig.ident')
raw_dir <- '/home/uqfzha11/Feng/SkinCancerAtlas/process-data/flat_file_format_raw_data/'
dir(raw_dir)

## B18_BCC_13 ----------------------
patient = "B18_BCC"
print(patient)
fov_id=13
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples[[patient]]@meta.data
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

fov_community_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='cluster10',cols=community_colour,is_legend=T) +
                        coord_flip() + theme(legend.position = "none")
ggsave(filename = paste0(patient,'_',fov_id,'_community.pdf'),plot=fov_community_p,width = 5472,height = 3648,units='px',dpi=200)

## B18_SCC_14 ----------------------
patient = "B18_SCC"
print(patient)
fov_id=14
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples[[patient]]@meta.data
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

fov_community_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='cluster10',cols=community_colour,is_legend=T) +
                        coord_flip() + theme(legend.position = "none")
ggsave(filename = paste0(patient,'_',fov_id,'_community.pdf'),plot=fov_community_p,width = 5472,height = 3648,units='px',dpi=200)


## 48974-2B_Melanoma_12 ----------------------
dir(raw_dir)
patient = "48974-2B"
print(patient)
fov_id=12
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples[[patient]]@meta.data
row.names(metadata) <- str_remove(metadata$cell_id,paste0(patient,'_'))
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)
nano.obj$cluster10 <- as.factor(nano.obj$cluster10)

# tumor community
community_tumor_colour <- c(rep('gray',6),'red',rep('gray',3))
names(community_tumor_colour) <- as.character(0:9)
fov_community_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='cluster10',cols=community_tumor_colour,is_legend=T) +
                        coord_flip() #+ theme(legened.position = "none")
ggsave(filename = paste0(patient,'_',fov_id,'_tumor_community.pdf'),plot=fov_community_p,width = 10,height =7)

# community
fov_community_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='cluster10',cols=community_colour,is_legend=T) +
                        coord_flip() + theme(legend.position = "none")
ggsave(filename = paste0(patient,'_',fov_id,'_community.pdf'),plot=fov_community_p,width = 5472,height = 3648,units='px',dpi=200)


# xenium ----------------
raw_dir <- "/home/uqfzha11/Feng/SkinCancerAtlas/process-data/Xenium_data"
patient_folder <- 'output-XETG00114__0017002__6475-07FC__20231129__074525'
patient <- '6475-07FC'
xenium.obj <- LoadXenium(file.path(raw_dir,patient_folder), fov = "fov")
metadata <- all_metadata[all_metadata$patient_cancer==paste0(patient,'_Mel'),]
row.names(metadata) <- metadata$cell_id
xenium.obj <- AddMetaData(xenium.obj,metadata=metadata)
xenium.obj <- xenium.obj[,!is.na(xenium.obj$cell_type)]
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj$cluster10 <- as.factor(xenium.obj$cluster10)

# community
fov_community_p <- plot_cell_boundary(xenium.obj,fov_name='fov',group_by='cluster10',cols=community_colour,is_legend=T) +
                        coord_flip() + scale_x_reverse()+ 
                        theme(legend.position = "none")
ggsave(filename = paste0(patient,'_community.pdf'),plot=fov_community_p,width =14993 ,height = 4299,units='px')

# tumor community
community_tumor_colour <- c(rep('gray',2),'red',rep('gray',4),'red',rep('gray',2))
names(community_tumor_colour) <- as.character(0:9)
fov_community_p <- plot_cell_boundary(xenium.obj,fov_name='fov',group_by='cluster10',is_legend=T,cols=community_tumor_colour,flip_xy = FALSE)+
                        coord_flip() + scale_x_reverse()
                        #labs(fill="Community") + theme( legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 60), legend.title = element_text(size = 70) )
ggsave(filename = paste0(patient,'_tumor_community.pdf'),plot=fov_community_p,width = 14993,height = 4299,units='px')

# plot legend
col_names <- names(community_colour)
order_cell <- order(col_names)
colors <- community_colour[order_cell]
col_names <- names(colors)
pdf(file.path(output_dir,'community_legend.pdf'), width = 5, height = 5)
par(mar = c(0, 0, 0, 0),mgp=c(0,0,0))
plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = NULL)
legend("center", legend = col_names, col = colors, cex=0.5,pch=15, pt.cex=1,pt.lwd = 1, bty = "n",border = "white",
         x.intersp = 1, y.intersp = 0.8,ncol=1
        )
dev.off()
