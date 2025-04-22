# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
raw_dir <- '/home/uqfzha11/Feng/SkinCancerAtlas/process-data/flat_file_format_raw_data/'
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


# cosmx ----------------
samples <- readRDS(file.path(clean_dir, "cosmx_spatial_all.rds"))

## B18_BCC fov 13 ------------------
# cell type
patient = "B18_BCC"
print(patient)
fov_id=13
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples@meta.data[samples$orig.ident==patient,]
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

fov_cell_type_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='predicted.id',molecules=c("KRT17",'S100A8','CD63'),
                            cols=unique_colour,is_legend=F,mols.size = 2,mols.cols=c('blue','green','red'))+
                        coord_flip()
ggsave(filename = paste0(patient,'_',fov_id,'_celltype.pdf'),plot=fov_cell_type_p,width = 5472,height = 3648,units='px',dpi=200)

# marker
# deg <- FindMarkers(samples[[patient]], group.by="predicted.id",ident.1 = "Melanocyte" )
# seu_obj <- cosmx_seu[, cosmx_seu$patient_cancer_fov==paste0(patient,"_",fov_id)]
# table(seu_obj$annotation,useNA = 'ifany')
# deg <- FindMarkers(seu_obj, group.by="annotation",ident.1 = "Tumor" )
# deg[deg$p_val_adj<0.01 & deg$avg_log2FC>0 & deg$pct.1>0.3,]
# marker <- 'KRT17' #'LGALS1'
# fov_marker_p <- plot_marker_boundary(nano.obj,fov_name='global',features=marker)+
#                         coord_flip()+
#                         theme( legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 60), legend.title = element_text(size = 70) )+
#                         scale_fill_gradient(low = "gray", high = "red",breaks=c(0,1,2,3))
# ggsave(filename = paste0(patient,'_',fov_id,'_',marker,'.pdf'),plot=fov_marker_p,width = 6000,height = 3648,units='px',dpi=200)

# # note the clouds of TPSAB1 molecules denoting mast cells
# ImageDimPlot(nano.obj, fov = "global", molecules = c("KRT17"), coord.fixed = FALSE) + 
#                 theme_void() + labs(title = NULL) 

## B18_SCC_14 ----------------------
# cell type
patient = "B18_SCC"
print(patient)
fov_id=14
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples@meta.data[samples$orig.ident==patient,]
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

# cell type
fov_cell_type_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='predicted.id',molecules=c("KRT17",'S100A8','CD63'),
                            cols=unique_colour,is_legend=F,mols.size = 2,mols.cols=c('blue','green','red'))+
                        coord_flip()
ggsave(filename = paste0(patient,'_',fov_id,'_celltype.pdf'),plot=fov_cell_type_p,width = 5472,height = 3648,units='px',dpi=200)

# marker
# deg <- FindMarkers(samples[[patient]], group.by="predicted.id",ident.1 = "Melanocyte" )
# seu_obj <- cosmx_seu[, cosmx_seu$patient_cancer_fov==paste0(patient,"_",fov_id)]
# table(seu_obj$annotation,useNA = 'ifany')
# deg <- FindMarkers(seu_obj, group.by="annotation",ident.1 = "Tumor" )
# deg[deg$p_val_adj<0.01 & deg$avg_log2FC>0 & deg$pct.1>0.3,]
# marker <- 'S100A8' #'IGFBP7'
# fov_marker_p <- plot_marker_boundary(nano.obj,fov_name='global',features=marker)+
#                 coord_flip()+ 
#                 theme( legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 60), legend.title = element_text(size = 70) )+
#                         scale_fill_gradient(low = "gray", high = "red")
# ggsave(filename = paste0(patient,'_',fov_id,'_',marker,'.pdf'),plot=fov_marker_p,width = 6000,height = 3648,units='px',dpi=200)

## 48974-2B_Melanoma_12 ----------------------
patient = "48974-2B"
print(patient)
fov_id=12
nano.obj <- LoadNanostring(file.path(raw_dir,patient),fov.filter=fov_id)
metadata <- samples@meta.data[samples$orig.ident==patient,]
row.names(metadata) <- paste0(metadata$cell_ID,'_',metadata$fov)
nano.obj <- add_meta_SCT(nano.obj,metadata = metadata)

# cell type
fov_cell_type_p <- plot_cell_boundary(nano.obj,fov_name='global',group_by='predicted.id',molecules=c("KRT17",'S100A8','CD63'),
                            cols=unique_colour,is_legend=F,mols.size = 2,mols.cols=c('blue','green','red'))+
                        coord_flip()
ggsave(filename = paste0(patient,'_',fov_id,'_celltype.pdf'),plot=fov_cell_type_p,width = 5472,height = 3648,units='px',dpi=200)

# marker
# deg <- FindMarkers(samples[[patient]], group.by="predicted.id",ident.1 = "Melanocyte" )
# seu_obj <- cosmx_seu[, cosmx_seu$patient_cancer_fov==paste0(patient,"_Melanoma_",fov_id)]
# table(seu_obj$annotation,useNA = 'ifany')
# deg <- FindMarkers(seu_obj, group.by="annotation",ident.1 = "Tumor" )
# deg[deg$p_val_adj<0.01 & deg$avg_log2FC>0 & deg$pct.1>0.3,]
# marker <- 'CD63' #'GPNMB'
# fov_marker_p <- plot_marker_boundary(nano.obj,fov_name='global',features=marker) +                        
#                         coord_flip() + 
#                         theme( legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 60), legend.title = element_text(size = 70) )+
#                         scale_fill_gradient(low = "gray", high = "red")
# ggsave(filename = paste0(patient,'_',fov_id,'_',marker,'.pdf'),plot=fov_marker_p,width = 6000,height = 3648,units='px',dpi=200)

# plot the legend -------------
# col_names <- c("Normal epidermis", "Tumor", "Immune", "Stromal")
# colors <- c("#00FF00", "#FF0000", "#0000FF", "#008000")
# pdf(file.path(output_dir,'IF_legend.pdf'))
# par(mar = c(0, 0, 0, 0),mgp=c(0,0,0))
# plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = NULL)
# legend("center", legend = col_names, col = colors, cex=0.8,pch=15, pt.cex=1.5,pt.lwd = 1, bty = "n",border = "white",
#          x.intersp = 1, y.intersp = 0.85,
#         )
# dev.off()

# col_names <- names(unique_colour)
# order_cell <- order(col_names)
# colors <- unique_colour[order_cell]
# col_names <- names(colors)
# pdf(file.path(output_dir,'celltype_legend.pdf'), width = 10, height = 10)
# par(mar = c(0, 0, 0, 0),mgp=c(0,0,0))
# plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = NULL)
# legend("center", legend = col_names, col = colors, cex=0.5,pch=15, pt.cex=1,pt.lwd = 1, bty = "n",border = "white",
#          x.intersp = 1, y.intersp = 0.8,ncol=1
#         )
# dev.off()

# col_names <- c('Adipose','Exocrine','Immune','Stroma','Tumor','Vascular','Other','Keratinocytes','Follicle')
# colors <- c('#00ffff','#ffff00','#0000ff','#b6e4b6','#000000','#ff0000','#ffc800','#ffccff','#ff00ff')
# pdf(file.path(output_dir,'xenium_annot_legend.pdf'))
# par(mar = c(0, 0, 0, 0),mgp=c(0,0,0))
# plot(NULL, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, main = NULL)
# legend("center", legend = col_names, col = colors, cex=0.5,pch=15, pt.cex=1,pt.lwd = 1, bty = "n",border = "white",
#          x.intersp = 1, y.intersp = 0.8
#         )
# dev.off()
