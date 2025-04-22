# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
input_dir <- '/scratch/project/stseq/Feng/SkinCancerAtlas/clean-data'
output_dir <- '/home/uqfzha11/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/fig6'
setwd(output_dir)

# library and set default parameter  ----------------
library(tidyverse)
library(Seurat)
library(hoodscanR)
library(SpatialExperiment)
library(patchwork)
library(speckle)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
set.seed(111)

cell_type_col <- c('DC','Endothelial Cell','Fibroblast','KC Basal','KC Cornified','KC Granular','KC Differentiating','KC Hair','LC',
            'Macrophage','Melanocytes','NK','T Cell','Mast Cell','B Cell','Pericytes','Sweat gland related','Plasma','KC IFN','Monocytes','KC Dysplastic')
n_neightbor <- 10
unique_cluster <- 0:9
unique_colour <- brewer.pal(length(unique_cluster),'Paired')
names(unique_colour) <- unique_cluster
community_colour <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#6A3D9A", "#FDBF6F", "#FF7F00", "#CAB2D6", "#E31A1C")
names(community_colour) <- as.character(0:9)

# load data ----------------
all_metadata <- read.csv(file.path(input_dir,'all_metadata.csv'), check.names = FALSE)
cluster_prop_all <- NULL
for(i in c('cosmx','visium','xenium')){
    # i <- 'cosmx'
    metadata <- all_metadata[all_metadata$dataset==i,]
    cluster_prop <- metadata %>% group_by(cluster10) %>% 
        summarise(across(all_of(cell_type_col),mean))
    cluster_prop$platform_type <- i
    cluster_prop$cluster10 <- paste0(cluster_prop$platform_type,'_',cluster_prop$cluster10)
    cluster_prop_all <- rbind(cluster_prop_all,cluster_prop)
}

# plot by ComplexHeatmap combining the correlation and cell type combination ----------
celltype_m = data.frame(cluster_prop_all[,cell_type_col],row.names = cluster_prop_all$cluster10,check.names = F) %>% as.matrix()
cor_m = cor(t(celltype_m))

hc <- hclust(dist(cor_m))
hc_num <- data.frame(cutree(hc,k=10),row.names = row.names(celltype_m))
colnames(hc_num) <- "cluster"
Heatmap(celltype_m,row_split=hc_num$cluster, right_annotation = rowAnnotation(Community=as.character(hc_num$cluster)))
hc_num$anno <- case_match(hc_num$cluster, c(2) ~ 'Immune',c(7)~'Tumor',c(5,10)~'Stromal',c(1,3,4,6,8,9)~'KC')
hc_num['visium_2','anno'] <- 'Tumor'
hc_num$anno_order <- factor(hc_num$anno,levels = c('KC','Tumor','Stromal','Immune'))
ha_cor <- rowAnnotation(Community=hc_num$anno,col=list(Community=c('Immune'='red','Tumor'='blue','Stromal'='green','KC'='orange'))
            # annotation_legend_param = list(#Community = list(title_gp = gpar(fontsize = 20)),
            #                     title_gp = gpar(fontsize = 20))  
  )
ht_cor = Heatmap(cor_m, name = "Correlation", right_annotation=ha_cor,
                row_split=hc_num[,"anno_order"],cluster_row_slices = FALSE,
                row_title=NULL,column_names_rot=45,
                show_row_dend = FALSE, show_column_dend = FALSE
                #column_names_gp = gpar(fontsize = 20)
                )

cell_df <- data.frame(cell_type=colnames(celltype_m),row.names = colnames(celltype_m),check.names = F)
cell_df$anno <- case_match(cell_df$cell_type, c('B Cell','DC','Macrophage','T Cell','Mast Cell','NK','LC','Plasma','Monocytes') ~ 'Immune',
                                            c('Melanocytes')~'Tumor',
                                            c('Fibroblast','Endothelial Cell','Pericytes')~'Stromal',
                                            c('KC Basal','KC Differentiating','KC Hair','KC Cornified',"KC Granular",'Sweat gland related','KC IFN','KC Dysplastic')~'KC')
cell_df$anno <- factor(cell_df$anno,levels = c('KC','Tumor','Stromal','Immune'))
ht_prop = Heatmap(celltype_m, name = "Proportion",column_split = cell_df$anno,column_names_rot=-45,
                  cluster_column_slices = FALSE,column_gap = unit(0, "mm"),show_column_dend=F,column_title=NULL,
                 # row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20),
                  )

pdf('combined_community_celltype.pdf',height =6,width = 12)
par(family='serif')
draw(ht_cor+ht_prop,padding = unit(c(3, 15, 5, 3), "mm"))
dev.off()
