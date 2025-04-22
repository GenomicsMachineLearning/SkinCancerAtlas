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

# load data ----------------
all_metadata <- read.csv(file.path(input_dir, 'all_metadata.csv'), check.names = FALSE)

for( i in c('cosmx','visium')){
    # i = 'cosmx'
    print(i)
    metadata <- all_metadata[all_metadata$dataset==i,]
    prop_logit <- getTransformedProps(clusters = metadata$cluster10, sample=metadata$batch, transform = "logit") 
    sample_info <- unique(metadata[,c('batch','tissue_type')])
    prop_dat <- t(prop_logit$Proportions) %>% as.data.frame() %>% left_join(sample_info,by=c('sample'='batch'))
    colnames(prop_dat)[2] <- 'Community'
    prop_dat$tissue_type <- factor(prop_dat$tissue_type, levels = c("BCC", "SCC", "Melanoma"))
    p <- ggboxplot(prop_dat, x = "tissue_type", y = "Freq",color = "tissue_type",
                add = "jitter", facet.by = "Community",
                short.panel.labs = FALSE,palette = "lancet", ncol = 5)
    # Use only p.format as label. Remove method name.
    p_point_cluster <- p + stat_compare_means(label = "p.format", vjust = 1.5,hjust= -0.5,size=5) +
                        labs(x=NULL,color=tools::toTitleCase(i))+
                        theme(
                            text = element_text(size = 30),          # Global text size
                            axis.text = element_text(size = 20),    # Axis text size
                            axis.title = element_text(size = 30),   # Axis title size
                            strip.text = element_text(size = 20),   # Facet label size
                            legend.text = element_text(size = 20),  # Legend text size
                            legend.title = element_text(size = 30),  # Legend title size
                            axis.text.x = element_text(angle = 45, hjust = 1)
                        )
    ggsave(filename = paste0(i,'_tissue_community_proportion.pdf'),plot=p_point_cluster,height=1800,width=5000,dpi=300,units='px')
}


# # barplot where patients as x and cluster proportion as y ----------
# prop_dat <- all_metadata %>% group_by(dataset,patient_cancer, cluster10) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   group_by(dataset,patient_cancer) %>%
#   mutate(proportion = count / sum(count)) %>%
#   ungroup()
# # Plot the barplot where cancer type as x and cluster proportion as y
# for(i in unique(prop_dat$dataset)){
#     # i = 'visium'
#     print(i)
#     sub_metadata <- prop_dat %>% filter(dataset==i)
#     sub_metadata$tissue_type <- str_split_fixed(sub_metadata$patient_cancer,"_",2)[,2] 
#     if (i == 'visium') sub_metadata$tissue_type <- ifelse(sub_metadata$tissue_type == "Mel", "Melanoma", sub_metadata$tissue_type)
#     if (i %in% c('cosmx','visium')) sub_metadata$tissue_type <- factor(sub_metadata$tissue_type, levels = c('BCC','SCC','Melanoma'))
#     sub_metadata <- sub_metadata %>% group_by(tissue_type,cluster10) %>% 
#         summarise(count = sum(count)) %>%
#         mutate(proportion = count / sum(count)) %>%
#         ungroup()
#     p <- sub_metadata %>% ggplot(aes(x = tissue_type, y = proportion, fill = as.factor(cluster10))) +
#             geom_bar(stat = "identity", position = "stack") +
#             scale_y_continuous(expand = c(0, 0.01)) +
#             theme_void() + #theme_minimal(base_size=16) +
#             theme(axis.text.x = element_text(angle = -45,hjust = 0,size = 30),panel.grid=element_blank(),
#                     plot.background = ggplot2::element_rect(fill = "white"),legend.margin = margin(l= -30),
#                     plot.margin = margin(r = 20),
#                     legend.text = element_text(size = 30),                          # Increase legend text size
#                     legend.title = element_text(size = 35),
#                     panel.border = element_blank())+ 
#             xlab(label = NULL)+labs(fill=paste0(tools::toTitleCase(i),"\nCommunity"))+
#             scale_fill_manual(values = community_colour)
#     ggsave(filename = paste0(i,'_cluster_tissue_type_proportion.png'),plot=p,dpi=300,width=3500,height=3000,units='px')
# }
