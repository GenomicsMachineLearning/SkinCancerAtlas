# library --------------
rm(list=ls());options(stringsAsFactors=FALSE)
library(tidyverse)
library(Seurat)
library(CellChat)

# set path ------------------
input_dir <- '/scratch/project/stseq/Feng/SkinCancerAtlas/clean-data'
output_dir <- '/home/uqfzha11/Feng/projects/SkinCancerAtlas/Revision/manuscript_code/fig6'
setwd(output_dir)
set.seed(111)

# cell chat analysis ----------------------
load(file = file.path(input_dir,"cosmx_c6_cellchat.RData"))

# 1. visualization the number of intersection and intersection weight ----------------------
# circle plot for each cancer type
pdf(file="cosmx_c6_interaction_number.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
for( i in c("BCC", "Melanoma")){
    # i <- "BCC"
    groupSize <- as.numeric(table(cellchat_res_list[[i]]@idents))
    netVisual_circle(cellchat_res_list[[i]]@net$count, vertex.weight = groupSize,
        weight.scale = T, label.edge= F, title.name = paste0("Number of interactions in ",i),
        vertex.label.cex = 1, edge.label.cex = 2)
}
dev.off()

pdf(file="cosmx_c6_interaction_strength.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
for( i in c("BCC", "Melanoma")){
    # i <- "BCC"
    groupSize <- as.numeric(table(cellchat_res_list[[i]]@idents))
    netVisual_circle(cellchat_res_list[[i]]@net$weight, vertex.weight = groupSize,
        weight.scale = T, label.edge= F, title.name = paste0("Interaction weights/strength in ",i))
}
dev.off()

# visualization the major sources and targets
group.cellType <- c('Melanocytes',"Fibroblast",'LC','NK','Endothelial Cell') %>% as.factor()
object.list <- lapply(cellchat_res_lift_list[c("BCC","Melanoma")], function(x) {mergeInteractions(x, group.cellType)})
cellchat_sub <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))

pdf(file='cosmx_c6_interaction_sig_number.pdf',width=10,height=10)
par(xpd=TRUE, mar=c(0,5,0,5))
p1 = netVisual_diffInteraction(cellchat_sub, weight.scale = T, measure = "count.merged", label.edge = T,
        vertex.label.cex = 3, edge.label.cex = 2,title.name=NA) #vertex.label.color='white',
dev.off()


pdf(file='cosmx_c6_interaction_sig_strength.pdf',width=10,height=10)
par( xpd=TRUE, mar=c(0,5,0,5))
p2 = netVisual_diffInteraction(cellchat_sub, weight.scale = F, measure = "weight.merged",label.edge = T,
        vertex.label.cex = 3, edge.label.cex = 2,title.name=NA)
dev.off()


# 3 Comparing pathway between cancer type ------------------
# based on the communication probability
lr_bubble_p <- netVisual_bubble(cellchat_res, sources.use = c('Melanocytes',"Fibroblast",'LC','NK','Endothelial Cell'), 
                targets.use = c('Melanocytes',"Fibroblast",'LC','NK','Endothelial Cell'), 
                comparison = c(2,1),  angle.x = 45, remove.isolate = T,thresh = 0.01,font.size=15,
                sort.by.source = F,sort.by.target = F, sort.by.source.priority = T,return.data = T)
df <- lr_bubble_p$communication
df$source.target
x_levels <- levels(df$source.target)
new_levels <- x_levels[order(str_split_fixed(x_levels,'\\(',n = 2)[,2])]
df$source.target <- factor(df$source.target, levels = new_levels)
df$source.target
y_levels <- levels(df$interaction_name_2)[c(5,7,10,12,14,16,6,8,9,11,13,15,17,1:4,18:28)]
df$interaction_name_2 <- factor(df$interaction_name_2, levels = y_levels)
netVisual_bubble2 <- function (object, df, sources.use = NULL, targets.use = NULL, signaling = NULL, 
          pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
          sort.by.source.priority = TRUE, color.heatmap = c("Spectral", 
                                                            "viridis"), n.colors = 10, direction = -1, thresh = 0.05, 
          comparison = NULL, group = NULL, remove.isolate = FALSE, 
          max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
          max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
          color.text = NULL, dot.size.min = NULL, dot.size.max = NULL, 
          title.name = NULL, font.size = 10, font.size.title = 10, 
          show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
          angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE){
  color.heatmap <- match.arg(color.heatmap)
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
                      color = prob, size = pval)) + geom_point(pch = 16) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = angle.x, hjust = hjust.x, 
                                     vjust = vjust.x), axis.title.x = element_blank(), 
          axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  if (is.null(dot.size.max)) {
    dot.size.max = max(df$pval)
  }
  if (is.null(dot.size.min)) {
    dot.size.min = min(df$pval)
  }
  g <- g + scale_radius(range = c(dot.size.min, dot.size.max), 
                        breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                                 sort(unique(df$pval))], name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white", limits = c(quantile(df$prob, 
                                                                            0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                                                                                         1, na.rm = T)), labels = c("min", "max")) + 
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  }
  else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                         title = "Commun. Prob."))
  }
  g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5 + length(dataset.name[comparison]), 
                       length(group.names0) * length(dataset.name[comparison]), 
                       by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                          color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      }
      else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, 
                                               "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                             2, stringr::str_length(dataset.name.order) - 
                                               1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "none")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  }
  else {
    return(g)
  }
} 
lr_bubble_p2 <- netVisual_bubble2(cellchat_res, df=df)

lr_bubble_p2 <- lr_bubble_p2 +
              theme(
                    legend.text = element_text(size = 15),    # Adjust legend text size
                    legend.title = element_text(size = 20),   # Adjust legend title size              
                    plot.margin = margin(t = 10, r = 10, b = 10, l = 30))  
ggsave(filename = 'cosmx_c6_LR_prob.pdf',plot=lr_bubble_p2,width=10,height=7)


# 4 visualize the certain pathway ------------------
intersect(cellchat_res@netP$Melanoma$pathways,cellchat_res@netP$BCC$pathways)
pathways.show <- c("COLLAGEN") 
# weight.max <- getMaxWeight(cellchat_res_lift_list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste0('cosmx_c6_chord_pathway_BCC.pdf'),width=7,height=7)
netVisual_aggregate(cellchat_res_lift_list[['BCC']], signaling = pathways.show, layout = "chord") 
                      #signaling.name = paste(pathways.show, names(cellchat_res_lift_list)[i]))
dev.off()
pdf(paste0('cosmx_c6_chord_pathway_Mel.pdf'),width=7,height=7)
netVisual_aggregate(cellchat_res_lift_list[['Melanoma']], signaling = pathways.show, layout = "chord", 
                     )
                      # paste(pathways.show, names(cellchat_res_lift_list)[i])
dev.off()
