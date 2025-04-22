
library(ComplexHeatmap)
library(circlize)


data <- read.csv("~/Downloads/LR_for_plot",sep="\t",header = TRUE)

# Combine Ligand and Receptor into a single axis
data$Ligand_Receptor <- paste(data$Ligand, "→", data$Receptor)

# Create the heatmap data matrix
heatmap_data <- table(data$Ligand_Receptor, paste(data$Ligand_Cell_Type, "→", data$Receptor_Cell_Type))

# Ensure all samples are represented in annotations
unique_samples <- unique(data$Sample)
sample_colors <- setNames(c("red", "blue", "green", "purple","orange","pink","brown","yellow"), unique_samples)

# Map rows to samples
row_sample_map <- data.frame(
  Ligand_Receptor = unique(data$Ligand_Receptor),
  Sample = sapply(unique(data$Ligand_Receptor), function(x) {
    sample <- data[data$Ligand_Receptor == x, "Sample"]
    if (length(sample) > 0) sample[1] else NA
  })
)
row_colors <- sample_colors[row_sample_map$Sample]

# Create the heatmap
Heatmap(
  as.matrix(heatmap_data),
  name = "Count",
  row_title = "Ligand → Receptor",
  column_title = "Ligand Cell Type → Receptor Cell Type",
  col = colorRamp2(c(0, max(heatmap_data)), c("white", "darkblue")),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  row_names_side = "left",
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  top_annotation = HeatmapAnnotation(
    Sample = anno_simple(factor(row_sample_map$Sample), col = sample_colors)
  ),
  row_split = row_sample_map$Sample  # Split rows by sample
)

# Add a legend for the sample annotations
legend <- Legend(labels = names(sample_colors), legend_gp = gpar(fill = sample_colors), title = "Samples")
draw(legend, x = unit(1, "npc") - unit(1, "cm"), just = "right")







library(ggplot2)
library(dplyr)
library(ggtext)

# Read data
data <- read.csv("~/Downloads/LR_for_plot2", sep = "\t", header = TRUE)
data <- read.csv("~/Downloads/LR_for_plots_pos_only", sep = "\t", header = TRUE)
data <- read.csv("~/Downloads/LR_plot_new", sep = "\t", header = TRUE)
# Create combined columns
data <- data %>%
  mutate(
    Ligand_Receptor = paste(Ligand, "→", Receptor),
    Ligand_Receptor_CellType = paste(Ligand_Cell_Type, "→", Receptor_Cell_Type)
  )

# Create a table for heatmap values
heatmap_data <- as.data.frame(table(data$Sample, data$Ligand_Receptor, data$Ligand_Receptor_CellType))
colnames(heatmap_data) <- c("Sample", "Ligand_Receptor", "Ligand_Receptor_CellType", "Count")

# Adjust the Count to presence/absence
heatmap_data$Count <- ifelse(heatmap_data$Count > 0, 1, 0)

# Add `ref` and `ref_color` information
# for old
#heatmap_data <- heatmap_data %>%
#  left_join(data %>% select(Ligand_Receptor, ref) %>% distinct(), by = "Ligand_Receptor") %>%
#  mutate(ref_color = ifelse(ref == "mel", "#74c0fc", "#ffa94d")) # Map `ref` to colors
# for new with three ref
heatmap_data <- heatmap_data %>% 
  left_join(data %>% select(Ligand_Receptor, ref) %>% distinct(), by = "Ligand_Receptor") %>% 
  mutate(
    ref_color = case_when(
      ref == "mel" ~ "#74c0fc",   # Blue for "mel"
      ref == "scc" ~ "maroon",   # Orange for "scc"
      ref == "bcc" ~ "#37b24d",   # Green for "bcc"
      TRUE ~ NA_character_        # Handle any unexpected values
    )
  )


# Order `Ligand_Receptor` by `ref`
heatmap_data <- heatmap_data %>%
  arrange(ref) %>%
  mutate(Ligand_Receptor = factor(Ligand_Receptor, levels = unique(Ligand_Receptor)))

# Order samples
ordered_samples <- c("mel48974", "mel66487", "mel6767", "B18_SCC", "E15_SCC", "P30_SCC", "P13_SCC",
                     "B18_BCC","F21_BCC")
heatmap_data <- heatmap_data %>%
  mutate(Sample = factor(Sample, levels = ordered_samples))

# Define custom colors for samples
sample_colors <- c(
  "mel48974" = "#d0ebff", # light blue
  "mel66487" = "#74c0fc", # medium blue
  "mel6767"  = "#1c7ed6", # dark blue
  "B18_BCC"      = "#82c91e",
  "F21_BCC" = "lightgreen",  # green
  "B18_SCC"      = "pink", # light red-orange
  "E15_SCC"      = "maroon", # medium red-orange
  "P30_SCC"      = "#e03132", # dark red-orange
  "P13_SCC"      = "darkred", # darkest red-orange
  "grey"         = "#d3d3d3"  # light grey for absent data
)

# Create HTML-styled y-axis labels
heatmap_data <- heatmap_data %>%
  mutate(Ligand_Receptor_label = paste0(
    "<span style='color:", ref_color, "'>", Ligand_Receptor, "</span>"
  ))

# Ensure Y-axis (`Ligand_Receptor_label`) is ordered by ref: mel, scc, bcc
heatmap_data <- heatmap_data %>%
  mutate(
    Ligand_Receptor_label = factor(
      Ligand_Receptor_label,
      levels = heatmap_data %>%
        filter(!is.na(ref)) %>%           # Exclude rows without `ref`
        arrange(
          factor(ref, levels = c("mel", "scc", "bcc")), # Explicitly order `ref`
          Ligand_Receptor_label                           # Secondary ordering within `ref`
        ) %>%
        pull(Ligand_Receptor_label) %>%
        unique()
    )
  )

# Plot the heatmap
ggplot(heatmap_data %>% filter(Count == 1), aes( # Only include rows where Count == 1
  x = Ligand_Receptor_CellType, 
  y = Ligand_Receptor_label, 
  fill = Sample
)) +
  geom_tile(color = "black") + # Add gridlines for separation
  scale_fill_manual(
    values = sample_colors, 
    name = "Sample", 
    na.translate = TRUE,  # Ensure NA is mapped
    na.value = "#d3d3d3"  # Light grey for missing data
  ) +
  facet_grid(. ~ Sample, scales = "free_x", space = "free") + # Separate samples in one row
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = ggtext::element_markdown(size = 10), # Use `element_markdown` for styled text
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    title = "Ligand-Receptor Interactions Heatmap",
    fill = "Sample"
  )


##############################

### cell types HM

library(dplyr)
library(patchwork)
df<-read.csv("~/Downloads/gsmap_cell_types_plot",sep="\t",header=TRUE)
# Normalize P_Cauchy for better visualization (optional)
df <- df %>%
  mutate(P_Cauchy_scaled = -log10(P_Cauchy))  # Transform to -log10 scale

# Plot heatmap
ggplot(df, aes(x = Sample, y = Annotation, fill = P_Cauchy_scaled)) +
  geom_tile(color = "black") +  # Add gridlines
  scale_fill_gradient(
    low = "blue", high = "red", 
    name = "-log10(P Cauchy)", 
    na.value = "grey90"
  ) +
  facet_wrap(~ SNPs, scales = "free", nrow = 1) +  # Facet by SNPs
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Heatmap of Samples vs. Annotations",
    x = "Sample",
    y = "Annotation",
    fill = "-log10(P Cauchy)"
  )


# each group instead of facet
# Subset data where SNPs == "ALL"
df_filtered <- df %>% filter(SNPs == "ALL")

# Normalize P_Cauchy within each sample
df_filtered <- df_filtered %>%
  group_by(Sample) %>%
  mutate(Scaled_P_Cauchy = (P_Cauchy - min(P_Cauchy)) / (max(P_Cauchy) - min(P_Cauchy))) %>%
  ungroup()

# Order x-axis labels based on Cancer
df_filtered <- df_filtered %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(Cancer)])))

# Plot the heatmap
ggplot(df_filtered, aes(x = Sample, y = Annotation, fill = -Scaled_P_Cauchy*10)) +
  geom_tile(color = "white") + # White gridlines
  scale_fill_gradient(low = "blue", high = "red", name = "Scaled P Cauchy") + # Color scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Heatmap of Sample vs Annotation", fill = "Scaled P Cauchy")


# for best distinguishing
# Define custom color scale
custom_colors <- c("blue", "lightblue", "yellow","red")

# Logarithmic scale breaks
log_breaks <- c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5) # Adjusted for finer control

# Plot heatmap
ggplot(df_filtered, aes(x = Sample, y = Annotation, fill = Scaled_P_Cauchy)) +
  geom_tile(color = "white") + # White gridlines for clarity
  scale_fill_gradientn(
    colors = custom_colors,
    trans = "log10", # Apply logarithmic transformation
    breaks = log_breaks, # Set custom breaks
    labels = signif(log_breaks, digits = 2), # Format labels for clarity
    name = "Scaled P Cauchy"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Heatmap with Logarithmic Scale", fill = "Scaled P Cauchy")





library(ggplot2)
library(dplyr)

# Identify the annotation with the lowest P_Cauchy for each sample
df_highlighted <- df_filtered %>%
  group_by(Sample) %>%
  mutate(Is_Lowest = Scaled_P_Cauchy == min(Scaled_P_Cauchy)) %>% # Flag the lowest value
  ungroup()

# Define custom colors and highlight settings
custom_colors <- c("blue", "lightblue", "yellow", "red")
highlight_color <- "black" # Color for highlighting the lowest value
  
# Logarithmic scale breaks
log_breaks <- c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5) 

# Plot heatmap
ggplot(df_highlighted, aes(x = Sample, y = Annotation, fill = Scaled_P_Cauchy)) +
  geom_tile(aes(color = Is_Lowest), size = 0.8) + # Highlight lowest values
  scale_fill_gradientn(
    colors = custom_colors,
    trans = "log10",
    breaks = log_breaks,
    labels = signif(log_breaks, digits = 2),
    name = "Scaled P Cauchy"
  ) +
  scale_color_manual(
    values = c("TRUE" = highlight_color, "FALSE" = "white"), 
    guide = "none" # Hide legend for highlighting
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Heatmap Highlighting Lowest P_Cauchy",
    fill = "Scaled P Cauchy"
  )

# professional color 
library(ggplot2)
library(dplyr)
library(viridis) # For professional color scales

df_filtered <- df %>% filter(SNPs == "SUG")
df_filtered <- df_filtered %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(Cancer)])))
cancer_order <- c("Mel", "SCC", "BCC")

# Reorder the Cancer column and Sample column based on Cancer
df_filtered <- df_filtered %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),  # Set custom Cancer order
    Sample = reorder(Sample, as.numeric(Cancer))    # Reorder Sample by Cancer
  )
# Scale each sample so that the lowest P_Cauchy becomes 1
df_scaled <- df_filtered %>%
  group_by(Sample) %>%
  mutate(
    Scaled_P_Cauchy = P_Cauchy / min(P_Cauchy), # Scale relative to the minimum
    Highlight = P_Cauchy == min(P_Cauchy)      # Flag the lowest value
  ) %>%
  ungroup()

# Plot heatmap with professional color scale
ggplot(df_scaled, aes(x = Sample, y = Annotation, fill = Scaled_P_Cauchy)) +
  geom_tile(color = "black", size = 0.5) + # Add gridlines for clarity
  geom_text(aes(label = ifelse(Highlight, "*", "")), color = "white", size = 5) + # Mark the lowest values
  scale_fill_viridis(
    option = "mako",          # Choose professional color scale (e.g., "viridis", "cividis", "mako", or "plasma")
    trans = "log10",          # Use log transformation for better distinction
    breaks = c(1, 2, 5, 10),  # Adjust breaks for clarity
    labels = signif,
    name = "Scaled P Cauchy"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top"
  ) +
  labs(
    title = "Heatmap with Scaled P_Cauchy (Professional Colors)",
    fill = "Scaled P_Cauchy"
  )

