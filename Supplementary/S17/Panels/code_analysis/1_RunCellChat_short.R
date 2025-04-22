# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 1st December 2021
# Title: RunCellChat.R
# Goal: To Run Cellchat
# Usage: Rscript RunCellChat.R {sampleID} {counts} {metadata} {species} {outdir}
# species must be "human" or "mouse"
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

print("step 1: setting up")
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example:  Rscript RunCellChat.R {sampleID} {counts} {metadata} {species} {outdir}\n")
  quit()
}

sampleID <- args[1]
counts <- read.delim(args[2], row.names = 1)
metadata <- read.delim(args[3], row.names = 1)
species <- args[4]
outdir <- args[5]

# check that species is right, and set database
if (species == "human" | species == "mouse") {
  print("species good to go!")
} else {
  cat("ERROR: species must be `human` or `mouse`\n")
  quit()
} 

# load packages
library(CellChat) # must already be installed
bioc_packages <- c("")
r_packages <- c("patchwork", "NMF", "circlize", "ComplexHeatmap", "ggalluvial", "Matrix")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
# ## function to load bioconductor packages
# biocondpkgTest <- function(x) {
#   if (!require(x,character.only = TRUE)) {
#     source("http://www.bioconductor.org/biocLite.R")
#     BiocManager::install(x)
#     if(!require(x,character.only = TRUE)) stop (paste0(x, "bioconductor package not found"))
#   }
# }
# ## load packages
# for (b_pkg in bioc_packages) {
#   biocondpkgTest(b_pkg)
# }
for (r_pkg in r_packages) {
  baseRpkgTest(r_pkg)
}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

if (species == "human") {
  CellChatDB <- CellChatDB.human
  PPI <- PPI.human
} else {
  CellChatDB <- CellChatDB.mouse
  PPI <- PPI.mouse
} 
print(CellChatDB)

# ------------------------------------------------------------------
# MODIFY METADATA
# ------------------------------------------------------------------
print("step 2: modifying metadata")
colnames(metadata) <- c("metadata")
metadata <- metadata[rownames(metadata) %in% colnames(counts), , drop = FALSE]
metadata$metadata <- gsub("\\+", "", metadata$metadata)

# ------------------------------------------------------------------
# PREPARE CELLCHAT DATA AND DATABASE
# ------------------------------------------------------------------
print("step 3a: prepare cellchat data")
cellchat <- createCellChat(object = as.matrix(counts), meta = metadata, group.by = "metadata", do.sparse = FALSE)

print("step 3b: prepare cellchat database")
#CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
ggsave(filename = paste0(outdir, "databaseComposition.jpeg"))
# Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
write.table(CellChatDB$interaction, file = paste0(outdir, "databaseComposition.txt"), sep = "\t", quote = FALSE, col.names = NA)

# ------------------------------------------------------------------
# DATA PRE-PROCESSING
# ------------------------------------------------------------------

print("step 4: data pre-processing")
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI)

# ------------------------------------------------------------------
# CCI inference
# ------------------------------------------------------------------

print("step 5: CCI inference")
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
df.net_AllInferredLRs <- subsetCommunication(cellchat)
write.table(df.net_AllInferredLRs, file = paste0(outdir, "AllPredictedInteractions_LR.txt"), sep = "\t", quote = FALSE, col.names = NA)
df.net_SigPWs <- subsetCommunication(cellchat, slot.name = "netP")
write.table(df.net_SigPWs, file = paste0(outdir, "AllPredictedInteractions_SigPW.txt"), sep = "\t", quote = FALSE, col.names = NA)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------
print("step 6: done")
sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()