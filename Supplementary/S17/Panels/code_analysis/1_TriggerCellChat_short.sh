#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 1st Dec 2021
# Title: TriggerCellChat.sh
# Goal: To run CellChat for different samples
# Usage: nohup TriggerCellChat.sh > nohup.out 2>&1&

# Usage: Rscript RunCellChat_short.R {sampleID} {counts} {metadata} {species} {outdir}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Run code
# ------------------------------------------------------------------

/Volumes/flan/Archive/SCC/1_sandbox/20240322_StartSingleCellAgain/13_CellChat/Level2/2_RunCellChat/

workdir=/Volumes/flan/Archive/SCC/1_sandbox/20240322_StartSingleCellAgain/13_CellChat/

# healthy
mkdir -p P1_healthy
Rscript RunCellChat_short.R P1_healthy "$workdir"/1_PrepareGenes/outdir/P1_healthy_counts.txt "$workdir"/1_PrepareGenes/outdir/P1_healthy_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P1_healthy

mkdir -p P2_healthy
Rscript RunCellChat_short.R P2_healthy "$workdir"/1_PrepareGenes/outdir/P2_healthy_counts.txt "$workdir"/1_PrepareGenes/outdir/P2_healthy_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P2_healthy

mkdir -p P3_healthy
Rscript RunCellChat_short.R P3_healthy "$workdir"/1_PrepareGenes/outdir/P3_healthy_counts.txt "$workdir"/1_PrepareGenes/outdir/P3_healthy_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P3_healthy

mkdir -p P4_healthy
Rscript RunCellChat_short.R P4_healthy "$workdir"/1_PrepareGenes/outdir/P4_healthy_counts.txt "$workdir"/1_PrepareGenes/outdir/P4_healthy_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P4_healthy

mkdir -p P5_healthy
Rscript RunCellChat_short.R P5_healthy "$workdir"/1_PrepareGenes/outdir/P5_healthy_counts.txt "$workdir"/1_PrepareGenes/outdir/P5_healthy_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P5_healthy

# cancer
mkdir -p P1_cancer
Rscript RunCellChat_short.R P1_cancer "$workdir"/1_PrepareGenes/outdir/P1_cancer_counts.txt "$workdir"/1_PrepareGenes/outdir/P1_cancer_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P1_cancer

mkdir -p P2_cancer
Rscript RunCellChat_short.R P2_cancer "$workdir"/1_PrepareGenes/outdir/P2_cancer_counts.txt "$workdir"/1_PrepareGenes/outdir/P2_cancer_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P2_cancer

mkdir -p P3_cancer
Rscript RunCellChat_short.R P3_cancer "$workdir"/1_PrepareGenes/outdir/P3_cancer_counts.txt "$workdir"/1_PrepareGenes/outdir/P3_cancer_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P3_cancer

mkdir -p P4_cancer
Rscript RunCellChat_short.R P4_cancer "$workdir"/1_PrepareGenes/outdir/P4_cancer_counts.txt "$workdir"/1_PrepareGenes/outdir/P4_cancer_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P4_cancer

mkdir -p P5_cancer
Rscript RunCellChat_short.R P5_cancer "$workdir"/1_PrepareGenes/outdir/P5_cancer_counts.txt "$workdir"/1_PrepareGenes/outdir/P5_cancer_meta_L2.txt human "$workdir"/Level2/2_RunCellChat/P5_cancer

