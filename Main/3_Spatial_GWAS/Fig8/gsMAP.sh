#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=250GB #250
#SBATCH -o out_%x_%j.txt
#SBATCH -e error_%x_%j.txt
#SBATCH --job-name=gsmap
#SBATCH --time=15:00:00
#SBATCH --partition=general
#SBATCH --account=a_nguyen_quan

module load miniconda3/4.12.0
#source activate /home/s4716765/.conda/envs/scTAR_cellranger
source activate /home/s4716765/.conda/envs/edger
cd /scratch/project/stseq/Prakrithi/skin_atlas/GWAS/gsmap_test
export LD_LIBRARY_PATH=/home/s4716765/.conda/envs/edger/lib:$LD_LIBRARY_PATH
export NUMBA_CACHE_DIR=/tmp/numba_cache

# Format summary stats if Z score is not given
gsmap format_sumstats  --sumstats 'melanoma_formatted' --out 'mel'


## Quick mode
gsmap quick_mode --workdir '/scratch/project/stseq/Prakrithi/skin_atlas/GWAS/gsmap_test' \
    --sample_name 'mel48974_all_raw_new_manh' --gsMap_resource_dir '/scratch/project/stseq/Prakrithi/skin_atlas/GWAS/gsMap_resource' \
    --hdf5_path 'mel48974_raw.h5ad' --annotation 'cell_type' --data_layer 'count' \
    --sumstats_file 'mel.sumstats.gz' --trait_name 'mel'
