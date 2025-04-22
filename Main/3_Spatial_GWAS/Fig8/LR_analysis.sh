awk -F',' 'NR==FNR {pairs[$1" "$2]; next} {if (($1" "$2) in pairs) print $0}' mel_up_LRpairs FS=',' mel48974_all_raw/report/mel/48974_all_sig_raw_mel_Gene_Diagnostic_Info.csv > mel48974_all_sig_raw_LRpairs


awk '
BEGIN {
    FS = ","; OFS = "\t";  # Set delimiter for file2
}
NR == FNR {
    # Read file2: Gene -> Cell Type mapping
    gsub(/"/, "", $1);  # Remove any quotes
    gene_celltype[$1] = $2;  # Map gene -> cell_type
    next;
}
{
    # Process file1 (tab-separated)
    ligand = $1; receptor = $2;

    # Get cell types for ligand and receptor
    ligand_celltype = (ligand in gene_celltype) ? ligand "_" gene_celltype[ligand] : "NA";
    receptor_celltype = (receptor in gene_celltype) ? receptor "_" gene_celltype[receptor] : "NA";

    # Output matched pairs
    print ligand_celltype, receptor_celltype;
}'  B18_SCC_all_eff_n_raw/report/scc_eff_n/B18_SCC_all_eff_n_raw_scc_eff_n_Gene_Diagnostic_Info.csv   FS='\t' connection_LRdb_LR_pairs > B18_SCC_LRs


grep "CCL5" mel*LRs | awk '!/NA/ {print}' | sed 's/ //g'


grep -f scc_up_ligands B18_SCC_LRs | awk '!/NA/ {print}' | grep -f scc_up_receptors | sed 's/ //g' | less

grep -f mel_up_ligands mel66487_LRs | awk '!/NA/ {print}' | grep -f mel_up_receptors | sed 's/ //g' | less




### loop pos corr only
find . -path "*/report/*/*_Gene_Diagnostic_Info.csv" -type f | while read file; do
    awk -F, '$4 >= 0.3' "$file" > "${file%.csv}_pos.csv"
done



# Loop to find LR pairs across samples that have PCC>0.3
#!/bin/bash
# Loop over all files matching the pattern */report/*/*_pos.csv
for input2 in */report/*/*_pos.csv; do
    # Extract the sample name (basename without extension)
    sample_name=$(basename "$input2" .csv)
    
    # Debugging: print the file being processed
    echo "Processing file: $input2"
    
    # Run the awk command for each input2 and process the output
    awk '
    BEGIN {
        FS = ","; OFS = "\t";  # Set delimiter for file2
    }
    NR == FNR {
        # Read file2: Gene -> Cell Type mapping
        gsub(/"/, "", $1);  # Remove any quotes
        gene_celltype[$1] = $2;  # Map gene -> cell_type
        next;
    }
    {
        # Process file1 (tab-separated)
        ligand = $1; receptor = $2;

        # Get cell types for ligand and receptor
        ligand_celltype = (ligand in gene_celltype) ? ligand "_" gene_celltype[ligand] : "NA";
        receptor_celltype = (receptor in gene_celltype) ? receptor "_" gene_celltype[receptor] : "NA";

        # Output matched pairs
        print ligand_celltype, receptor_celltype;
    }
    ' "$input2"  FS='\t'  mel_up_LR_pairs | awk '!/NA/ {print}' | sed 's/ //g' > mel/"${sample_name}_LRs.txt"

    # Indicate processing is done for the current file
    echo "Processed $input2 -> ${sample_name}_LRs.txt"
done
