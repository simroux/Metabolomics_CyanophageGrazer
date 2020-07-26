### Metabolomics_CyanophageGrazer
## Scripts used in the analysis of metabolomics data from the Cyanophage - Cyanobacteria - Grazer experiment
### First to get initial matrices (default uses peak area, option "-ph" used to get matrices of peak height instead)
./get_full_matrices.pl -ph run
### Then moving on to sample-sample correlation to identify outliers, and compound-compound correlation to identify duplicates
# This is done in R, with the script check_correlation.R
### After manual curation, should have 2 files: Correlated_compounds_decision_clean.tsv & Correlated_samples.tsv
### Now generate some clean matrices (same as get_full_matrices.pl, option "-ph" used to get matrices of peak height instead of area)
./get_final_matrices.pl -ph run
### Then back to R, using the script process_clean_matrices.R
# In R, FC and FDR of individual compounds in individual experiments will be calculated
### These will be used to generate the input for pheatmap with the following scripts (first for Media, second for Pellet)
./get_heatmap_matrices.pl -c Media_ph_s
./get_heatmap_matrices.pl -c Pellet_cnorm_ph_s
### Heatmaps are then generated in R using the script Plot_heatmaps.R
