# BRIEF NOTES ON DEPENDENCIES ======================================
# > Python: version 3.4.5

# PIPELINE =========================================================
# This pipeline does not include sample data as it takes all compressed sensing
# files in the NCBI SRA repository (i.e. *_cs_*).

# > Get read counts <-----------------------------------------------
# For each compressed sensing file, do the same thing as done in 
# 2_phenotype_inference to get to count files.

# > First pass address assignment <---------------------------------
# Accepted barcodes in first_pass_spore_addresses.txt in format BC1, BC2, batch,
# 	set, plate, well.
# Barcodes passed to lap in bcs_for_lap.txt in format BC1, BC2, batch count array,
# 	set count array, plate count array, row count array, col count array.
python cs_address_first_pass.py

# > Create cost matrix <--------------------------------------------
# Takes files blanks_all.txt and first_pass_spore_addresses.txt to generate wells
# 	to match, takes bcs_for_lap.txt to get bcs to match.
# Stores cost matrix as lap_cost_matrix.txt, of format:
# 	line 1: barcode list
# 	line 2: well index list
# 	then cost matrix (bc=row,index=column).
python lap_cost_matrix.py

# > Assign remaining wells <---------------------------------------
# Get LAP assignments for remaining wells.
python lap_solve.py

# END ==========================================================================
