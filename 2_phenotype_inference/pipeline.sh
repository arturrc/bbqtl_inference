# BRIEF NOTES ON DEPENDENCIES ======================================
# > Python: version 3.4.5

# COMPILATION ======================================================
g++ error_correct_from_dict.cpp -o bin/error_correct_from_dict -O3 -std=c++0x -I$HOME/local/include -lpthread -fopenmp -D_GLIBCXX_PARALLEL
g++ PCG_BFA_with_prior_5.cpp -o bin/PCG_BFA_with_prior_5 -std=c++11 -O3 -I$HOME/local/include -fopenmp

# PIPELINE =========================================================
FILEHANDLE=('4NQO')

# > Parse reads <--------------------------------------------------
# Clean reads from both lanes and join them
# INPUT: 
# OUPUT: 
bash clean_reads.sh data/fastq/SEPSIS/p${seq_pool}_${timepoint}_S*_L001_R1_001.fastq.gz data/fastq/SEPSIS/p${seq_pool}_${timepoint}_S*_L001_R2_001.fastq.gz > ${scratch}/cleaned_p${seq_pool}_tp${timepoint}.txt
bash clean_reads.sh data/fastq/SEPSIS/p${seq_pool}_${timepoint}_S*_L002_R1_001.fastq.gz data/fastq/SEPSIS/p${seq_pool}_${timepoint}_S*_L002_R2_001.fastq.gz >> ${scratch}/cleaned_p${seq_pool}_tp${timepoint}.txt
cat ${scratch}/cleaned_p${seq_pool}_tp${timepoint}.txt | python read_parser.py > ${scratch}/parsed_p${seq_pool}_tp${timepoint}.txt

# > Error-correct BC sequences <------------------------------------
# INPUT: 
# OUPUT: 
. ./error_correct_from_dict -bc 5 -err 3 -dict data/BY_dictionary.txt ${scratch}/parsed_p${seq_pool}_tp${timepoint}.txt > ${scratch}/corrected_p${seq_pool}_tp${timepoint}.tmp
. ./error_correct_from_dict -bc 6 -err 3 -dict data/RM_dictionary.txt ${scratch}/corrected_p${seq_pool}_tp${timepoint}.tmp > ${scratch}/corrected_p${seq_pool}_tp${timepoint}.txt
rm ${scratch}/corrected_p${seq_pool}_tp${timepoint}.tmp

# > Remove duplicate reads <---------------------------------------
# INPUT: 
# OUPUT: 
grep -P ${p5_id}'\t'${p7_id}'\t'${col_id}'\t'${row_id}'\t' ${scratch}/corrected_p${seq_pool}_tp${timepoint}.txt | sort | uniq  > ${scratch}/indexed_p${seq_pool}_tp${timepoint}.txt

# > Generate read counts <-----------------------------------------
# INPUT: 
# OUPUT: 
cut -f 5,6 ${scratch}/indexed_p${seq_pool}_tp${timepoint}.txt | sort | uniq -c | sort -nr | awk "{print \$2,\$3,\$1}" OFS="\\t" > ${scratch}/counts_p${seq_pool}_tp${timepoint}.txt

# > Join read counts <--------------------------------------------
# INPUT: 
# OUPUT: 
python join_counts.py ${seq_pool} # generates ${scratch}/joined_counts_p${seq_pool}.txt

# > Infer fitness <------------------------------------------------
# INPUT: 
# OUPUT: 
prior_mean=0
prior_var=0
for i in {1..3}
do
	echo it $i mean = ${prior_mean} var = ${prior_var}
	. ./PCG_BFA_with_prior_5 ${prior_mean} ${prior_var} ${outdir}/filtered_counts_${env_label}.txt > ${scratch}/fitness_${env_label}_it${i}.txt
	prior_mean=$(Rscript calculate_fitness_mean.R ${scratch}/fitness_${env_label}_it${i}.txt)
	prior_var=$(Rscript calculate_fitness_var.R ${scratch}/fitness_${env_label}_it${i}.txt)
done

# END ==========================================================================
