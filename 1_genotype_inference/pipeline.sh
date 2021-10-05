# BRIEF NOTES ON DEPENDENCIES ======================================
# > Trimmomatic: version 0.39
# > Python: version 3.4.5
# > bowtie2
# > samtools

# PIPELINE =========================================================
FILEHANDLE=('geno_B2T1_subset')

# > Parse reads <--------------------------------------------------
# INPUT: *.fastq.gz
# OUPUT: *.trim.fastq.gz
python parse_genotyping.py \
	--R1=sample_data/${FILEHANDLE}_R1.fastq.gz \
	--R2=sample_data/${FILEHANDLE}_R2.fastq.gz

# > Trim <----------------------------------------------------------
# INPUT: *.trim.fastq.gz
# OUPUT: *.R[12].full_trim.fastq.gz
java -jar $HOME/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	-threads $SLURM_CPUS_ON_NODE \
	-phred33 ${FILEHANDLE}_R1.trim.fastq.gz \
	${FILEHANDLE}_R2.trim.fastq.gz \
	${FILEHANDLE}.R1.full_trim.fastq.gz \
	${FILEHANDLE}.R1.unpaired.full_trim.fastq.gz \
	${FILEHANDLE}.R2.full_trim.fastq.gz \
	${FILEHANDLE}.R2.unpaired.full_trim.fastq.gz \
	ILLUMINACLIP:$HOME/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:TRUE \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# > Bowtie2 to BY <----------------------------------------------------
# INPUT: *.R[12].full_trim.fastq.gz
# OUTPUT: *.sam
bowtie2 -t -p 8 --local  -L 20 --ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals -i S,1,0.25 \
	--score-min L,1,0.75  --reorder --no-unal -x indexes/BY4742_fixed \
	-U ${FILEHANDLE}.R1.full_trim.fastq.gz,${FILEHANDLE}.R2.full_trim.fastq.gz,${FILEHANDLE}.R1.unpaired.full_trim.fastq.gz,${FILEHANDLE}.R2.unpaired.full_trim.fastq.gz \
	-S ${FILEHANDLE}.sam

# > Bowtie2 to RM <---------------------------------------------------------
# INPUT: *.R[12].full_trim.fastq.gz
# OUTPUT: *.RM.sam
bowtie2 -t -p 8 --local  -L 20 --ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals -i S,1,0.25 \
	--score-min L,1,0.75  --reorder --no-unal -x indexes/bowtie2_index/RM11a \
	-U ${FILEHANDLE}.R1.full_trim.fastq.gz,${FILEHANDLE}.R2.full_trim.fastq.gz,${FILEHANDLE}.R1.unpaired.full_trim.fastq.gz,${FILEHANDLE}.R2.unpaired.full_trim.fastq.gz \
	-S ${FILEHANDLE}.RM.sam

# > Annotate BY sam <------------------------------------------------------
# INPUT: *.sam
# OUTPUT: *.annotated.sam
batch=2
set=1

python annotate_SAM.py -x 1 -b ${batch} -s ${set} \
	--read 1_trimmed/${FILEHANDLE}_R1_001.trim.fastq.gz \
	--sam ${FILEHANDLE}.sam > ${FILEHANDLE}.annotated.sam

# > Annotate RM sam <----------------------------------------------------------
# INPUT: *.RM.sam
# OUTPUT: *.RM.annotated.sam
batch=2
set=1

python annotate_SAM.py -x 1 -b ${batch} -s ${set} \
	--read 1_trimmed/${FILEHANDLE}_R1_001.trim.fastq.gz \
	--sam ${FILEHANDLE}.RM.sam >${FILEHANDLE}.RM.annotated.sam
mv ${FILEHANDLE}.RM.sam 4_bowtie_RM/

# > Index BY sam file <---------------------------------------------------------
# INPUT: *.annotated.sam
# OUTPUT: *.BY.bam
samtools view -@ $(($SLURM_CPUS_ON_NODE-1)) -b ${FILEHANDLE}.annotated.sam \
	| samtools sort -@ $(($SLURM_CPUS_ON_NODE-1)) -T ${FILEHANDLE} \
	> ${FILEHANDLE}.BY.bam

# > Index RM sam file <---------------------------------------------------------
# INPUT: *.RM.annotated.sam
# OUTPUT: *.RM.bam
samtools view -@ $(($SLURM_CPUS_ON_NODE-1)) -b ${FILEHANDLE}.RM.annotated.sam \
	| samtools sort -@ $(($SLURM_CPUS_ON_NODE-1)) -T ${FILEHANDLE} \
	> ${FILEHANDLE}.RM.bam

# > Generate read counts at each SNP <------------------------------------------
# INPUT: *.BY.bam *.RM.bam
# OUTPUT: ${array_prefix}_genotyping_RM.txt.gz ${array_prefix}_genotyping_BY.txt.gz
array_param=${FILEHANDLE}
array_prefix=${array_param%_*}

python parse_annotated_SAM.py \
	${FILEHANDLE}.BY.bam \
	${FILEHANDLE}.RM.bam

# > Run HMM for genotype inference
python initialize_mmap.py # to initialize the memmaps that will store the reads
python reformat_geno.py # to parse sequencing data into memmap
python hmm_infer_genomes.py # to infer genomes and create recombination map

# END ==========================================================================
