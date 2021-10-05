import getopt, sys
from pathlib import Path
from collections import OrderedDict
sys.path.append('/n/home00/nnguyenba/lib/python')

import regex
import pysam
import gzip



SNPdict = {}
genotyping = {}
total_counts = {}

# Parse the SNP file
SNPfile = open("/n/desai_lab/users/nnguyenba/BBQ/pilot_1/anc_reads/RM_snps.gd", "r")

for i in range(17):
	SNPdict["chr" + str(i+1).zfill(2)] = {}
	genotyping["chr" + str(i+1).zfill(2)] = {}
	total_counts["chr" + str(i+1).zfill(2)] = {}

for line in SNPfile:
	fields = line.rstrip().split("	")
	if(fields[0] == "SNP"):
		# Convert to 0 index for pysam
		SNPdict[fields[1]][int(fields[2])-1] = fields[3]
		genotyping[fields[1]][int(fields[2])-1] = {}
		total_counts[fields[1]][int(fields[2])-1] = {}

	if(fields[0] == "SUB"):
		if(len(fields[4]) == int(fields[3])):
			for i in range(0, int(fields[3])):
				SNPdict[fields[1]][int(fields[2])+i-1] = fields[4][i]
				genotyping[fields[1]][int(fields[2])+i-1] = {}
				total_counts[fields[1]][int(fields[2])+i-1] = {}

		

# Open the reference file
BYref = open("/n/home00/nnguyenba/genome_references/bowtie2_index/BY4742_fixed.fasta", "r")
current_ID = ""
references = {}
for line in BYref:
	field = line.rstrip()
	if(len(field) > 0):
		if(field[0] == ">"):
			current_ID = field[1:]
			references[current_ID] = ""
		else:
			references[current_ID] += field

wells = {}

BY_samfile = pysam.AlignmentFile(sys.argv[1])
RM_samfile = pysam.AlignmentFile(sys.argv[2])

samfile_prefix = Path(sys.argv[1]).stem.split('.')[0]

samfile_prefix = samfile_prefix.split("_S")[0]

# Loop through all the SNPs
# It would be best to print to stderr some progress. This script is slow.
for chr in SNPdict:
	if(chr == "chr17"):
		continue
	for pos in SNPdict[chr]:
		print (chr + "	" + str(pos) + "	" + SNPdict[chr][pos], file=sys.stderr)
		
		for (pileupcolumn, pileupcolumn_RM) in zip(BY_samfile.pileup(chr, start=pos, stop=pos+1, truncate=1), RM_samfile.pileup(chr, start=pos, stop=pos+1, truncate=1)):
			#print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
			#print("Reference base at pos: " + str(pileupcolumn.pos) + " is " + references[chr][pileupcolumn.pos])

			# Pileupcolumn is the array of reads we'll be analyzing.
			# Only analyze reads that also show up in pileupcolumn_RM
			RM_query_names = pileupcolumn_RM.get_query_names()
			#print(RM_query_names)
			#print(pileupcolumn.get_query_names())
			for pileupread in pileupcolumn.pileups:
				# query position is None if is_del or is_refskip is set.
				
				if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_name in RM_query_names:
					# Grab the well id
					fields = pileupread.alignment.query_name.split(":")
					well_ID = fields[-1]
					wells[well_ID] = 1
					
					base_in_read = pileupread.alignment.query_sequence[pileupread.query_position]
					
					if(base_in_read == references[chr][pos]):
						# genotype to BY
						genotyping[chr][pos][well_ID] = genotyping[chr][pos].get(well_ID,0) + 0
						total_counts[chr][pos][well_ID] = total_counts[chr][pos].get(well_ID,0) + 1
					elif(base_in_read == SNPdict[chr][pos]):
						# genotype to RM
						genotyping[chr][pos][well_ID] = genotyping[chr][pos].get(well_ID,0) + 1
						total_counts[chr][pos][well_ID] = total_counts[chr][pos].get(well_ID,0) + 1
					
					#print ('\tbase in read %s = %s' %
					#	(pileupread.alignment.query_name,
					#	pileupread.alignment.query_sequence[pileupread.query_position]))
			#exit()

# Done getting all the genotypes
# now list in a nice table format
# Transpose of the previous format

# Open the output files.
fileBY_output = gzip.open(samfile_prefix + "_genotyping_BY.txt.gz","wt")
fileRM_output = gzip.open(samfile_prefix + "_genotyping_RM.txt.gz","wt")

# Print header
for chr in sorted(SNPdict):
	for pos in sorted(SNPdict[chr]):
		fileBY_output.write("	" + chr + "_" + references[chr][pos] + "_" + str(pos+1) + "_" + SNPdict[chr][pos])
		fileRM_output.write("	" + chr + "_" + references[chr][pos] + "_" + str(pos+1) + "_" + SNPdict[chr][pos])
		#print ("	" + chr + "_" + references[chr][pos] + "_" + str(pos+1) + "_" + SNPdict[chr][pos], end='')
	#print("	DUMMY", end='')
	fileBY_output.write("	DUMMY")
	fileRM_output.write("	DUMMY")

#print()
fileBY_output.write("\n")
fileRM_output.write("\n")

# Print data
for well in wells:
	#print(well, end='')
	fileBY_output.write(well)
	fileRM_output.write(well)

	for chr in sorted(SNPdict):
		for pos in sorted(SNPdict[chr]):
			if(well in genotyping[chr][pos]):
				#print("	" + str(genotyping[chr][pos][well]/total_counts[chr][pos][well]), end='')
				fileBY_output.write("	" + str(total_counts[chr][pos][well] - genotyping[chr][pos][well]))
				fileRM_output.write("	" + str(genotyping[chr][pos][well]))
			else:
				#print("	nan", end='')
				fileBY_output.write("	nan")
				fileRM_output.write("	nan")
		#print("	" + str(2.0), end='')
		fileBY_output.write("	" + str(2.0))
		fileRM_output.write("	" + str(2.0))
	print()
	fileBY_output.write("\n")
	fileRM_output.write("\n")

sys.exit()


