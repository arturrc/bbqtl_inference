#! /usr/bin/python3 -u

# ##############################################################################
# Performs parsing of fastq files to parsed barcoded wgs from the trial run 4 of 
# BBQ.
# Modified by Alex N Nguyen Ba
# nnguyenba@fas.harvard.edu
# 2018
#
# Version 1
#
# LICENCE 
#
# ##############################################################################

# Let python unzip. We then output as a zipped file with the mention "trim" in the file name.
# python parse_genotyping.py --R1=XXXXX.R1.gz --R2=XXXXX.R2.gz
# The strategy is to parse the well ID using the tagmentation barcodes.
# Typical read is:
# barcode - ME:
# We then clip the barcode/ME from the quality and the read itself.
# There is a small chance that the ME is found at the end of the read, we must clip them.
# Use trimmomatic? Try that for now, we'll see later if this is sufficient.

import getopt, sys
from pathlib import Path
from collections import OrderedDict
sys.path.append('/n/home00/nnguyenba/lib/python')

import regex
import gzip


complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))

def make_corrector(options, num_mismatch=0):
	checkers = [regex.compile("("+o+"){e<=" + str(num_mismatch) + "}") for o in options]
	def corrector(match):
		current_error_count = 10000
		current_best_i = ""
		for (i,c) in enumerate(checkers):
			m = c.fullmatch(match)
			if m:
				error_count = m.fuzzy_counts[0] + m.fuzzy_counts[1] + m.fuzzy_counts[2]
				if(error_count == 0):
					return i+1
				elif (error_count < current_error_count and error_count <= num_mismatch):
					current_best_i = i+1
					current_error_count = error_count
		if (current_best_i != ""):
			return current_best_i
		else:
			return 0
	return corrector

def OR(xs):
	return "(" + "|".join(["(?:"+x+")" for x in xs]) + ")"

Tn5Lindex = [
	"CACAGT",
	"GTCGTAT",
	"CCTAAGAC",
	"CGAGTTAGT",
	"ATGTCCTAGC",
	"ACGTGTCTGAT",
	"AGAAGCATGTTG",
	"AGTGAG",
	"CAATGGT",
	"CTGAACTG",
	"AAGCCTTGC",
	"TCGACTGGAC",
	"ACGCAAGAGTT",
	"GTAACCAGAGGT",
	"TGGAAG",
	"CTTGGTT",
	"TCTTCCTG",
	"ATGGTAACG",
	"CAGTCATCAA",
	"TAGTAACGTAC",
	"ATTCACAACACG",
	"CCTTGA",
	"AACTCCA",
	"GCACATTC",
	"AGGCATCTG",
	"GTTCTGTGCA",
	"TAACACCGTCG",
	"ATGCAGCTTGTA",
	"ACTACC",
	"TAACAGG",
	"TGCTACCA",
	"GCAAGTCTA"]

Tn5Rindex = [
	"TGTCTC",
	"CTTGGTC",
	"GCAGGTAA",
	"GACCTACAT",
	"CCGATCTTCT",
	"TGGCCTAGTAT",
	"AGGTGTGAACCT",
	"CCGTAT",
	"GAGGCTT",
	"ATCCTCGC",
	"TTGGAGCTC",
	"ATGCATACTG",
	"TTCTTGCGATT",
	"ACCGTGAAGCAC",
	"TGGAGT",
	"GAATCCT",
	"GATGTCTA",
	"TTGAGAAGG",
	"TCACCTGAGC",
	"GTAGCTCACGA",
	"CACAGCGTTCCA",
	"GGAACT",
	"GCATATG",
	"CTGCGTTC",
	"AAGCGCCAA",
	"TGTGAACGGT",
	"TTCACCACGTC",
	"TCTTCCTGTCTC",
	"CTGCTA",
	"TGCAAGC",
	"ACTTCAGA",
	"AACGTTGCT",
	"GTGCACTCAG|GTGCACCTCC",
	"CAAGGTATGGA",
	"TGCTACGATCTC",
	"AGAGCG",
	"TGATCCA",
	"TACCACCA",
	"CTATGCTCC",
	"GACATGCTAA",
	"TTGTAGGTTAT",
	"CTGAGGCAGTGT",
	"AATACG",
	"TCAGTAC",
	"CACGACCG",
	"AGGATAAGT",
	"TGCCTCCTGA",
	"CCATTAGATCG"]

ME = "AGATGTGTATAAGAGACAG";

correct_L = make_corrector(Tn5Lindex[0:32],1)
correct_R = make_corrector(Tn5Rindex[0:48],1)

read_1 = ("(?e)" + OR(Tn5Lindex[0:32]) + "{e<=1}"
	+ "(?e)" + ME + "{e<=2}")

read_4 = ("(?e)" + OR(Tn5Rindex[0:48]) + "{e<=1}"
	+ "(?e)" + ME + "{e<=2}")

read_1_re = regex.compile(read_1)
read_4_re = regex.compile(read_4)

optlist, args = getopt.getopt(sys.argv[1:],'',['R1=','R2='])

fileR1_name = ""
fileR2_name = ""

fileR1_prefix = ""
fileR2_prefix = ""

for o, a in optlist:
	if o == "--R1":
		fileR1_name = a
	elif o == "--R2":
		fileR2_name = a

if(fileR1_name == "" or fileR2_name == ""):
	print("Requires two input gzipped files")
	sys.exit()

# Check that the files exist
fileR1_file = Path(fileR1_name)
if(not(fileR1_file.is_file())):
	print("Cannot find R1 file")
	sys.exit()

fileR2_file = Path(fileR2_name)
if(not(fileR2_file.is_file())):
	print("Cannot find R2 file")
	sys.exit()

fileR1_prefix = Path(fileR1_name).stem.split('.')[0]
fileR2_prefix = Path(fileR2_name).stem.split('.')[0]

# Should return just the XXXXXX part of XXXXXX.R1.fastq

# Now process the gzip file and loop through
fileR1_file = gzip.open(fileR1_name,"rt")
fileR2_file = gzip.open(fileR2_name,"rt")

# Open the output files.
fileR1_output = gzip.open(fileR1_prefix + ".trim.fastq.gz","wt")
fileR2_output = gzip.open(fileR2_prefix + ".trim.fastq.gz","wt")

niterator = 0
mem_desc_line_1 = ""
mem_desc_line_2 = ""
mem_qual_line_1 = ""
mem_qual_line_2 = ""
mem_read_line_1 = ""
mem_read_line_2 = ""

L_index = ""
R_index = ""

for lineR1, lineR2 in list(zip(fileR1_file, fileR2_file)):
	lineR1 = lineR1.rstrip()
	lineR2 = lineR2.rstrip()

	niterator += 1
	if niterator == 1:
		# These contain the header line.
		mem_desc_line_1 = lineR1
		mem_desc_line_2 = lineR2
		
	if niterator == 2:
		mem_read_line_1 = lineR1
		mem_read_line_2 = lineR2
		
	if niterator == 3:
		continue

	if niterator == 4:
		# These contain the quality scores.
		mem_qual_line_1 = lineR1
		mem_qual_line_2 = lineR2

		# These contain the reads
		niterator = 0

		# Attempt to match the index
		m1 = read_1_re.match(mem_read_line_1)
		if m1:
			L_index = correct_L(m1.groups()[0])
			#print(m1.groups()[0] + "	" + str(R_index))

		m2 = read_4_re.match(mem_read_line_2)
		if m2:
			R_index = correct_R(m2.groups()[0])
			#print(m2.groups()[0] + "	" + str(R_index))

		

		if(L_index != "" and R_index != ""):
			# New description line
			mem_desc_line_1 += ":L" + str(L_index) + "R" + str(R_index)
			mem_desc_line_2 += ":L" + str(L_index) + "R" + str(R_index)

			# New quality scores
			mem_qual_line_1 = mem_qual_line_1[m1.end():]
			mem_qual_line_2 = mem_qual_line_2[m2.end():]

			# Output everything
			fileR1_output.write(mem_desc_line_1 + "\n")
			fileR2_output.write(mem_desc_line_2 + "\n")
			
			fileR1_output.write(mem_read_line_1[m1.end():] + "\n")
			fileR2_output.write(mem_read_line_2[m2.end():] + "\n")			

			fileR1_output.write("+" + "\n")
			fileR2_output.write("+" + "\n")

			fileR1_output.write(mem_qual_line_1 + "\n")
			fileR2_output.write(mem_qual_line_2 + "\n")


		

		# Reset all the memory
		mem_desc_line_1 = ""
		mem_desc_line_2 = ""
		mem_qual_line_1 = ""
		mem_qual_line_2 = ""
		mem_read_line_1 = ""
		mem_read_line_2 = ""
		L_index = ""
		R_index = ""


	