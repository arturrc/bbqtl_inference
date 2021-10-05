# Read in SAM file, and add to the read_ID the well ID (X1B?L?R?)
# First, read in the R2.fastq file
# Parse it, and store as dict the read name.
# Then read the SAM file, parse out the read name, and append to it the well id


import getopt, sys
from pathlib import Path
from collections import OrderedDict
sys.path.append('/n/home00/nnguyenba/lib/python')

import regex
import gzip

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def rc(seq):
	return "".join(complements.get(base, base) for base in reversed(seq))



optlist, args = getopt.getopt(sys.argv[1:],"x:b:s:",["read=","sam="])

batch_ID = ""
cross_ID = ""
read_file_name = ""
sam_file_name = ""
set_ID = ""

for o, a in optlist:
	if o == "-b":
		batch_ID = a
	elif o == "-x":
		cross_ID = a
	elif o == "-s":
		set_ID = a
	elif o == "--read":
		read_file_name = a
	elif o == "--sam":
		sam_file_name = a

if(batch_ID == "" or cross_ID == "" or read_file_name == "" or sam_file_name == "" or set_ID == ""):
	print("Requires a batch ID, a set ID, a cross ID, a read file, and a sam file")
	sys.exit()

# Now deal with the read input file
read_file = Path(read_file_name)
if(not(read_file.is_file())):
	print("Cannot find read file")
	sys.exit()

sam_file = Path(sam_file_name)
if(not(sam_file.is_file())):
	print("Cannot find sam file")
	sys.exit()


current_L = ""
current_R = ""
current_ID = ""

IDs = {}

well_ID = ("@([^\s]+).+?(L[0-9]+R[0-9]+)")
well_ID_re = regex.compile(well_ID)

read_file = gzip.open(read_file_name,"rt")
niterator = 0
nwells = 0
for line in read_file:
	niterator += 1
	if niterator == 1:
		m = well_ID_re.match(line)
		if m:
			current_ID = m.groups()[0]
			current_well_ID = m.groups()[1]
			IDs[current_ID] = "X" + str(cross_ID) + "B" + str(batch_ID) + "S" + str(set_ID) + current_well_ID
			nwells += 1
				
	if niterator == 4:
		niterator = 0
		current_ID = ""
		current_well_ID = ""

read_file.close()

#print(str(nwells))

# Done parsing this.
# Now open sam file and modify the read ID
sam_file = open(sam_file_name,"r")
for line in sam_file:
	if("XS:" in line):
		# Read contains a duplicate position
		continue
	if(line[0:1] == "@"):
		print(line.strip())
		continue

	fields = line.strip().split("	")
	if(fields[1] == "4" or fields[1] == "8"):
		# Read is not aligned
		#print(line)
		continue

	if(fields[0] in IDs):
		print (fields[0] + ":" + IDs[fields[0]] + "	" + '	'.join(fields[1:]))
		#test = 1
	#else:
	#	print (line)

sam_file.close()
