import string
import numpy as np
import sys
import csv
import itertools
import gzip
import re

sys.path.append('/n/desai_lab/users/klawrence/BBQ/alldata')
from spore_defs import *


def get_address(location,tagset):
	# parses Alex's location info into address format
	address = "X(\d+)B(\d+)S(\d+)L(\d+)R(\d+)"
	address_re = re.compile(address)
	m1 = address_re.match(location)
	if m1:
		batch = int(m1.groups()[1])
		set = int(m1.groups()[2])
		left = int(m1.groups()[3])
		right = int(m1.groups()[4])
				
		if tagset == 'AA':
			if left > 16 or right > 24: 
				return([])				
			if left <= 8:
				if right <= 12: plate = 1
				elif right > 12: plate = 2
			elif left > 8:
				if right <= 12: plate = 3
				elif right > 12: plate = 4
		elif tagset == 'AB':
			if left > 16 or right <= 24: 
				return([])				
			if left <= 8:
				if right <= 36: plate = 1
				elif right > 36: plate = 2
			elif left > 8:
				if right <= 36: plate = 3
				elif right > 36: plate = 4
		elif tagset == 'BA':
			if left <= 16 or right > 24: 
				return([])				
			if left <= 24:
				if right <= 12: plate = 1
				elif right > 12: plate = 2
			elif left > 24:
				if right <= 12: plate = 3
				elif right > 12: plate = 4		
		elif tagset == 'BB':
			if left <= 16 or right <= 24: 
				return([])				
			if left <= 24:
				if right <= 36: plate = 1
				elif right > 36: plate = 2
			elif left > 24:
				if right <= 36: plate = 3
				elif right > 36: plate = 4
		
		
		well_dict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8':'H'}
		
		left_adjusted = left-8*np.int((left-1)/8)
		right_adjusted = right-12*np.int((right-1)/12)
		well = well_dict[str(left_adjusted)] + str(right_adjusted)
					
		return([batch,set,plate,well])
	else:
		return([])
	


#######	

# Read SNP map
SNP_reader = csv.reader(open('/n/desai_lab/users/klawrence/BBQ/alldata/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')
SNP_list = genome_to_chroms(genome_str_to_int(next(SNP_reader)))
num_chroms = len(SNP_list)
num_SNPs = [len(x) for x in SNP_list]
num_SNPs_total = sum(num_SNPs)
print(num_SNPs,file=sys.stdout,flush=True)
print(num_SNPs_total,file=sys.stdout,flush=True)

# Define memmap to store read counts
num_spores = 384*12*23
spore_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_novaseq_spores_reads',dtype=[('address','U3',4),('reads_RM','f8',num_SNPs_total+15),('reads_BY','f8',num_SNPs_total+15) ],mode='r+',shape=num_spores)

# Find the spore indices relevant for this batch
batch = int(sys.argv[1]) # batch in range 2 to 23
section_start = 384*12*(batch-1)
section_stop = 384*12*batch

this_spore = 384*12*(batch-1)
num_both_reads = 0

# set of tagmentase indices to use for each batch
tagsets = ['AA','AB','AA','BB','AA','AB','BA','AA','BB','AB','BA','AA','BB','AA','BB','AB','AB','BA','AA','AA','BB','AB','BA']
tagset = tagsets[batch-1]

# Loop over sets in each batch, writing address & genome for each spore
print('Batch '+str(batch),file=sys.stdout,flush=True)	
for i in range(1,13):
	print('Set '+str(i),file=sys.stdout,flush=True)
	addresses = {}
	
	# RM reads
	with gzip.open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/data_novaseq/geno_B'+str(batch)+'T'+str(i)+'_genotyping_RM.txt.gz','rt') as readfile:	
		spore_reader = csv.reader(readfile,delimiter='\t')
		header = next(spore_reader)
		split_indices = [l for l in range(len(header)) if header[l] == 'DUMMY'][:15]
		#print(split_indices)
		for row in spore_reader:
			name = row.pop(0)
			address = get_address(name,tagset)
				
			if len(address) > 0:
				addresses[tuple(address)] = this_spore
				spore_mmap[this_spore]['address'] = address
				genome = genome_str_to_float(row)[:num_SNPs_total+15]
				for index in split_indices:
					genome[index-1] = -1.0
				spore_mmap[this_spore]['reads_RM'] = genome
				this_spore += 1
				
		readfile.close()

	# BY reads
	with gzip.open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/data_novaseq/geno_B'+str(batch)+'T'+str(i)+'_genotyping_BY.txt.gz','rt') as readfile:	
		spore_reader = csv.reader(readfile,delimiter='\t')
		header = next(spore_reader)
		split_indices = [l for l in range(len(header)) if header[l] == 'DUMMY'][:15]
		#print(split_indices)
		for row in spore_reader:
			name = row.pop(0)
			address = get_address(name,tagset)
			if len(address) > 0 and tuple(address) in addresses.keys():
				sporenum = addresses[tuple(address)]
				genome = genome_str_to_float(row)[:num_SNPs_total+15]
				for index in split_indices:
					genome[index-1] = -1.0
				spore_mmap[sporenum]['reads_BY'] = genome
				num_both_reads += 1
		readfile.close()
		

spore_mmap.flush()
print(section_start,section_stop,file=sys.stdout,flush=True)
print(this_spore,num_both_reads,file=sys.stdout,flush=True)

