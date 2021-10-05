import string
import numpy as np
import sys
import csv
import itertools

from vars import *

seq_pool = sys.argv[1]

bc_dict = {}
for i in range(len(timepoints)):
	with open(scratch+'/counts_p'+str(seq_pool)+'_tp'+str(timepoints[i])+'.txt','r') as readfile:
		bc_reader = csv.reader(readfile,delimiter='\t')
		print(timepoints[i])
		for row in bc_reader:
			bc1,bc2,count = row
			bc = bc1+'_'+bc2
			if (bc) in bc_dict:
				bc_dict[(bc)][i] = int(count)
			else:
				bc_dict[(bc)] = np.zeros(len(timepoints),dtype=int)
				bc_dict[(bc)][i] = int(count)
		readfile.close()


writefile = open(scratch+'/joined_counts_p'+str(seq_pool)+'.txt','w')
bc_writer = csv.writer(writefile,delimiter='\t')
bc_writer.writerow(['BC',timepoints[0],timepoints[1],timepoints[2],timepoints[3],timepoints[4]])

for bc,counts in bc_dict.items():
	bc_writer.writerow([bc,counts[0],counts[1],counts[2],counts[3],counts[4]])

writefile.close()
