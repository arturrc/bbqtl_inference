import string
import numpy as np
import sys
import csv
import itertools
from collections import Counter

def count_list_not_empty(bc):
	if bc not in bc_pool_counts:
		return False
	else:
		address = bc_pool_counts[bc]
		if len(address[0]) > 0 and len(address[1]) > 0 and len(address[2]) > 0 and len(address[3]) > 0 and len(address[4]) > 0:
			return True
		else:
			return False
			
def get_counts(bc):
	if bc not in bc_pool_counts:
		return 0
	else:
		address = bc_pool_counts[bc]
		counts = 0
		for i in range(5):
			for j in range(len(address[i])):
				counts += address[i][j][1]
		return counts
		
def get_address(bc):
	address = []
	for i in range(5):
		address.append(bc_pool_counts[bc][i][0][0])
	return(address)

def get_num_pools(bc):
	address = bc_pool_counts[bc]
	num_pools = [0,0,0,0,0]
	for i in range(5):
		num_pools[i] = len(address[i])
	return num_pools
	
def get_best_pool(bc):
	address = bc_pool_counts[bc]
	#print(address)
	for i in range(5):
		if len(address[i]) > 1:
			counts = []
			pools = []
			for j in range(len(address[i])):
				counts.append(address[i][j][1])
				pools.append(address[i][j][0])
			#print(counts)
			maxcount = np.max(np.array(counts))
			maxpool = pools[np.argmax(np.array(counts))]
			#print(maxpool)
			address[i] = [(maxpool,maxcount)]
	#print(address)
	return address

def rowcol_to_well(row,col):
	row_names = ['A','B','C','D','E','F','G','H']
	well = row_names[int(row)-1]+str(col)
	return well

bc_pool_counts = {}
n_pcol = 12
n_prow = 8
n_batch = 23
n_set = 12
n_plate = 4
threshold = 3

for i in range(1,n_batch+1):
	batch_reader = csv.reader(open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/counts_batch_'+str(i)+'.txt','r'))
	for row in batch_reader:
		bc1,bc2,count = row[0].split()
		count = int(count)
		if count >= threshold:
			if (bc1,bc2) not in bc_pool_counts:
				bc_pool_counts[(bc1,bc2)] = [[(i,count)],[],[],[],[]]
			else:
				bc_pool_counts[(bc1,bc2)][0].append((i,count))

for i in range(1,n_set+1):
	set_reader = csv.reader(open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/counts_set_'+str(i)+'.txt','r'))
	for row in set_reader:
		bc1,bc2,count = row[0].split()
		count = int(count)
		if count >= threshold:
			if (bc1,bc2) not in bc_pool_counts:
				bc_pool_counts[(bc1,bc2)] = [[],[(i,count)],[],[],[]]
			else:
				bc_pool_counts[(bc1,bc2)][1].append((i,count))

for i in range(1,n_plate+1):
	plate_reader = csv.reader(open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/counts_plate_'+str(i)+'.txt','r'))
	for row in plate_reader:
		bc1,bc2,count = row[0].split()
		count = int(count)
		if count >= threshold:
			if (bc1,bc2) not in bc_pool_counts:
				bc_pool_counts[(bc1,bc2)] = [[],[],[(i,count)],[],[]]
			else:
				bc_pool_counts[(bc1,bc2)][2].append((i,count))

for i in range(1,n_prow+1):
	prow_reader = csv.reader(open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/counts_prow_'+str(i)+'.txt','r'))
	for row in prow_reader:
		bc1,bc2,count = row[0].split()
		count = int(count)
		if count >= threshold:
			if (bc1,bc2) not in bc_pool_counts:
				bc_pool_counts[(bc1,bc2)] = [[],[],[],[(i,count)],[]]
			else:
				bc_pool_counts[(bc1,bc2)][3].append((i,count))

for i in range(1,n_pcol+1):
	pcol_reader = csv.reader(open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/counts_pcol_'+str(i)+'.txt','r'))
	for row in pcol_reader:
		bc1,bc2,count = row[0].split()
		count = int(count)
		if count >= threshold:
			if (bc1,bc2) not in bc_pool_counts:
				bc_pool_counts[(bc1,bc2)] = [[],[],[],[],[(i,count)]]
			else:
				bc_pool_counts[(bc1,bc2)][4].append((i,count))


print('Finished reading files.')

n=0
n_above_threshold = 0
n_not_empty = 0
n_unique = 0
n_unique_post = 0
n_took_best_pool = 0
address_dict = {}
lap_bcs = {}

for bc in bc_pool_counts.keys():
	n += 1
	count_list = bc_pool_counts[bc]
	
	if get_counts(bc) > threshold:
		n_above_threshold += 1
	
		if count_list_not_empty(bc):
			#print(bc_pool_counts[bc])
			n_not_empty += 1
		
			num_pools = get_num_pools(bc)
		
			if num_pools == [1,1,1,1,1]:
				n_unique += 1
				final_count_list = str(count_list[0][0][0])+'\t'+str(count_list[1][0][0])+'\t'+str(count_list[2][0][0])+'\t'+str(count_list[3][0][0])+'\t'+str(count_list[4][0][0])
				if final_count_list not in address_dict.keys():
					address_dict[final_count_list] = [bc]
				else:
					address_dict[final_count_list].append(bc)
		
			elif sum(num_pools) > 5:
				
				lap_bcs[bc] = count_list

				# new_count_list = get_best_pool(bc)
				# final_count_list = str(new_count_list[0][0][0])+'\t'+str(new_count_list[1][0][0])+'\t'+str(new_count_list[2][0][0])+'\t'+str(new_count_list[3][0][0])+'\t'+str(new_count_list[4][0][0])
				# n_took_best_pool += 1
				# if final_count_list not in address_dict.keys():
				# 	address_dict[final_count_list] = [bc]
				# else:
				# 	address_dict[final_count_list].append(bc)

		elif sum(get_num_pools(bc)) >= 3:
			lap_bcs[bc] = count_list

				
		
print('Total barcodes: '+str(n))
print('Barcodes above '+str(threshold)+' reads: '+str(n_above_threshold))
print('Barcodes with no missing pools: '+str(n_not_empty))
print('Barcodes with unique addresses: '+str(n_unique))
print('Barcodes where best pool was taken: '+str(n_took_best_pool))
print('Total located barcodes: '+str(len(address_dict)))

n_duplicates = 0
duplicates = []
for address,bc_list in address_dict.items():
	if len(bc_list) > 1:
		n_duplicates += 1
		duplicates.append(address)
				
		#for j in range(len(bc_list)):
			#lap_bcs[bc_list[j]] = bc_pool_counts[bc_list[j]]

print('Overlapping assignments: '+str(n_duplicates))
final_bcs = {k:v for k,v in address_dict.items() if k not in duplicates}
print('Final fist-pass barcodes: '+str(len(final_bcs)))
print('Barcodes passed to LAP solver: '+str(len(lap_bcs)))

duplicate_bcs = {k:v for k,v in address_dict.items() if k in duplicates}
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/duplicate_addresses.txt','w') as writefile:
	loc_writer = csv.writer(writefile,delimiter='\t')
	for address,bc_list in duplicate_bcs.items():
		address = address.split('\t')
		well = rowcol_to_well(address[3],address[4])
		loc_writer.writerow([address[0],address[1],address[2],well,len(bc_list),bc_list])
writefile.close()


with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/first_pass_spore_addresses.txt','w') as writefile:
	loc_writer = csv.writer(writefile,delimiter='\t')
	for address,bc in final_bcs.items():
		bc1 = bc[0][0]
		bc2 = bc[0][1]
		address = address.split('\t')
		well = rowcol_to_well(address[3],address[4])
		loc_writer.writerow([bc1,bc2,address[0],address[1],address[2],well])
writefile.close()

with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/bcs_for_lap.txt','w') as writefile:
	loc_writer = csv.writer(writefile,delimiter='\t')
	for bc,count_list in lap_bcs.items():
		bc1 = bc[0]
		bc2 = bc[1]

		total_pools = n_pcol+n_prow+n_batch+n_set+n_plate
		pool_counts = list(np.full(total_pools,0))

		batch_list = count_list[0]
		if len(batch_list) > 0:
			for i in range(len(batch_list)):
				pool_counts[batch_list[i][0]-1] = batch_list[i][1]
				
		set_list = count_list[1]
		if len(set_list) > 0:
			for i in range(len(set_list)):
				pool_counts[23+set_list[i][0]-1] = set_list[i][1]
				
		plate_list = count_list[2]
		if len(plate_list) > 0:
			for i in range(len(plate_list)):
				pool_counts[23+12+plate_list[i][0]-1] = plate_list[i][1]

		prow_list = count_list[3]
		if len(prow_list) > 0:
			for i in range(len(prow_list)):
				pool_counts[23+12+4+prow_list[i][0]-1] = prow_list[i][1]
				
		pcol_list = count_list[4]
		if len(pcol_list) > 0:
			for i in range(len(pcol_list)):
				pool_counts[23+12+4+8+pcol_list[i][0]-1] = pcol_list[i][1]
		
		pool_counts.insert(0,bc2)
		pool_counts.insert(0,bc1)

		loc_writer.writerow(pool_counts)
writefile.close()

