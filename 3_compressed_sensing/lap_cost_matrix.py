import string
import numpy as np
import sys
import csv
import itertools
from collections import Counter

from scipy.misc import logsumexp
import time


def address_not_empty(bc):
	if bc not in cs_bc_wells:
		return False
	else:
		address = cs_bc_wells[bc]
		if len(address[0]) > 0 and len(address[1]) > 0 and len(address[2]) > 0 and len(address[3]) > 0:
			return True
		else:
			return False
			
def get_counts(bc):
	if bc not in cs_bc_wells:
		return 0
	else:
		address = cs_bc_wells[bc]
		counts = 0
		for i in range(4):
			for j in range(len(address[i])):
				counts += address[i][j][1]
		return counts
		
def get_address(bc):
	address = []
	for i in range(4):
		address.append(cs_bc_wells[bc][i][0][0])
	return(address)

	
def log_sum_of_exp(a,b):
	print(np.abs(b-a))
	
	return a + np.log(1+np.exp(np.abs(b-a)))



def well_to_index(batch,sett,plate,col,row):
	index = (int(batch)-1)*4608+(int(sett)-1)*(384)+(int(plate)-1)*96+(int(row)-1)*12+(int(col)-1)+1
	return(index)

def well_to_index(batch,sett,plate,well):
	row_letter = well[0]
	col = well[1:]
	row_dict = {'A':'1', 'B':'2', 'C':'3', 'D':'4', 'E':'5', 'F':'6', 'G':'7', 'H':'8'}
	if row_letter not in row_dict.keys():
		print(batch,sett,plate,well)
	row = row_dict[row_letter]
	index = (int(batch)-1)*4608+(int(sett)-1)*(384)+(int(plate)-1)*96+(int(row)-1)*12+(int(col)-1)+1
	
	#print(index)
	return(index)
	
def index_to_well(index):
	batch = int(np.floor(float(index)/4608.01)+1)
	remainder = index % 4608
	if remainder == 0: remainder = 4608
	#print(batch,remainder)

	sett = int(np.floor(float(remainder)/384.01)+1)
	remainder = remainder % 384
	if remainder == 0: remainder = 384
	#print(sett,remainder)

	plate = int(np.floor(float(remainder)/96.01)+1)
	remainder = remainder % 96
	if remainder == 0: remainder = 96
	#print(plate,remainder)

	row = int(np.floor(float(remainder)/12.01)+1)
	remainder = remainder % 12
	if remainder == 0: remainder = 12
	#print(row,remainder)

	col = int(remainder)
	#print(col)
	
	# row_dict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8':'H'}
	# well = row_dict[str(row)]+str(col)
	
	return(batch,sett,plate,row,col)


n_batch = 23
n_sett = 12
n_plate = 4
n_row = 8
n_col = 12


# Define the list of wells

# Start with all wells
total_well_number = 23*4608
index_list = np.arange(1,total_well_number+1)
# for i in range(len(index_list)):
# 	well = index_to_well(index_list[i])
# 	print(well[0],well[1],well[2],well[3],'\t',index_list[i])



# Remove wells that were already assigned
indices_to_remove = []
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/first_pass_spore_addresses_5thresh.txt','r') as readfile:
	address_reader = csv.reader(readfile,delimiter='\t')
	for row in address_reader:
		bc1,bc2,batch,sett,plate,well = row
		index = well_to_index(batch,sett,plate,well)
		if index not in indices_to_remove:
			indices_to_remove.append(index)
	readfile.close()

#print(len(indices_to_remove))
# Remove blank wells
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/blanks_all.txt','r') as readfile:
	address_reader = csv.reader(readfile,delimiter='\t')
	header = next(address_reader)
	for row in address_reader:
		if len(row) != 4:
			print(row)
		batch,sett,plate,well = row
		index = well_to_index(batch,sett,plate,well)
		if index not in indices_to_remove:
			indices_to_remove.append(index)
			
	readfile.close()
#print(len(indices_to_remove))	
# Remove wells with duplicate addresses
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/duplicate_addresses_5thresh.txt','r') as readfile:
	address_reader = csv.reader(readfile,delimiter='\t')
	for row in address_reader:
		batch,sett,plate,well,num,list = row
		index = well_to_index(batch,sett,plate,well)
		if index not in indices_to_remove:
			indices_to_remove.append(index)
	readfile.close()

#print(np.sort(np.array(indices_to_remove[:100])))
indices_to_remove = np.sort(np.array(indices_to_remove))
#print(indices_to_remove[:100])

indices_to_match = np.delete(index_list,np.array(indices_to_remove)-1)
print('Wells to match: ',len(indices_to_match),file=sys.stdout,flush=True)
#print(indices_to_match[:100])

# for i in range(len(indices_to_match)):
# 	well = index_to_well(indices_to_match[i])
# 	print(well[0],'\t',well[1],'\t',well[2],'\t',well[3],'\t',indices_to_match[i])
	

# Create masks to add costs to certain wells based on their pools
batch_masks = np.full((23,len(indices_to_match)),0)
sett_masks = np.full((12,len(indices_to_match)),0)
plate_masks = np.full((4,len(indices_to_match)),0)
row_masks = np.full((8,len(indices_to_match)),0)
col_masks = np.full((12,len(indices_to_match)),0)

batch =0
sett =0
plate =0
row =0
col = 0
for i in range(len(indices_to_match)): 
	batch,sett,plate,row,col = index_to_well(indices_to_match[i])
	batch_masks[batch-1][i] = 1
	sett_masks[sett-1][i] = 1
	plate_masks[plate-1][i] = 1
	row_masks[row-1][i] = 1
	col_masks[col-1][i] = 1






# Now, get list of barcodes with their counts
bc_list = []

bc_counts = {}
bc_total_counts = {}
batch_total_counts = np.full(23,0)
sett_total_counts = np.full(12,0)
plate_total_counts = np.full(4,0)
row_total_counts = np.full(8,0)
col_total_counts = np.full(12,0)

batch_counts = np.full(23,0)
sett_counts = np.full(12,0)
plate_counts = np.full(4,0)
row_counts = np.full(8,0)
col_counts = np.full(12,0)


with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/bcs_for_lap_5thresh.txt','r') as readfile:
	bc_reader = csv.reader(readfile,delimiter='\t')
	for row in bc_reader:
		bc1 = row.pop(0)
		bc2 = row.pop(0)
		bc_list.append((bc1,bc2))
		row = [int(x) for x in row]
		batch_counts = np.array(row[0:23])+1
		sett_counts = np.array(row[23:35])+1
		plate_counts = np.array(row[35:39])+1
		row_counts = np.array(row[39:47])+1
		col_counts = np.array(row[47:59])+1
	
			
		if (bc1,bc2) in bc_counts.keys():
			print('Error! double barcode')
			
		bc_counts[(bc1,bc2)] = [batch_counts,sett_counts,plate_counts,row_counts,col_counts]
		bc_total_counts[(bc1,bc2)] = np.sum(batch_counts)+np.sum(sett_counts)+np.sum(plate_counts)+np.sum(row_counts)+np.sum(col_counts)
		
		batch_total_counts += batch_counts
		sett_total_counts += sett_counts
		plate_total_counts += plate_counts
		row_total_counts += row_counts
		col_total_counts += col_counts
		

	readfile.close()

print(batch_total_counts)
print(sett_total_counts)
print(plate_total_counts)
print(row_total_counts)
print(col_total_counts)

print('Barcodes to match: ',str(len(bc_counts)),file=sys.stdout,flush=True)
	

# Now, calculate the cost matrix

cost_matrix = np.empty((len(bc_counts),len(indices_to_match)))
for i in range(len(bc_list)):
	
	bc = bc_list[i]
	#print(bc,file=sys.stdout,flush=True)
	
	#print(bc_counts[bc])
	
	# batch posterior
	posterior_batch = np.empty(n_batch)
	for j in range(n_batch):		
		n1 = bc_counts[bc][0][j] # counts for this barcode in batch pool j
		N1 = batch_total_counts[j] # total counts in batch pool j
		p1 = float(n1)/float(N1)
		sum_n2 = np.sum(bc_counts[bc][0])-n1 # counts for this barcode in other batch pools
		sum_N2 = np.sum(batch_total_counts)-N1 # total counts in other batch pools
		p2 = float(sum_n2)/float(sum_N2)
		if p2 > p1: p2 = p1
		posterior_batch[j] = n1*np.log(p1) + (N1-n1)*np.log(1-p1) + sum_n2*np.log(p2) + (sum_N2-sum_n2)*np.log(1-p2)		
	ll_sum_batch = logsumexp(posterior_batch) # total post for normalization
	cost_batch = -1*(posterior_batch - ll_sum_batch) # normalized cost for each batch pool	
	#print('batch: ',cost_batch,file=sys.stdout,flush=True)
	cost_vector = np.tensordot(cost_batch,batch_masks,axes=1)
	cost_matrix[i] += cost_vector
	
	
	# sett posterior
	posterior_sett = np.empty(n_sett)
	for j in range(n_sett):		
		n1 = bc_counts[bc][1][j] # counts for this barcode in sett pool j
		N1 = sett_total_counts[j] # total counts in sett pool j
		p1 = float(n1)/float(N1)
		sum_n2 = np.sum(bc_counts[bc][1])-n1 # counts for this barcode in other sett pools
		sum_N2 = np.sum(sett_total_counts)-N1 # total counts in other sett pools
		p2 = float(sum_n2)/float(sum_N2)
		if p2 > p1: p2 = p1
		posterior_sett[j] = n1*np.log(p1) + (N1-n1)*np.log(1-p1) + sum_n2*np.log(p2) + (sum_N2-sum_n2)*np.log(1-p2)		
	ll_sum_sett = logsumexp(posterior_sett) # total post for normalization
	cost_sett = -1*(posterior_sett - ll_sum_sett) # normalized cost for each sett pool	
	cost_vector = np.tensordot(cost_sett,sett_masks,axes=1)
	#print('sett: ',cost_sett,file=sys.stdout,flush=True)
	cost_matrix[i] += cost_vector
	
	
	# plate posterior
	posterior_plate = np.empty(n_plate)
	for j in range(n_plate):		
		n1 = bc_counts[bc][2][j] # counts for this barcode in plate pool j
		N1 = plate_total_counts[j] # total counts in plate pool j
		p1 = float(n1)/float(N1)
		sum_n2 = np.sum(bc_counts[bc][2])-n1 # counts for this barcode in other plate pools
		sum_N2 = np.sum(plate_total_counts)-N1 # total counts in other plate pools
		p2 = float(sum_n2)/float(sum_N2)
		if p2 > p1: p2 = p1
		posterior_plate[j] = n1*np.log(p1) + (N1-n1)*np.log(1-p1) + sum_n2*np.log(p2) + (sum_N2-sum_n2)*np.log(1-p2)		
	ll_sum_plate = logsumexp(posterior_plate) # total post for normalization
	cost_plate = -1*(posterior_plate - ll_sum_plate) # normalized cost for each plate pool	
	cost_vector = np.tensordot(cost_plate,plate_masks,axes=1)
	#print('plate: ',cost_plate,file=sys.stdout,flush=True)	
	cost_matrix[i] += cost_vector
	
	
	# row posterior
	posterior_row = np.empty(n_row)
	for j in range(n_row):		
		n1 = bc_counts[bc][3][j] # counts for this barcode in row pool j
		N1 = row_total_counts[j] # total counts in row pool j
		p1 = float(n1)/float(N1)
		sum_n2 = np.sum(bc_counts[bc][3])-n1 # counts for this barcode in other row pools
		sum_N2 = np.sum(row_total_counts)-N1 # total counts in other row pools
		p2 = float(sum_n2)/float(sum_N2)
		if p2 > p1: p2 = p1
		posterior_row[j] = n1*np.log(p1) + (N1-n1)*np.log(1-p1) + sum_n2*np.log(p2) + (sum_N2-sum_n2)*np.log(1-p2)		
	ll_sum_row = logsumexp(posterior_row) # total post for normalization
	cost_row = -1*(posterior_row - ll_sum_row) # normalized cost for each row pool	
	cost_vector = np.tensordot(cost_row,row_masks,axes=1)
	#print('row: ',cost_row,file=sys.stdout,flush=True)	
	cost_matrix[i] += cost_vector
	
	
	# col posterior
	posterior_col = np.empty(n_col)
	for j in range(n_col):		
		n1 = bc_counts[bc][4][j] # counts for this barcode in col pool j
		N1 = col_total_counts[j] # total counts in col pool j
		p1 = float(n1)/float(N1)
		sum_n2 = np.sum(bc_counts[bc][4])-n1 # counts for this barcode in other col pools
		sum_N2 = np.sum(col_total_counts)-N1 # total counts in other col pools
		p2 = float(sum_n2)/float(sum_N2)
		if p2 > p1: p2 = p1
		posterior_col[j] = n1*np.log(p1) + (N1-n1)*np.log(1-p1) + sum_n2*np.log(p2) + (sum_N2-sum_n2)*np.log(1-p2)		
	ll_sum_col = logsumexp(posterior_col) # total post for normalization
	cost_col = -1*(posterior_col - ll_sum_col) # normalized cost for each col pool	
	cost_vector = np.tensordot(cost_col,col_masks,axes=1)
	#print('col: ',cost_col,file=sys.stdout,flush=True)	
	cost_matrix[i] += cost_vector
	
print('Done calculating cost matrix.',file=sys.stdout,flush=True)


# Now store cost matrix & lists
with open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/lap_cost_matrix_5thresh_bayes.txt','w') as writefile:
	cost_writer = csv.writer(writefile,delimiter='\t')
	cost_writer.writerow(bc_list)
	cost_writer.writerow(indices_to_match)
	for i in range(len(cost_matrix)):
		cost_writer.writerow(cost_matrix[i])

	writefile.close()