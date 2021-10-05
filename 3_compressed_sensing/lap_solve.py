import string
import numpy as np
import sys
import csv
import itertools
import lap


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
	
	row_dict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', '7':'G', '8':'H'}
	well = row_dict[str(row)]+str(col)
	
	return(batch,sett,plate,well)

	

##Read lists & cost matrix
with open('/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/cs/count_files/lap_cost_matrix_5thresh_bayes.txt','r') as readfile:
	cost_reader = csv.reader(readfile,delimiter='\t')
	bc_list = next(cost_reader)
	index_list = next(cost_reader)

	cost_matrix = np.empty((len(bc_list),len(index_list)))
	for i in range(len(cost_matrix)):
		cost_matrix[i] = next(cost_reader)

	readfile.close()

print('Finished reading cost matrix.',file=sys.stdout,flush=True)



##Solve!
cost,x,y = lap.lapjv(cost_matrix,extend_cost=True)

print('Finished solving.',file=sys.stdout,flush=True)

##Print well assigments to file
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/lap_assignments_5thresh_bayes.txt','w') as writefile:
	well_writer = csv.writer(writefile,delimiter='\t')
	
	if len(bc_list) < len(index_list):
		for i in range(len(x)):
			bc = bc_list[i]
			index = index_list[x[i]]
			well = index_to_well(int(index))
			bc1 = bc[2:18]
			bc2 = bc[22:38]
			well_writer.writerow([bc1,bc2,well[0],well[1],well[2],well[3]])
			
		
	elif len(index_list) < len(bc_list):
		for j in range(len(y)):
			index = index_list[j]
			bc = bc_list[y[j]]
			well = index_to_well(int(index))
			bc1 = bc[2:18]
			bc2 = bc[22:38]
			well_writer.writerow([bc1,bc2,well[0],well[1],well[2],well[3]])

	writefile.close()
