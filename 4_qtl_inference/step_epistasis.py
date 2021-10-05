# Here's an attempt to recode the perl script that threads the QTL finding wrapper into python.
# Instead of having a wrapper to call python scripts, we'll use a single script to launch everything. This avoids having to reparse the data (even though it is fast).

# Ok, so now we're going to try a heuristic to accelerate the QTL addition step.
# The heuristic will be to scan every X QTLs instead of every single one. Once we find a good one, we only scan the x*2 positions around the top hit. I am hoping that this will give at least 2 times faster searches.

import string
import numpy as np
from scipy import linalg
import sys
import csv
import itertools
import time
import random
import argparse
import os
cwd = os.getcwd()
import psutil
process = psutil.Process(os.getpid())

import multiprocessing as mp
from multiprocessing import Pool

#sys.path.append('/n/desai_lab/users/klawrence/BBQ/alldata')
#sys.path.append('/n/home00/nnguyenba/scripts/BBQ/alldata')
try:
	sys.path.append('/n/home00/nnguyenba/scripts/BBQ/alldata')
except:
	sys.path.append('/n/holyscratch01/desai_lab/nnguyenba/BBQ/all_data')
	pass

from spore_defs import *

# Read SNP map
#SNP_reader = csv.reader(open('/n/desai_lab/users/klawrence/BBQ/alldata/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')
#SNP_reader = csv.reader(open('/n/home00/nnguyenba/scripts/BBQ/alldata/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')
SNP_reader = csv.reader(open('/n/holyscratch01/desai_lab/nnguyenba/BBQ/all_data/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')

genome_str = genome_str_to_int(next(SNP_reader))
SNP_list = genome_to_chroms(genome_str)
num_chroms = len(SNP_list)
num_SNPs = [len(x) for x in SNP_list]
num_SNPs_total = sum(num_SNPs)
#print(num_SNPs,file=sys.stdout,flush=True)
#print(num_SNPs_total,file=sys.stdout,flush=True)
chrom_startpoints = get_chrom_startpoints(genome_str)
chrom_endpoints = get_chrom_endpoints(genome_str)

# print(chrom_startpoints) [0, 996, 4732, 5291, 9327, 11187, 12476, 16408, 18047, 20126, 23101, 26341, 30652, 33598, 35398, 39688]
# print(chrom_endpoints) [994, 4730, 5289, 9325, 11185, 12474, 16406, 18045, 20124, 23099, 26339, 30650, 33596, 35396, 39686, 41608]
# print(num_SNPs) [995, 3735, 558, 4035, 1859, 1288, 3931, 1638, 2078, 2974, 3239, 4310, 2945, 1799, 4289, 1921]
#exit()

# Systematically check every positions

from argparse import ArgumentParser, SUPPRESS
# Disable default help
parser = ArgumentParser(add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Add back help 
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)
optional.add_argument('--fit', help='Plain text two-column file containing the fitnesses and the standard errors.', default="/n/holyscratch01/desai_lab/nnguyenba/BBQ/all_data/fitnesses/pheno_data_30C.txt")
optional.add_argument('--log', help='Plain text file logging the progress of the QTL search.', default="epistasis.txt")
optional.add_argument('--oCV', help='Outside cross-validation value (k = 0-9)', type=int, default=0)
optional.add_argument('--iCV', help='Inside cross-validation value (l = 0-8)', type=int, default=0)
optional.add_argument('--model', help='Whether to fit on the training set (m = 0), on the train+test set (m = 1) or on the complete data (m = 2)', type=int, default=0)
optional.add_argument('--dir', help='Directory where intermediate files are found.', default=cwd)
optional.add_argument('--scratch', help='Local scratch directory', default='/n/holyscratch01/desai_lab/nnguyenba/BBQ/all_data/genomes/')
optional.add_argument('--init', help='File of initial positions from single non-epistatic QTLs')
optional.add_argument('--unweighted', help='Only run the forward search on unweighted data.', default=0, type=int)
optional.add_argument('--steps', help='Max number of epistatic terms to search for.', default=300, type=int)
optional.add_argument('--cpu', help='Number of threads to run on.', default=16, type=int)
optional.add_argument('--nosave', help='Set to 1 to avoid saving the npy progress files.', default=0, type=int)
optional.add_argument('--downsample', help='Number of segregants to downsample.', default=0, type=int)
optional.add_argument('--sporelist', help='Restrict searches to a list of spores.')

args = parser.parse_args()

outside_CV = args.oCV # Goes from 0 to 9 # k = 10
inside_CV = args.iCV # Goes from 0 to 8 # l = 9

if(outside_CV > 9 or outside_CV < 0):
	print("--oCV must be [0,9]")
	exit()

if(inside_CV > 8 or inside_CV < 0):
	print("--iCV must be [0,8]")
	exit()

if(~np.isin(args.model , range(3))):
	print("--model must be [0,2]")
	exit()

# Read in the fitness data
fitnesses_data = np.loadtxt(args.fit)
# Parse and see if it has standard errors

if(len(fitnesses_data.shape) != 2 or args.unweighted == 1):
	# No errors found, assume all errors the same.
	if(len(fitnesses_data.shape) == 1):
		fitnesses_data = np.reshape(fitnesses_data,(-1,1))

	fitnesses = fitnesses_data[:,0]
	#fitnesses = np.reshape(fitnesses,(len(fitnesses_data,1)))
	errors = np.ones(len(fitnesses_data))
else:
	fitnesses = fitnesses_data[:,0]
	errors = fitnesses_data[:,1]

errors = np.square(errors)
errors = np.reciprocal(errors)

seed = 100000
np.random.seed(seed) # This allows us to keep the same cross validation sets.

sporelist = np.array(range(len(fitnesses)))
if(args.sporelist):
	sporelist = np.loadtxt(args.sporelist, dtype=int)

if(args.downsample > 0 and args.downsample < len(sporelist)):
	#fitnesses = fitnesses[0:args.downsample]
	#errors = errors[0:args.downsample]
	sporelist = sporelist[0:args.downsample]


# First let's take care of the outside CV
perm = np.random.permutation(sporelist)
train_perm = perm.copy()

if(args.model != 2):
	train_perm = np.delete(train_perm, np.r_[outside_CV/10 * len(sporelist):(outside_CV + 1)/10 * len(sporelist)].astype(int),axis=0)
	validation_perm = np.take(perm, np.r_[outside_CV/10 * len(sporelist):(outside_CV + 1)/10 * len(sporelist)].astype(int))

	if(args.model != 1):
		# Ok now let's take care of the inside CV
		# To do this, we split the train_perm into a train/test permutation
		test_perm = np.take(train_perm, np.r_[inside_CV/9 * len(train_perm):(inside_CV + 1)/9 * len(train_perm)].astype(int))
		train_perm = np.delete(train_perm, np.r_[inside_CV/9 * len(train_perm):(inside_CV + 1)/9 * len(train_perm)].astype(int))


# We're doing a k*l fold validation procedure, where l = k-1.
# This allows us to only create 10 test sets, and only 10 validation sets, so the cross validation loops do not explode.
# For example, let the 80 - 10 - 10 (train - test - validation) split
# We can use the same validation for the following split: 10 - 80 -10 (test - train - validation)
# Now looking at that split, we can use the same test to do the following: 10 - 10 - 80 (test - validation - train)

# We will only 'train' on a subset of the data
train_set = np.take(fitnesses,train_perm) # This is 80% of the fitness data
errors = np.take(errors,train_perm)

phenotypes = train_set[~np.isnan(train_set)] # Is a numpy.ndarray
mean_phenotypes = np.mean(phenotypes)
TSS = np.sum((phenotypes-mean_phenotypes)**2)
errors = errors[~np.isnan(train_set)]
num_usable_spores = len(phenotypes)

# Open all the genotype files
genotypes_file = []
num_lines_genotypes = []
chr_to_scan = []
start = time.clock()
for i in range(16):
	genotypes_file.append(np.load(str(args.scratch) + "/chr"+str(i+1)+"_pos_major.npy", mmap_mode="r")) # Uses 30 gb. Need to load once to cache into memory. Then subsequent searches are near instant.
	#genotypes_file.append(np.load(str(args.scratch) + "/chr"+str(i+1)+"_pos_major.npy"))
	num_lines_genotypes.append(genotypes_file[i].shape[0])
	chr_to_scan.append(i)
	print(str(i) + "	" + str(time.clock() - start) + "	" + str(process.memory_info().rss/1024/1024),file=sys.stderr)

# Here we handle whether the script has to restart or whether we are starting from scratch
current_pos_line = ""
current_beta_line = ""
current_progress_line = ""
flag_refined_pos = 0

geno_file = ""
Q_file = ""
R_file = ""
num_QTLs = 0

if(os.path.isfile(args.dir + "/" + args.log)):
	with open(args.dir + "/" + args.log,'r') as readfile:
		linecount = 0
		for line in readfile:
			line = line.rstrip()
			if(linecount % 4 == 0):
				current_likelihood = line
			elif(linecount % 4 == 1):
				current_pos_line = line
			elif(linecount % 4 == 2):
				current_beta_line = line
			elif(linecount % 4 == 3):
				current_progress_line = line
			linecount = linecount + 1

			if("Done" in line):
				exit()

		# split the progress_line into the relevant flags
		if(linecount > 0):
			arr = current_progress_line.split("\t")
			geno_file = arr[0]
			Q_file = arr[1]
			R_file = arr[2]
			num_QTLs = int(arr[3])


# Read in the file of previous computations if we have found QTLs before. Otherwise, generate them.
prev_pos = []
prev_genotypes = []
prev_pos = np.array(prev_pos)
prev_genotypes = np.array(prev_genotypes)
q = []
r = []
flag_begin = 1

if(num_QTLs != 0):
	# This is restarting, reload all the previous computations
	print(current_pos_line)
	#prev_pos = np.fromstring(current_pos_line, dtype='str',sep="	") # This line includes the epistatic terms, which will be labeled as [x,x]
	prev_pos = np.array(current_pos_line.split("	"))
	# Epistatic terms are ALWAYS after the single effect terms.

	flag_load_prev = 0
	try:
		prev_genotypes = np.load(args.dir + "/" + geno_file)
	except:
		flag_load_prev = 1
		pass

	size_of_prev_genome = (prev_pos.size)
	
	# Consistent prev_pos and prev_genotypes?
	if(flag_load_prev == 1 or prev_genotypes.shape[1] != size_of_prev_genome):
		# Looks like something went wrong. We have to remake the dataframes, we'll base on the positions identified in the previous output file.
		prev_genotypes = np.ones((num_usable_spores,size_of_prev_genome))
		for pos_index in range(len(prev_pos)):
			# is the position an epistatic term?
			if("," in prev_pos[pos_index]):
				# Epistatic term
				# Remove bracket term.
				epistatic_term = prev_pos[pos_index][1:]
				epistatic_term = epistatic_term[:-1]

				# Now split into the two positions
				epistatic_term_arr = epistatic_term.split(",")
				pos_1 = int(epistatic_term_arr[0])
				pos_2 = int(epistatic_term_arr[1])

				chr_qtl_1 = np.searchsorted(np.array(chrom_startpoints), pos_1+0.5)
				chr_qtl_2 = np.searchsorted(np.array(chrom_startpoints), pos_2+0.5)

				start_of_chr_1 = chrom_startpoints[chr_qtl_1-1]
				start_of_chr_2 = chrom_startpoints[chr_qtl_2-1]

				pos_in_chr_1 = pos_1 - start_of_chr_1
				pos_in_chr_2 = pos_2 - start_of_chr_2

				pos_line_1 = genotypes_file[chr_qtl_1-1][pos_in_chr_1]
				pos_line_2 = genotypes_file[chr_qtl_2-1][pos_in_chr_2]

				pos_line_1 = np.take(pos_line_1,train_perm)
				pos_line_2 = np.take(pos_line_2,train_perm)
				
				pos_line_1 = pos_line_1[~np.isnan(train_set)]
				pos_line_2 = pos_line_2[~np.isnan(train_set)]
				
				epis_pos_line = pos_line_1 * pos_line_2
				prev_genotypes[:,pos_index] = epis_pos_line.copy()

			else:
				pos = int(prev_pos[pos_index])
				chr_qtl = np.searchsorted(np.array(chrom_startpoints), pos+0.5)
				start_of_chr = chrom_startpoints[chr_qtl-1]
				pos_in_chr = pos - start_of_chr

				pos_line = genotypes_file[chr_qtl-1][pos_in_chr]
				pos_line = np.take(pos_line, train_perm)
				pos_line = pos_line[~np.isnan(train_set)]
			
				prev_genotypes[:,pos_index] = pos_line.copy()
	
		base_genotypes = np.ones((num_usable_spores,1+size_of_prev_genome))
		base_genotypes[:,1:] = prev_genotypes # First index is the intercept.
		q,r = np.linalg.qr(base_genotypes * np.sqrt(np.reshape(errors,(num_usable_spores,1))))
	else:
		# Ok, so we've successfully saved the previous genotype information. now check if we also successfully saved q,r
		# Do we have q,r?
		flag_remake = 0
		if(os.path.isfile(args.dir + "/" + Q_file) and os.path.isfile(args.dir + "/" + R_file)):
			try:
				q = np.load(args.dir + "/" + Q_file)
			except:
				flag_remake = 1
				pass

			try:
				r = np.load(args.dir + "/" + R_file)
			except:
				flag_remake = 1
				pass
		else:
			flag_remake = 1

		if(flag_remake == 1):
			# Remake
			base_genotypes = np.ones((num_usable_spores,1+size_of_prev_genome))
			base_genotypes[:,1:] = prev_genotypes # First index is the intercept.
			q,r = np.linalg.qr(base_genotypes * np.sqrt(np.reshape(errors,(num_usable_spores,1))))

	flag_begin = 0
	pass

if(flag_begin == 1):
	# Ok, no terms yet. So we'll rebuild all the information from the initial positions
	prev_pos = np.genfromtxt(args.init,dtype='str', max_rows=1)
	size_of_prev_genome = (prev_pos.size)
	# Build the genotype matrix
	prev_genotypes = np.ones((num_usable_spores,size_of_prev_genome))
	for pos_index in range(len(prev_pos)):
		pos = int(prev_pos[pos_index])
		chr_qtl = np.searchsorted(np.array(chrom_startpoints), pos+0.5)
		start_of_chr = chrom_startpoints[chr_qtl-1]
		pos_in_chr = pos - start_of_chr

		pos_line = genotypes_file[chr_qtl-1][pos_in_chr]
		pos_line = np.take(pos_line, train_perm)
		pos_line = pos_line[~np.isnan(train_set)]
			
		prev_genotypes[:,pos_index] = pos_line.copy()
		
	#base_genotypes = np.ones((num_usable_spores,1+size_of_prev_genome))
	#base_genotypes[:,1:] = prev_genotypes # First index is the intercept.
	#q,r = np.linalg.qr(base_genotypes * np.sqrt(np.reshape(errors,(num_usable_spores,1))))


# Let's count the number of single effect qtls
count_epistasis = 0
count_single_effects = 0
for i in range(len(prev_pos)):
	if("," in str(prev_pos[i])):
		count_epistasis = count_epistasis + 1
	else:
		count_single_effects = count_single_effects + 1


# ok we reloaded all the previous computations.
# set up the computation settings:
poolcount = args.cpu*2

def find_epistasis(num):
	lowest_RSS = np.Infinity
	genome_at_lowest_RSS = []
	pos_index_at_lowest_RSS = 0
	last_q = []

	for i in range(count_single_effects-1):
		if(i % poolcount == num): # stupid, but works
			for j in range(i+1,count_single_effects):
				
				if(np.all(np.isin("[" + str(int(prev_pos[i])) + "," + str(int(prev_pos[j])) + "]", prev_pos))):
					# Already found this interaction
					continue;

				epistasis_at_pos_ij = prev_genotypes[:,i] * prev_genotypes[:,j]
				epistasis_at_pos_ij = np.reshape(epistasis_at_pos_ij, (num_usable_spores,1))
				
				WX = epistasis_at_pos_ij * np.sqrt(np.reshape(errors,(num_usable_spores,1)))
				QtX = np.dot(np.transpose(q),WX) # Gets the scale for each vectors in Q.
				QtX_Q = np.einsum('ij,j->i',q,np.ravel(QtX))
				orthogonalized = WX-np.reshape(QtX_Q,(num_usable_spores,1)) # Orthogonalize
				new_q = orthogonalized/np.linalg.norm(orthogonalized) # Orthonormalize
				# This gets the last column of Q.
				# We only need the last column of Q to get the new residuals. We'll assemble the full Q or the full R if we need it (i.e. to obtain betas).

				q_upTy = np.einsum('i,i', np.ravel(new_q), phenotypes * np.sqrt(errors))
				q_upq_upTy = np.ravel(new_q) * q_upTy
				predicted_fitnesses = initial_predicted_fitnesses + q_upq_upTy/np.sqrt(errors)

				# Scale the intercept term
				mean_predicted_fitnesses = np.mean(predicted_fitnesses)
	
				# RSS
				RSS = np.sum((phenotypes - mean_phenotypes - predicted_fitnesses + mean_predicted_fitnesses)**2) # This is the RSS for 1:1 line.
				#print(str(RSS) + "	" + "[" + str(int(prev_pos[i])) + "," + str(int(prev_pos[j])) + "]" + "	" + str(num_usable_spores * math.log(RSS/num_usable_spores)))

				if(RSS < lowest_RSS):
					lowest_RSS = RSS
					genome_at_lowest_RSS = epistasis_at_pos_ij.copy()
					pos_index_at_lowest_RSS = ("[" + str(int(prev_pos[i])) + "," + str(int(prev_pos[j])) + "]")

					last_q = np.ravel(new_q)

	# Now return the values
	return lowest_RSS,genome_at_lowest_RSS,pos_index_at_lowest_RSS, last_q


while(count_epistasis < args.steps):
	
	# Put the intercept
	base_genotypes = np.ones((num_usable_spores,1+size_of_prev_genome))
	base_genotypes[:,1:] = prev_genotypes # First index is the intercept.

	if(flag_begin == 1): # If num_QTLs is not zero, then it was loaded previously.
		q,r = np.linalg.qr(base_genotypes * np.sqrt(np.reshape(errors,(num_usable_spores,1))))

	# Ok, ready to infer!
	start_find_new = time.time()
	# Obtain the initial predicted fitnesses:
	initial_beta = linalg.solve_triangular(r,np.dot(np.transpose(q), phenotypes * np.sqrt(errors)), check_finite=False) # 3.49s for 10000 loops
	# first beta index is the intercept term.
	initial_predicted_fitnesses = np.dot(q,np.dot(r,initial_beta))*1/np.sqrt(errors) # Optimal multiplication order

	mean_initial_predicted_fitnesses = np.mean(initial_predicted_fitnesses)
	
	# RSS
	initial_RSS = np.sum((phenotypes - mean_phenotypes - initial_predicted_fitnesses + mean_initial_predicted_fitnesses)**2) # This is the RSS for 1:1 line.
	# Now search for epistasis
	p = Pool(poolcount)
	results = p.map(find_epistasis, range(poolcount))
	p.close()
	p.join()
	#exit()
	
	# Now parse the results:
	lowest_RSS = np.Infinity
	genome_at_lowest_RSS = []
	pos_index_at_lowest_RSS = ""
	last_q = []

	for i in range(len(results)):
		RSS = results[i][0]
		if(RSS < lowest_RSS):
			lowest_RSS = RSS
			genome_at_lowest_RSS = results[i][1]
			pos_index_at_lowest_RSS = results[i][2] # This is a pair
			last_q = results[i][3]

	likelihood = num_usable_spores * math.log(lowest_RSS/num_usable_spores)

	# Did we improve?
	if(lowest_RSS > initial_RSS):
		# No epistatic component improves fit.
		print("Done finding epistatic terms. (" + str(count_epistasis) + ") . Took : " + str(time.time() - start_find_new) + " seconds. Likelihood: " + str(likelihood), file=sys.stderr)
		# Update the logs
		with open(args.dir + "/" + args.log, "a+") as logfile:
			print("Done", file=logfile)

		exit()


	# update the QR
	q_up = np.zeros([q.shape[0],q.shape[1]+1])
	q_up[:,:-1] = q
	q_up[:,-1] = last_q
	# Compute R
	last_r = np.ravel(np.dot(np.transpose(q_up),genome_at_lowest_RSS * np.sqrt(np.reshape(errors,(num_usable_spores,1)))))
	r_up = np.zeros([r.shape[0]+1,r.shape[1]+1])
	r_up[:-1,:-1] = r
	r_up[:,-1] = last_r
	beta = linalg.solve_triangular(r_up,np.dot(np.transpose(q_up), phenotypes * np.sqrt(errors)), check_finite=False) 

	# Update the logs
	with open(args.dir + "/" + args.log, "a+") as logfile:
		print(likelihood, file=logfile)
		if(size_of_prev_genome > 0):
			print(*prev_pos, sep="\t", end="", file=logfile)
			print("\t", end="", file=logfile)
		print(pos_index_at_lowest_RSS, file=logfile)
		print(*beta[1:len(beta)], file=logfile)
		print("pickle_geno_epis.npy" + "	" + "pickle_q_epis.npy" + "	" + "pickle_r_epis.npy" + "	" + str(num_QTLs + 1), file=logfile)
		
	# Update the values
	# Update the values (remove intercept)
	base_genotypes[:,:-1] = prev_genotypes
	base_genotypes[:,size_of_prev_genome] = np.matrix.flatten(genome_at_lowest_RSS)
	q = q_up
	r = r_up
	if(args.nosave == 0):
		np.save(args.dir + "/" + "pickle_geno_epis", base_genotypes) # No longer has the intercept term
		np.save(args.dir + "/" + "pickle_q_epis", q)
		np.save(args.dir + "/" + "pickle_r_epis", r)

	num_QTLs = num_QTLs + 1
	count_epistasis = count_epistasis + 1
	prev_pos = np.append(prev_pos, pos_index_at_lowest_RSS)
	prev_genotypes = base_genotypes
	size_of_prev_genome = (prev_pos.size)
	
	print("Found new epistasis term (" + str(count_epistasis) + ") . Took : " + str(time.time() - start_find_new) + " seconds. Likelihood: " + str(likelihood), file=sys.stderr)

		
