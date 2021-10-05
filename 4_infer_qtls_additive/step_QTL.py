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
required.add_argument('--fit', help='Plain text two-column file containing the fitnesses and the standard errors.')
optional.add_argument('--log', help='Plain text file logging the progress of the QTL search.', default="output.txt")
optional.add_argument('--oCV', help='Outside cross-validation value (k = 0-9)', type=int, default=0)
optional.add_argument('--iCV', help='Inside cross-validation value (l = 0-8)', type=int, default=0)
optional.add_argument('--model', help='Whether to fit on the training set (m = 0), on the train+test set (m = 1) or on the complete data (m = 2)', type=int, default=0)
optional.add_argument('--dir', help='Directory where intermediate files are found.', default=cwd)
optional.add_argument('--scratch', help='Local scratch directory', default='/n/holyscratch01/desai_lab/nnguyenba/BBQ/all_data/genomes/')
optional.add_argument('--refine', help='Refine every X QTLs, default is 5. 0 means never refine.', default=5, type=int)
optional.add_argument('--unweighted', help='Only run the forward search on unweighted data.', default=0, type=int)
optional.add_argument('--cpu', help='Number of threads to run on.', default=16, type=int)
optional.add_argument('--nosave', help='Set to 1 to avoid saving the npy progress files.', default=0, type=int)
optional.add_argument('--maxqtl', help='Number of QTLs to find.', default=300, type=int)
optional.add_argument('--downsample', help='Number of segregants to downsample.', default=0, type=int)
optional.add_argument('--sporelist', help='Restrict searches to a list of spores.')

args = parser.parse_args()

print(args, file=sys.stderr)

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

if(args.refine == 0):
	args.refine = np.Infinity

# Read in the fitness data
fitnesses_data = np.loadtxt(args.fit)
# Parse and see if it has standard errors

if(len(fitnesses_data.shape) != 2 or args.unweighted == 1):
	# No errors found, assume all errors the same.

	if(len(fitnesses_data.shape) == 1):
		fitnesses_data = np.reshape(fitnesses_data,(-1,1))

	fitnesses = fitnesses_data[:,0]
	errors = np.ones(len(fitnesses_data))
else:
	fitnesses = fitnesses_data[:,0]
	errors = fitnesses_data[:,1]

errors = np.square(errors)
errors = np.reciprocal(errors)


seed = 100000
np.random.seed(seed) # This allows us to keep the same cross validation sets.

# If we are restricting search to a list of spores, then need to parse the list of spores.
sporelist = np.array(range(len(fitnesses)))
if(args.sporelist):
	sporelist = np.loadtxt(args.sporelist, dtype=int)

# First let's take care of the outside CV

if(args.downsample > 0 and args.downsample < len(sporelist)):
	#fitnesses = fitnesses[0:args.downsample]
	#errors = errors[0:args.downsample]
	sporelist = sporelist[0:args.downsample]

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
start = time.perf_counter()
for i in range(16):
	#genotypes_file.append(np.load(str(args.scratch) + "/chr"+str(i+1)+"_pos_major.npy", mmap_mode="r")) # Uses 30 gb. Need to load once to cache into memory. Then subsequent searches are near instant.
	genotypes_file.append(np.load(str(args.scratch) + "/chr"+str(i+1)+"_pos_major.npy"))
	num_lines_genotypes.append(genotypes_file[i].shape[0])
	chr_to_scan.append(i)
	print(str(i) + "	" + str(time.perf_counter() - start) + "	" + str(process.memory_info().rss/1024/1024),file=sys.stderr)


# Here we will handle whether the script has been restart or whether we are starting from scratch.
# Open the log file.
current_likelihood = np.Infinity
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

		# split the progress_line into the relevant flags
		if(linecount > 0):
			arr = current_progress_line.split("\t")
			geno_file = arr[0]
			Q_file = arr[1]
			R_file = arr[2]
			if(arr[3] == "find_new"):
				flag_refined_pos = 1 # Need to refine
			num_QTLs = int(arr[4])


# Read in the file of previous computations if we have found QTLs before. Otherwise, generate them.
prev_pos = []
prev_genotypes = []
prev_pos = np.array(prev_pos, dtype=np.int32)
prev_genotypes = np.array(prev_genotypes)
q = []
r = []
if(num_QTLs != 0):
	# This is restarting.
	prev_pos = np.fromstring(current_pos_line, dtype=int, sep="	")
	flag_load_prev = 0

	try:
		prev_genotypes = np.load(args.dir + "/" + geno_file)
	except:
		flag_load_prev = 1
		pass

	size_of_prev_genome = (prev_pos.size)

	# Consistent prev_pos and prev_genotypes?
	if(flag_load_prev == 1 or prev_genotypes.shape[1] != size_of_prev_genome):
		# We have to remake it from the prev_pos line.
		prev_genotypes = np.ones((num_usable_spores,size_of_prev_genome))
		for pos_index in range(len(prev_pos)):
			pos = prev_pos[pos_index]
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
	
		# Do we have q,r?
		flag_remake = 0
		if(os.path.isfile(args.dir + "/" + Q_file) and os.path.isfile(args.dir + "/" + R_file)):
			#q = np.load(args.dir + "/" + Q_file) 
			#r = np.load(args.dir + "/" + R_file)

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

	
else:
	size_of_prev_genome = 0


# Ok, we've now reloaded all the previous computations. 
# Set up computation settings
poolcount = args.cpu*2
num_chrom_to_scan = len(genotypes_file)

def find_QTL(num):
	lowest_RSS = np.Infinity
	genome_at_lowest_RSS = []
	pos_index_at_lowest_RSS = 0
	last_q = []

	#start = time.clock()
	for chr in range(num_chrom_to_scan):
		loc = chrom_startpoints[chr_to_scan[chr]]
		for i in range(0 + num, num_lines_genotypes[chr_to_scan[chr]], poolcount):
			if(np.isin(loc+i, prev_pos)):
				continue
			
			genome_line = genotypes_file[chr_to_scan[chr]][i]
			# Remove genomes that have no phenotypes
			# We need to remove genomes that have no phenotypes and genomes that aren't in the train set
			genomes = np.take(genome_line,train_perm)
			genomes = genomes[~np.isnan(train_set)] 
			genomes = np.reshape(genomes,(num_usable_spores,1)) # A N row by 1 column matrix
			
			WX = genomes * np.sqrt(np.reshape(errors,(num_usable_spores,1)))  # X = X * sqrt(W)  -> N by 1
			QtX = np.dot(np.transpose(q),WX) # Gets the scale for each vectors in Q.   # Q^t * X -> k by 1
			QtX_Q = np.einsum('ij,j->i',q,np.ravel(QtX))   # Dot product of Q and Q^t * X, but shaped as a single vector. This is the sum of all the projections of the new genotype on Q
			orthogonalized = WX-np.reshape(QtX_Q,(num_usable_spores,1)) # Orthogonalize: Remove the projections from the real vector.
			new_q = orthogonalized/np.linalg.norm(orthogonalized) # Orthonormalize: Now do final conversion.
			# This gets the last column of Q.
			# We only need the last column of Q to get the new residuals. We'll assemble the full Q or the full R if we need it (i.e. to obtain betas).

			q_upTy = np.einsum('i,i', np.ravel(new_q), phenotypes * np.sqrt(errors))
			q_upq_upTy = np.ravel(new_q) * q_upTy
			predicted_fitnesses = initial_predicted_fitnesses + q_upq_upTy/np.sqrt(errors)

			# Scale the intercept term
			mean_predicted_fitnesses = np.mean(predicted_fitnesses)
	
			# RSS
			RSS = np.sum((phenotypes - mean_phenotypes - predicted_fitnesses + mean_predicted_fitnesses)**2) # This is the RSS for 1:1 line.
			#print(str(loc+i) + "	" + str(RSS))

			if(RSS < lowest_RSS):
				lowest_RSS = RSS
				genome_at_lowest_RSS = genomes.copy()
				pos_index_at_lowest_RSS = loc+i # Position is zero indexed.
				last_q = np.ravel(new_q)
				#last_r = np.ravel(np.dot(np.transpose(q_up),WX))


	# Now return the values
	return lowest_RSS,genome_at_lowest_RSS,pos_index_at_lowest_RSS, last_q


def refine_positions(num):
	lowest_RSS = np.Infinity
	genome_line_lowest_RSS = []
	pos_index_at_lowest_RSS = -1
	last_q = []
	start_index = chrom_startpoints[chr_of_snp-1]
	for scan_pos in range(left_bracket + num, right_bracket + 1, poolcount):
		genome_line = genotypes_file[chr_of_snp-1][scan_pos-start_index]

		# Remove the genomes that are not in the train set
		pos_line = np.take(genome_line,train_perm)
			
		# Remove the genomes that have no phenotypes
		pos_line = pos_line[~np.isnan(train_set)]
		pos_line = np.reshape(pos_line,(num_usable_spores,1)) # A N row by 1 column matrix

		WX = pos_line * np.sqrt(np.reshape(errors,(num_usable_spores,1)))
		QtX = np.dot(np.transpose(q_down),WX) # Gets the scale for each vectors in Q.
		QtX_Q = np.einsum('ij,j->i',q_down,np.ravel(QtX))
		orthogonalized = WX-np.reshape(QtX_Q,(num_usable_spores,1)) # Orthogonalize
		new_q = orthogonalized/np.linalg.norm(orthogonalized) # Orthonormalize
		# This gets the last column of Q.
		# Now let's get the last column of R.
		# Assemble q_up
		q_upTy = np.einsum('i,i', np.ravel(new_q), phenotypes * np.sqrt(errors))
		q_upq_upTy = np.ravel(new_q) * q_upTy
		predicted_fitnesses = initial_predicted_fitnesses + q_upq_upTy/np.sqrt(errors)
			
		# Scale the intercept term
		mean_predicted_fitnesses = np.mean(predicted_fitnesses)
	
		# RSS
		RSS = np.sum((phenotypes - mean_phenotypes - predicted_fitnesses + mean_predicted_fitnesses)**2) # This is the RSS for 1:1 line.
		#print(str(scan_pos) + "	" + str(RSS))
		if(RSS < lowest_RSS):
			lowest_RSS = RSS
			genome_line_lowest_RSS = pos_line.copy()
			pos_index_at_lowest_RSS = scan_pos
			last_q = np.ravel(new_q)


	# Return the values
	return lowest_RSS,genome_line_lowest_RSS,pos_index_at_lowest_RSS, last_q

# Let's code it to run a loop of fixed size instead for now, so it'll be easier to think about it.

while(num_QTLs < args.maxqtl):

	# Put the intercept
	base_genotypes = np.ones((num_usable_spores,1+size_of_prev_genome))
	base_genotypes[:,1:] = prev_genotypes # First index is the intercept.

	if(num_QTLs == 0): # If num_QTLs is not zero, then it was loaded previously.
		q,r = np.linalg.qr(base_genotypes * np.sqrt(np.reshape(errors,(num_usable_spores,1))))

	if(flag_refined_pos == 0):
		start_find_new = time.time()
		# Obtain the initial predicted fitness
		initial_beta = linalg.solve_triangular(r,np.dot(np.transpose(q), phenotypes * np.sqrt(errors)), check_finite=False) # 3.49s for 10000 loops
		# first beta index is the intercept term.
		initial_predicted_fitnesses = np.dot(q,np.dot(r,initial_beta))*1/np.sqrt(errors) # Optimal multiplication order

		# Time to search for a new QTL
		p = Pool(poolcount)
		results = p.map(find_QTL, range(poolcount))
		p.close()
		p.join()

		# Now parse the results
		lowest_RSS = np.Infinity
		genome_at_lowest_RSS = []
		pos_index_at_lowest_RSS = 0
		beta_at_lowest_RSS = []
		last_q = []

		for i in range(len(results)):
			RSS = results[i][0]

			if(RSS < lowest_RSS):
				lowest_RSS = RSS
				genome_at_lowest_RSS = results[i][1]
				pos_index_at_lowest_RSS = results[i][2] # Position is zero indexed
				last_q = results[i][3]	

		
		# Update Q/R to get the the new Beta
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
		likelihood = num_usable_spores * math.log(lowest_RSS/num_usable_spores)
		with open(args.dir + "/" + args.log, "a+") as logfile:
			print(likelihood, file=logfile)
			if(size_of_prev_genome > 0):
				print(*prev_pos, sep="\t", end="", file=logfile)
				print("\t", end="", file=logfile)
			print(pos_index_at_lowest_RSS, file=logfile)
			print(*beta[1:len(beta)], file=logfile)
			print("pickle_geno.npy" + "	" + "pickle_q.npy" + "	" + "pickle_r.npy" + "	" + "find_new" + "	" + str(num_QTLs + 1), file=logfile)
		

		# Update the values (remove intercept)
		base_genotypes[:,:-1] = prev_genotypes
		base_genotypes[:,size_of_prev_genome] = np.matrix.flatten(genome_at_lowest_RSS)
		q = q_up
		r = r_up
		if(args.nosave == 0):
			np.save(args.dir + "/" + "pickle_geno", base_genotypes) # No longer has the intercept term
			np.save(args.dir + "/" + "pickle_q", q)
			np.save(args.dir + "/" + "pickle_r", r)

		num_QTLs = num_QTLs + 1
		flag_refined_pos = 1
		prev_pos = np.append(prev_pos, pos_index_at_lowest_RSS)
		prev_genotypes = base_genotypes
		size_of_prev_genome = (prev_pos.size)
		chr_of_snp = np.searchsorted(np.array(chrom_startpoints),pos_index_at_lowest_RSS + 0.5) # Deal with ties

		print("Found new QTL (" + str(num_QTLs) + ") @ " + str(chr_of_snp) + ". Took : " + str(time.time() - start_find_new) + " seconds. Likelihood: " + str(likelihood), file=sys.stderr)
		

	#elif(flag_refined_pos == 1 and num_QTLs > 1):
	elif(flag_refined_pos == 1 and num_QTLs % args.refine == 0):
		# Must refine the positions
		# We need to sort all the positions and the genotyping array.
		start_refine = time.time()

		# What was the chr of the last index?
		chr_last_qtl = np.searchsorted(np.array(chrom_startpoints),prev_pos[len(prev_pos)-1]+0.5)

		sorted_indexes = np.argsort(prev_pos) 
		# prev_pos is a list of QTL positions. sorted_index gives the indexes that would sort prev_pos. So prev_pos has QTL positions (0 to 45000), such as 2,8,5,6
		# So sorted_indxes returns an array of size num_QTLs, where each value is the order you would have to take the index to sort the array. in the example, it would be 0, 2, 3, 1
		sorted_prev_pos = np.take(prev_pos,sorted_indexes) # sorted_prev_pos is a sorted prev_pos array. # We then sort prev_pos with the sorted index. Returns an array 2,5,6,8
		reverse_sorted_indexes = np.argsort(sorted_indexes) # This reverses the sorting procedure, which returns the 'order' of the original list, or the position where the value has ended up in. It would be: 0, 3, 1, 2

		# We sorted the positions so that we know where to refine the positions of a QTL.
		iterations_max = len(sorted_prev_pos)
		lowest_RSS = np.Infinity
		previous_refined = -1
		for iterations in range(iterations_max):
			sorted_pos_index = np.random.randint(0, len(sorted_prev_pos)) # Numpy is range exclusive (can never return len(sorted_prev_pos) in this case).
			chr_of_snp = np.searchsorted(np.array(chrom_startpoints),sorted_prev_pos[sorted_pos_index] + 0.5) # Deal with ties


			if(sorted_pos_index == previous_refined):
			#if(sorted_pos_index == previous_refined or chr_of_snp != chr_last_qtl): # Heuristic to only refine the last chromosome
				continue

			previous_refined = sorted_pos_index


			# Column downdating of QR.
			q_down,r_down = linalg.qr_delete(q,r, sorted_indexes[sorted_pos_index]+1,1,"col",check_finite=False)
			
			initial_beta = linalg.solve_triangular(r_down,np.dot(np.transpose(q_down), phenotypes * np.sqrt(errors)), check_finite=False) # 3.49s for 10000 loops # Beta for the WEIGHTED phenotypes.
			# first beta index is the intercept term.
			initial_predicted_fitnesses = np.dot(q_down,np.dot(r_down,initial_beta))*1/np.sqrt(errors) # Optimal multiplication order # Obtain the predicted fitnesses in the unweighted world.

			# Now we go through the bracketed positions, and obtain the likelihoods for every position
			# First, let's check the left and right side of the position of interest. If the likelihood is worse, we do not update.
	
			left_bracket = chrom_startpoints[chr_of_snp-1] # Minimally the beginning of the chromosome
			if(sorted_prev_pos[sorted_pos_index]-16 > left_bracket):
				left_bracket = sorted_prev_pos[sorted_pos_index]-16
			if(sorted_pos_index > 0):
				if(sorted_prev_pos[sorted_pos_index-1]+1 > left_bracket):
					left_bracket = sorted_prev_pos[sorted_pos_index-1]+1

			# Now find the right bracket
			right_bracket = chrom_endpoints[chr_of_snp-1]
			if(sorted_prev_pos[sorted_pos_index]+16 < right_bracket):
				right_bracket = sorted_prev_pos[sorted_pos_index]+16
			if(sorted_pos_index < len(sorted_prev_pos)-1):
				if(sorted_prev_pos[sorted_pos_index+1]-1 < right_bracket):
					right_bracket = sorted_prev_pos[sorted_pos_index+1]-1

			left_bracket = int(left_bracket)
			right_bracket = int(right_bracket)

			p = Pool(poolcount)
			results = p.map(refine_positions, range(poolcount))
			p.close()
			p.join()

			# Parse the results
			lowest_RSS = np.Infinity
			genome_line_lowest_RSS = []
			pos_index_at_lowest_RSS = -1
			for i in range(len(results)):
				RSS = results[i][0]
				if(RSS < lowest_RSS):
					lowest_RSS = RSS
					genome_line_lowest_RSS = results[i][1]
					pos_index_at_lowest_RSS = results[i][2]

			# We now have the best results.
			# Did the position change? If not, then we do nothing.
			if(pos_index_at_lowest_RSS != sorted_prev_pos[sorted_pos_index]):
				# Else, we update the position
				sorted_prev_pos[sorted_pos_index] = pos_index_at_lowest_RSS
		
				# We update the genome line
				prev_genotypes[:,sorted_indexes[sorted_pos_index]] = genome_line_lowest_RSS[:,0]

				# Update the QR
				q,r = linalg.qr_insert(q_down,r_down,genome_line_lowest_RSS * np.sqrt(np.reshape(errors,(num_usable_spores,1))),sorted_indexes[sorted_pos_index]+1,'col',check_finite=False) # Update the QR decomposition.

		# Done iterating.
		# Obtain betas

		beta = linalg.solve_triangular(r,np.dot(np.transpose(q), phenotypes * np.sqrt(errors)), check_finite=False)
		prev_pos = sorted_prev_pos[reverse_sorted_indexes]

		# Output to log
		# Update the logs
		likelihood = num_usable_spores * math.log(lowest_RSS/num_usable_spores)

		with open(args.dir + "/" + args.log, "a+") as logfile:
			print(num_usable_spores * math.log(lowest_RSS/num_usable_spores), file=logfile)
			print(*prev_pos, sep="\t", file=logfile)
			print(*beta[1:len(beta)], file=logfile)
			print("pickle_geno.npy" + "	" + "pickle_q.npy" + "	" + "pickle_r.npy" + "	" + "refine" + "	" + str(num_QTLs), file=logfile)

		# Update the values
		if(args.nosave == 0):
			np.save(args.dir + "/" + "pickle_geno", prev_genotypes) # No longer has the intercept term
			np.save(args.dir + "/" + "pickle_q", q)
			np.save(args.dir + "/" + "pickle_r", r)
		flag_refined_pos = 0
		print("Attempted to refine QTL. Took : " + str(time.time() - start_refine) + " seconds. Likelihood: " + str(likelihood), file=sys.stderr)
	else:
		flag_refined_pos = 0
exit()

			
