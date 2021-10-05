import string
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import itertools
import time
import scipy.stats as st
sys.path.append('/n/desai_lab/users/klawrence/BBQ/alldata')
from spore_defs import *



def haldane_chrom(d,chrom):
	# distance d is in basepairs
	# 1 cM = 2 kb
	#rates_by_chrom_old = [1384.04365247,2436.76980165,1693.43573485,2276.48993531,2311.91848742,1457.01815933,2471.6244856,2131.20772674,
								#2128.15505231,2376.75200613,2291.77716017,2438.58106786,2539.99772823,2352.42233483,2223.93613249,2297.38772424]
								
	rates_by_chrom = [1578.14349582,
							2508.96380702,
							1877.25014075,
							2312.4896111,
							2326.29062162,
							1939.5698633,
							2503.27060554,
							2173.41070953,
							2225.02620068,
							2389.99027447,
							2315.98587503,
							2470.58724189,
							2502.60184574,
							2467.5457676,
							2551.69201858,
							2412.52317237]
	
	return 0.5*(1-np.exp(-2*d*(0.01/rates_by_chrom[chrom])))

def haldane(d):
	# distance d is in basepairs
	# 1 cM = 2 kb
	rate = 2000.0
	
	return 0.5*(1-np.exp(-2*d*(0.01/rate)))

		
# def inv_haldane(p,d):
# 	return -1*np.log(1-2*p)/(0.01*d/2000)
#
# def inv_haldane_approx(SNP_map,p,chrom,i):
# 	l = SNP_map.get_SNP_dist_between(chrom,i+1,i)
# 	return l/(p*100.0)
	

def get_transition_matrices(SNP_list,chrom):
	length = len(SNP_list[chrom])-1
	T = np.full((length,4,4),np.nan)

	error_prob = 0.01
	return_prob = 0.3

	for i in range(length):
		dist = SNP_list[chrom][i+1]-SNP_list[chrom][i]
		pr  = haldane_chrom(dist,chrom)
		#pr  = haldane(dist)
		# Upper left: recombination probabilities between correct states
		T[i][0][0] = 1-pr-error_prob
		T[i][1][0] = pr
		T[i][0][1] = pr
		T[i][1][1] = 1-pr-error_prob
		# Upper right: prob to enter error states
		T[i][0][2] = error_prob
		T[i][0][3] = 0
		T[i][1][2] = 0
		T[i][1][3] = error_prob
		# Lower left: Prob to return from error states
		T[i][2][0] = return_prob
		T[i][2][1] = 0
		T[i][3][0] = 0
		T[i][3][1] = return_prob
		# Lower right: Prob to stay in error states
		T[i][2][2] = 1-return_prob
		T[i][2][3] = 0
		T[i][3][2] = 0
		T[i][3][3] = 1-return_prob

	return(T)
	
def get_observation_probs(RM_reads,BY_reads,total_avg_reads):
	O = np.full([len(RM_reads),4,4],0.0)

	nRM = np.nan_to_num(RM_reads)
	nBY = np.nan_to_num(BY_reads)
	N = nRM+nBY
	p_error = 0.01
	
	for i in range(len(N)):
		if total_avg_reads > 1 and N[i] > 10*int(total_avg_reads): 
			N[i] = 0
			nRM[i] = 0
			nBY[i] = 0
		
	pRM = st.binom.pmf(nRM,N,(1-p_error))
	pBY = st.binom.pmf(nRM,N,p_error)
		
	O[:,0,0] = pRM[:]
	O[:,1,1] = pBY[:]
	O[:,2,2] = pBY[:]
	O[:,3,3] = pRM[:]		
			
	return(O)

def forward(T,O):
	fprobs_norm = np.full((len(O),4),np.nan,dtype=np.longdouble)
	f_log_weights = np.full(len(fprobs_norm),np.nan,dtype=np.longdouble)
	fprobs_norm[0][0] = 0.49
	fprobs_norm[0][1] = 0.49
	fprobs_norm[0][2] = 0.01
	fprobs_norm[0][3] = 0.01
	f_log_weights[0] = np.log(sum(fprobs_norm[0]))
	for i in range(1,len(fprobs_norm)):
		fprobsi = np.dot(np.dot(fprobs_norm[i-1],T[i-1]),O[i])
		weight = sum(fprobsi)
		if weight == 0:
			print(fprobs_norm[i-1])
			print(T[i-1])
			print(O[i])
			print(i)
			#sys.exit()
		fprobs_norm[i] = fprobsi/weight
		f_log_weights[i] = np.log(weight)+f_log_weights[i-1]
		
	return(fprobs_norm,f_log_weights)
	
def backward(T,O):
	bprobs_norm = np.full((len(O),4),np.nan)
	b_log_weights = np.full(len(bprobs_norm),np.nan)
	bprobs_norm[-1][0] = 1.0
	bprobs_norm[-1][1] = 1.0
	bprobs_norm[-1][2] = 1.0
	bprobs_norm[-1][3] = 1.0
	b_log_weights[-1] = 0
	for i in range(-2,-1*len(bprobs_norm)-1,-1):
		bprobsi = np.dot(np.dot(T[i+1],O[i+1]),bprobs_norm[i+1])
		weight = sum(bprobsi)
		bprobs_norm[i] = bprobsi/weight
		b_log_weights[i] = np.log(weight)+b_log_weights[i+1]
			
	return(bprobs_norm,b_log_weights)
	
def posteriors(fprobs_norm,f_log_weights,bprobs_norm,b_log_weights,T,O):

	total_log_L = f_log_weights[-1]

	# initialize vectors for posteriors
	log_posteriors = np.full((len(fprobs_norm),4),np.nan)
	post_norms = np.full(len(fprobs_norm),np.nan)
	
	# intialize vector for recomb map quantities
	post_RM_norms = np.full(len(fprobs_norm),np.nan)
	post_BY_norms = np.full(len(fprobs_norm),np.nan)
	post_RMerr_norms = np.full(len(fprobs_norm),np.nan)
	post_BYerr_norms = np.full(len(fprobs_norm),np.nan)
	
	trans_RM_BY = np.full(len(T),np.nan)
	trans_BY_RM = np.full(len(T),np.nan)	
			
	for s in range(4):
		log_posteriors[:,s] = np.log(fprobs_norm[:,s])+f_log_weights[:]+np.log(bprobs_norm[:,s])+b_log_weights[:]-total_log_L
				
	post_norms[:] = np.exp(log_posteriors[:,0])+np.exp(log_posteriors[:,2])
	
	post_RM_norms[:] = np.exp(log_posteriors[:,0])
	post_BY_norms[:] = np.exp(log_posteriors[:,1])
	post_RMerr_norms[:] = np.exp(log_posteriors[:,2])
	post_BYerr_norms[:] = np.exp(log_posteriors[:,3])
	
	# add to recombination counts for recomb map
	trans_RM_BY[:] = np.exp(np.log(fprobs_norm[:-1,0])+f_log_weights[:-1]+np.log(T[:,0,1])+np.log(bprobs_norm[1:,1])+b_log_weights[1:]+np.log(O[1:,1,1])-total_log_L) 
	trans_BY_RM[:] = np.exp(np.log(fprobs_norm[:-1,1])+f_log_weights[:-1]+np.log(T[:,1,0])+np.log(bprobs_norm[1:,0])+b_log_weights[1:]+np.log(O[1:,0,0])-total_log_L) 
	
	
	return(log_posteriors,post_norms,post_RM_norms,post_BY_norms,post_RMerr_norms,post_BYerr_norms,trans_RM_BY,trans_BY_RM)
	

#######	

# Program to infer missing SNPs using HMM forward-backward algorithm
# Uses SNP location map and sequenced genome
t0 = time.time()


# Read SNP map
SNP_reader = csv.reader(open('/n/desai_lab/users/klawrence/BBQ/alldata/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')
SNP_list = genome_to_chroms(genome_str_to_int(next(SNP_reader)))
num_chroms = len(SNP_list)
num_SNPs = [len(x) for x in SNP_list]
num_SNPs_total = sum(num_SNPs)
print(num_SNPs,file=sys.stdout,flush=True)
print(num_SNPs_total,file=sys.stdout,flush=True)


# Get transition matrices
transition_matrices = []
for i in range(num_chroms):
	transition_matrices.append(get_transition_matrices(SNP_list,i))
print('Done finding transition matrices.',file=sys.stdout,flush=True)

	
# Memmaps to read from & write to
num_spores = 384*12*23
spore_mmap =  np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_allseq_spores_reads',dtype=[('address','U3',4),('reads_RM','f8',num_SNPs_total+15),('reads_BY','f8',num_SNPs_total+15) ],mode='r',shape=num_spores)
inferred_genome_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_allseq_inferred',dtype=[('address','U3',4),('genotype','f8',num_SNPs_total+15) ],mode='r+',shape=num_spores)

# vectors to store recombination map
trans_count_RM_BY = []
trans_count_BY_RM = []
state_count_RM = []
state_count_BY = []

for chrom in range(num_chroms):
	length = num_SNPs[chrom]-1
	trans_count_RM_BY.append(np.full(length,0.0))
	trans_count_BY_RM.append(np.full(length,0.0))
	state_count_RM.append(np.full(length,0.0))
	state_count_BY.append(np.full(length,0.0))


# Read list of verified barcode-associated wells, to use for recombination map
wells = {}
with open('/n/desai_lab/users/klawrence/BBQ/alldata/cs/count_files/complete_addresses_5thresh.txt','r') as readfile:
	cs_reader = csv.reader(readfile,delimiter='\t')
	for row in cs_reader:
		if len(row) == 6:
			bc1,bc2,batch,set,plate,well = row
			address = tuple([batch,set,plate,well])
			wells[address] = [bc1,bc2]
	readfile.close()
print('Number of verified wells: ',len(wells),file=sys.stdout,flush=True)


batch = int(sys.argv[1])
batch_start = 4608*(batch-1)


# Keep track of the number of good genomes (not blanks, no errors)
num_inferred_genomes = 0
inferred_index = batch_start
num_recomb_map_genomes = 0

# initialize some variables that will get re-calculated for each spore
this_spore_reads_RM = reads_to_chroms(spore_mmap[0]['reads_RM'])
this_spore_reads_BY = reads_to_chroms(spore_mmap[0]['reads_BY'])
all_reads_RM = np.array([item for chrom in range(num_chroms) for item in this_spore_reads_RM[chrom]])
all_reads_BY = np.array([item for chrom in range(num_chroms) for item in this_spore_reads_BY[chrom]])
total_avg_reads = (np.nansum(all_reads_RM)+np.nansum(all_reads_BY))/float(len(all_reads_RM))
all_reads = np.nansum([all_reads_RM,all_reads_BY],axis=0)
coverage = float(np.count_nonzero(all_reads))/float(len(all_reads))

# loop over spores
for i in range(batch_start,batch_start+4608):
	
	if i%100 == 0: print(str(i),str(inferred_index),file=sys.stdout,flush=True)
	
	# Check that this spore has reads
	this_spore_reads_RM = reads_to_chroms(spore_mmap[i]['reads_RM'])
	this_spore_reads_BY = reads_to_chroms(spore_mmap[i]['reads_BY'])
	if len(this_spore_reads_RM) != 16: 
		continue
		
	all_reads_RM = np.array([item for chrom in range(num_chroms) for item in this_spore_reads_RM[chrom]])
	all_reads_BY = np.array([item for chrom in range(num_chroms) for item in this_spore_reads_BY[chrom]])
	total_avg_reads = (np.nansum(all_reads_RM)+np.nansum(all_reads_BY))/float(len(all_reads_RM))
		
	all_reads = np.nansum([all_reads_RM,all_reads_BY],axis=0)
	coverage = float(np.count_nonzero(all_reads))/float(len(all_reads))
	if coverage > 0.26 and tuple(spore_mmap[i]['address']) in wells.keys():
		num_recomb_map_genomes += 1
		
	# Create array to store the newly inferred genome	
	inferred_genome = []
				
	# Loop over chromosomes
	for chrom in range(num_chroms):

		# Pick out reads for this chrom
		all_reads_RM = np.array(this_spore_reads_RM[chrom])
		all_reads_BY = np.array(this_spore_reads_BY[chrom])

		# Pick out transition matrices for this chrom
		T = transition_matrices[chrom]
				
		# Get observation matrices
		O = get_observation_probs(this_spore_reads_RM[chrom],this_spore_reads_BY[chrom],total_avg_reads)
		if np.isnan(O).any():
			print('Error! nan in observation probabilites, chrom '+str(chrom)+'. ',file=sys.stderr,flush=True)	
			print(i,inferred_index)
			inferred_genome = []		
			break	
				
		# Get forward & backward probabilities
		fprobs_norm,f_log_weights = forward(T,O) 	
		bprobs_norm,b_log_weights = backward(T,O)
		if np.isnan(fprobs_norm).any():
			print('Error! nan in fprobs, chrom '+str(chrom)+'. ',file=sys.stderr,flush=True)	
			print(i,inferred_index)		
			inferred_genome = []			
			break	
		if np.isnan(bprobs_norm).any():
			print('Error! nan in bprobs, chrom '+str(chrom)+'. ',file=sys.stderr,flush=True)
			print(i,inferred_index)	
			inferred_genome = []				
			break	

		# Get total log likelihood
		log_likelihood = f_log_weights[-1]

		# Get posterior probs 
		log_posteriors,post_probs,post_probs_RM,post_probs_BY,post_probs_RMerr,post_probs_BYerr,trans_chrom_RM_BY,trans_chrom_BY_RM = posteriors(fprobs_norm,f_log_weights,bprobs_norm,b_log_weights,T,O)
		
		if np.isnan(post_probs).any():
			print('Error! nan in inferred genome, chrom '+str(chrom)+'. ',file=sys.stderr,flush=True)	
			print(i,inferred_index)		
			inferred_genome = []			
			break	
		
		# If no errors, add inferred genome to list
		inferred_genome.append(list(post_probs))
		
		# If coverage was above median and spore barcode is verified, add this spore's info to recomb map
		if coverage > 0.26 and tuple(spore_mmap[i]['address']) in wells.keys():
			trans_count_RM_BY[chrom] += trans_chrom_RM_BY
			trans_count_BY_RM[chrom] += trans_chrom_BY_RM
			state_count_RM[chrom] += post_probs_RM[:-1]
			state_count_BY[chrom] += post_probs_BY[:-1]
		
		
	# If there were no errors, should have 16 inferred chrom genomes; otherwise go to next spore
	if len(inferred_genome) == 16:
		
		inferred_genome_mmap[inferred_index]['genotype'] = chroms_to_genome(inferred_genome)
		inferred_genome_mmap[inferred_index]['address'] = spore_mmap[i]['address']	
		num_inferred_genomes += 1
		inferred_index += 1


print(str(num_inferred_genomes),file=sys.stdout,flush=True)		
print(str(num_recomb_map_genomes),file=sys.stdout,flush=True)		
		

inferred_genome_mmap.flush()

# write recomb map to files
for chrom in range(num_chroms):
	
	with open('/n/desai_lab/users/klawrence/BBQ/alldata/geno/recomb_mapping/iter2_recomb_map_chrom'+str(chrom+1)+'_batch_'+str(batch)+'.txt','w') as writefile:
		recomb_writer = csv.writer(writefile,delimiter='\t')
			
		recomb_writer.writerow(['Spores used: ',num_recomb_map_genomes])	
		recomb_writer.writerow(['SNP1 loc','SNP2 loc','RM->BY transitions','BY->RM transitions','RM occupancy','BY occupancy'])
		
		for snp in range(len(trans_count_RM_BY[chrom])):
			recomb_writer.writerow([SNP_list[chrom][snp],SNP_list[chrom][snp+1],trans_count_RM_BY[chrom][snp],trans_count_BY_RM[chrom][snp],state_count_RM[chrom][snp],state_count_BY[chrom][snp]])
		
		writefile.close()

