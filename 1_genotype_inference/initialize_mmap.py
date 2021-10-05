import string
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import itertools
import time
sys.path.append('/n/desai_lab/users/klawrence/BBQ/alldata')
from spore_defs import *

#######	

# Read SNP map
SNP_reader = csv.reader(open('/n/desai_lab/users/klawrence/BBQ/alldata/BYxRM_nanopore_SNPs.txt','r'),delimiter='\t')
SNP_list = genome_to_chroms(genome_str_to_int(next(SNP_reader)))
num_chroms = len(SNP_list)
num_SNPs = [len(x) for x in SNP_list]
num_SNPs_total = sum(num_SNPs)
print(num_SNPs,file=sys.stdout,flush=True)
print(num_SNPs_total)

# Initialize memmaps for reads & inferred genotypes
num_spores = 384*12*23
# novareads_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_novaseq_spores_reads',dtype=[('address','U3',4),('reads_RM','f8',num_SNPs_total+15),('reads_BY','f8',num_SNPs_total+15) ],mode='w+',shape=num_spores)
# oldreads_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_oldseq_spores_reads',dtype=[('address','U3',4),('reads_RM','f8',num_SNPs_total+15),('reads_BY','f8',num_SNPs_total+15) ],mode='w+',shape=num_spores)
# allreads_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_allseq_spores_reads',dtype=[('address','U3',4),('reads_RM','f8',num_SNPs_total+15),('reads_BY','f8',num_SNPs_total+15) ],mode='w+',shape=num_spores)

inferred_genome_mmap = np.memmap(filename='/n/scratchlfs02/desai_lab/klawrence/BBQ/alldata/geno/memmaps/BYxRM_allseq_inferred',dtype=[('address','U3',4),('genotype','f8',num_SNPs_total+15) ],mode='w+',shape=num_spores)

