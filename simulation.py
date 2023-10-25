from region import sample_region_meth
from utils import set_seed, set_logger

import numpy as np
from scipy.stats import binom

from Bio import SeqIO

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
import pickle as pk
import os, argparse, random, shutil
import logging

set_logger()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
	
def _long_read_cpg_sample(r_start: int, r_seq: str, dmr):
	cpg_loci = list()
	before_dmr = 0
	after_dmr = 0
	for j in range(len(r_seq) - 1):
		if r_seq[j:j+2] == "CG":
			cpg_loci.append(j)
			if r_start + j < dmr["start"]: 
				before_dmr += 1
			elif r_start + j > dmr["end"]:
				after_dmr +=1
	return cpg_loci, before_dmr, after_dmr


def bulk_simulation(reads: pd.DataFrame, n_bulks: int, output_dir : str, std: float=0.0) -> None:
	'''
	Simulate tumour-normal pseudo-bulks from the simulated reads
	
	reads: pd.DataFrame
		pandas DataFrame containing simulated reads
	n_bulks: int
		Number of bulks to simulate
	output_dir: str (default: ./)
		Directory to save the results
	std: float (default: 0.0)
		Standard deviation of a Gaussian distribution to sample the number of reads in each region
	'''

	assert (n_bulks > 0) and (std > 0), f"n_bulks and std must be > 0 (given n_bulks={n_bulks}, std={std})"

	logger.info(f"Simulate {n_bulks} pseudo-bulk samples from the simulated reads")
	ctypes = reads["ctype"].unique()
	n_dmrs = len(reads["dmr_label"].unique())
	
	# Global cell-type proportions for n pseudo-bulk samples
	pri_ratios = np.random.uniform(size=n_bulks)
	bulk_names = ["bulk_%d"%(i+1) for i in range(n_bulks)]

	post_ratios = list()

	logger.info(f"Creating pseudo-bulk samples. You can find each bulk sample in {output_dir}/bulk_name.txt")
	for bulk_idx, bulk_name in tqdm(enumerate(bulk_names)):
		# cell-type proportions for each region
		local_proportions = np.random.normal(loc=pri_ratios[bulk_idx], scale=std, size=n_dmrs)
		bulk_reads = list()
		
		# Sample reads for each region
		for dmr_idx in reads["dmr_label"].unique():
			dmr_reads = reads.loc[reads["dmr_label"]==dmr_idx,].copy().reset_index().drop(columns=["index"])
			n_reads_per_region = int(dmr_reads.shape[0]/len(ctypes))
			local_prop = local_proportions[dmr_idx]

			for ctype in ctypes:
				ctype_ratio = local_prop if ctype=="T" else 1-local_prop
				sampled_idces = random.sample(dmr_reads.loc[dmr_reads["ctype"]==ctype].index.tolist(), 
											  int(n_reads_per_region * ctype_ratio))
				bulk_reads.append(dmr_reads.loc[sampled_idces, ])

		bulk_reads = pd.concat(bulk_reads)
		bulk_reads.to_csv(f"{output_dir}/{bulk_name}.txt", sep="\t", index=False)

		bulk_ctype_ratio = bulk_reads["ctype"].value_counts(normalize=True)
		bulk_ctype_ratio.name=bulk_name
		post_ratios.append(bulk_ctype_ratio)

	gt_ratios = pd.DataFrame(post_ratios)
	gt_ratios.to_csv(os.path.join(output_dir, "bulk_cell_type_proportions.csv"), sep="\t", index=True)

def read_simulation(f_ref: str, 
			   output_dir: str="./", 
			   f_region: str=None, 
			   n_region: int=1,
			   a: float=1e-6,
			   b: float=5.0, 
			   k: int=1, 
			   n_reads: int=120,
			   len_read: int=150,
			   save_img: bool=False,
			   seed: int=950410) -> pd.DataFrame:
	'''
	Simulate reguion-wise differentially methylated patterns and read-level methylomes 

	f_ref: str
		.fasta or .fa file containing reference genome (e.g. hg19.fa)
	output_dir: str (default: ./)
		Directory to save the results
	f_region: str
		a tab-separated .txt/.csv file containing region information. The file must have a header including "chr", "start", and "end" indicating chromosome, start position and end position of the region
	n_regions: int (default: 1)
		Number of regions to simulate read-level methylomes. Only applicable when f_region=None
	a, b: float (default: 1e-6, 5.0)
		alpha and beta parameters for the beta distribution
	k : int (default: 1)
		K to process the simulated read sequences into K-mer sequence
	n_reads: int (default: 120)
		Number of reads to simulate in each region
	len_read: int (default: 150)
		Read length (bp) for simulation
	save_img: bool (default: False)
		Whether save figures summarising the simulation or not
	seed: int (default: 950410)
		random seed

	return: pd.DataFrame
		pandas DataFrame containing simulated read-level methylomes
	'''

	# Value range handling
	assert (a > 0) and (b > 0) and (k > 0) and (n_reads >0) and (len_read>0), \
			f"a, b, k, n_reads, len_read must be > 0 (given a={a}, b={b}, k={k}, n_reads={n_reads}, len_read={len_read})"

	# Setup random seed 
	set_seed(seed)

	# Setup output files
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	if save_img:
		region_dir = os.path.join(output_dir, "regions/")
		if not os.path.exists(region_dir):
			os.mkdir(region_dir)

	# Reference genome
	record_iter = SeqIO.parse(f_ref, "fasta")

	dict_ref = dict()
	for r in record_iter:
		dict_ref[r.id] = str(r.seq.upper())
	del record_iter

	# Load regions
	if f_region is None:
		logger.info(f"Regions are not given. Top {n_region} with the largest number of CpGs will be chosen from hg19.")
		df_regions = pd.read_csv("data/CGIs.csv", sep="\t")
		df_regions = df_regions[df_regions["chrom"].str.fullmatch("chr\d+")].sort_values(by="perCpg", ascending=False)
		df_regions = df_regions.iloc[:n_region,]
		df_regions = df_regions.rename(columns={"chrom":"chr", "chromStart":"start", "chromEnd":"end"})
	else:
		df_regions = pd.read_csv(f_region, sep="\t")

	# Sample mean methylation level of the regions
	df_regions = sample_region_meth(regions=df_regions, a=a, b=b, plot=save_img, save_dir= output_dir)

	logger.info(f"Mean methyl values are sampled for the regions:\n{df_regions.head()}")
	if ("chr" not in df_regions.columns) or \
		("start" not in df_regions.columns) or \
		("end" not in df_regions.columns):
		logger.error("Regions must be a dataframe with chr, start and end columns")

	reads = list()
	for i in tqdm(range(df_regions.shape[0])):

		dmr = df_regions.iloc[i,:]
		dmr_len = dmr["end"] - dmr["start"]
		

		# Whether long-read seq or not
		if len_read > dmr_len:
			read_sample_start = dmr["start"] - len_read + dmr_len
			read_sample_end = dmr["end"] + len_read - dmr_len

			read_sample_start = max(0, read_sample_start)
			read_sample_end = min(read_sample_end, len(dict_ref[dmr["chr"]]))

			# Methylation patterns should be sample from the same distribution in non-dmr regions
			non_dmr_mean_methy = 0.95 if dmr["meanMethy2"] > 0.5 else 0.05

		else:
			read_sample_start = dmr["start"]
			read_sample_end = dmr["end"]

		# Sample start positions
		starts = np.random.randint(read_sample_start, read_sample_end-len_read, size=n_reads).tolist()

		if save_img:
			methyl_array = np.ones([n_reads, read_sample_end - read_sample_start]) * -1
	
		for r_idx, start in enumerate(starts):
			# simulate same number of reads for both cell types 
			r_ctype = "N" if r_idx > int(n_reads/2) else "T"

			# read position
			end = start + len_read + 1 # for finding "CG" context, this will make (len_read + 1) length of sequence
			r_seq = dict_ref[dmr["chr"]][start:end] # start considering 0-base 
			
			# Set probability based on mean methylation level for each cell type
			p=dmr["meanMethy1"] if r_ctype=="T" else dmr["meanMethy2"]
			
			# Collect CpGs
			if len_read > dmr_len:
				cpg_loci, before_dmr, after_dmr = _long_read_cpg_sample(start, r_seq, dmr)

				# Simulate methylation patterns
				m_patterns_1 = binom.rvs(n=1, p=non_dmr_mean_methy, size=before_dmr)
				m_patterns_2 = binom.rvs(n=1, p=non_dmr_mean_methy, size=after_dmr)
				methyl_patterns = binom.rvs(n=1, p=p, 
											size=len(cpg_loci)-before_dmr-after_dmr).tolist()        
				methyl_patterns = np.concatenate([m_patterns_1, methyl_patterns, m_patterns_2])
			else:
				cpg_loci = [j for j in range(len(r_seq)-2) if r_seq[j:j+2] == "CG"]
				
				# Simulate methylation patterns
				methyl_patterns = binom.rvs(n=1, p=p, size=len(cpg_loci)).tolist()        
				methyl_patterns = np.array(methyl_patterns)


			r_methyl = list(np.ones(len(r_seq), dtype=int)*2)
			for ri, rm in zip(cpg_loci, methyl_patterns):
				r_methyl[ri] = rm

			# Remove the last one included for finding "CG" context
			r_methyl = r_methyl[:-1]
			r_seq = r_seq[:-1]

			assert (len(r_methyl) == len_read) and (len(r_seq) == len_read), f"length of methyl_seq and dna_seq are {len(r_methyl)} and {len(r_seq)}, respectively, although the given read length is {len_read}"			
			if save_img:
				methyl_array[r_idx, start-read_sample_start:(start+len_read-read_sample_start)] = r_methyl[1:-1]
			
			# K-mer
			if k % 2 == 0:
				logger.warning(f"With an even number of K (K={k}), Cytosine methylation in CpG-context cannot be at the centre of each K-mer substring. CpG methylation pattern will given for the K-mer substring where CG is at the middle.")
			
			if k > 1:
				seq_start, seq_end =  k//2 - int(k%2 == 0), k//2
				r_methyl = "".join([str(r) for r in r_methyl[seq_start:-seq_end]])
				r_seq = [r_seq[i:i+k] for i in range(len(r_seq)-(k-1))]
				r_seq = " ".join(r_seq)
			else:
				r_methyl = "".join([str(r) for r in r_methyl])

			reads.append(pd.DataFrame({"dna_seq": [r_seq], 
								 "methyl_seq": [r_methyl],
								 "dmr_ctype":[dmr["ctype"]],
								 "ctype": [r_ctype],
								 "dmr_label":[i],
								 "ref_name":[dmr["chr"]],
								 "ref_pos":[start]}))

		if save_img:

			fig, ax=plt.subplots(1, figsize=(30,5))
			s=sns.heatmap(methyl_array, cmap=["white", "yellow", "black", "grey"], ax=ax)
			ax.set_title("Region %d"%(i))
			ax.hlines([int(n_reads/2)+1], colors="blue", xmin=0, xmax=methyl_array.shape[1])
			plt.tight_layout()
			s.get_figure().savefig(os.path.join(region_dir,f"region_{i}.png"))    
	
	reads = pd.concat(reads)
	reads.to_csv(os.path.join(output_dir, "simulated_reads.txt"), sep="\t", index=False)

	return reads