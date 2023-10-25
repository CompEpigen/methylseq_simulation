import pandas as pd
from scipy.stats import beta

import matplotlib.pyplot as plt
import seaborn as sns
import os

import logging

logger = logging.getLogger(__name__)

def _sample_from_beta(n_regions: int, a: float=1e-6, b: float=5.0):
	'''
	Sample methylation level using the beta distribution
	
	n_regions: int
		Number of regions to simulate mean methylation level
	a, b: float (default: 1e-6 and 5.0, each)
		Alpha and beta values for the beta distribution
	'''
	# Create a beta distribution
	rv = beta(a, b)
	meanMethy1 = rv.rvs(size=n_regions)

	return meanMethy1, rv

def sample_region_meth(regions: pd.DataFrame, a: float=1e-6, b: float=5.0, 
					   plot: bool=False, save_dir: str="./") -> pd.DataFrame:
	'''
	Sample methylation level for the regions using the beta distribution
	
	regions: pd.DataFrame
		pandas DataFrame containing the regions
	a, b: float (default: 1e-6, 5.0)
		alpha and beta parameters for the beta distribution
	plot: bool (default: False)
		Whether save a histogram of sampled methylation levels 
	save_dir: str (default: ./)
		Directory to save the created histogram. Only applicable when plot=True

	return: pd.DataFrame
		regions in pandas DataFrame in which meanMethy1, meanMethy2 and ctype=T (stands for Tumour) are included 
	'''

	# Sample mean methyl value
	meanMethy1, rv = _sample_from_beta(n_regions=regions.shape[0], a=a, b=b)
	regions["meanMethy1"] = meanMethy1
	regions["meanMethy2"] = 1 - meanMethy1
	regions["ctype"] = "T"

	if plot:
		# Plot the simulated results
		fig, ax = plt.subplots(1, figsize=(5, 3))
		sns.histplot(regions["meanMethy1"], bins=20, color="#e056fd", alpha=0.3, kde=True, line_kws={"linewidth":3}, label="cell type 1",  ax=ax)
		sns.histplot(regions["meanMethy2"], bins=20, color="#1abc9c", alpha=0.3, kde=True, line_kws={"linewidth":3}, label="cell type 2",  ax=ax)

		ax.set_title("Sampled with beta(a=%.2f, b=%.2f)"%(a, b))
		ax.set_ylim(0, regions["meanMethy1"].value_counts(bins=20).max()+5)
		ax.set_xlabel("Mean methyl level in regions")
		ax.set_xlim(0,1)
		ax.legend(loc='best', frameon=False)
		plt.tight_layout()
		fig.savefig(os.path.join(save_dir, "region_methyl_level_sampling.png"), format="png", dpi=300)
		plt.cla()
		plt.close(fig)
	else:
		if save_dir is not None:
			logger.warning(f"Although plot option is off, {save_dir} is given as a path to save the plot. It will be ignored.")

	return regions
