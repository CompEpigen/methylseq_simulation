from simulation import read_simulation, bulk_simulation
import pandas as pd
import os

def test_simulation(f_ref, output_dir, n_region):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1)
	assert len(reads["dmr_label"].unique()) == n_region

def test_k_mers(f_ref, output_dir, n_region, k):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, k=k)
	if k > 1:
		assert reads["dna_seq"].apply(lambda x: all([len(xx) == k for xx in x.split(" ")])).all()
		assert (reads["methyl_seq"].apply(lambda x: len(x)) == reads["dna_seq"].apply(lambda x: len(x.split(" ")))).all()
	else:
		assert reads["dna_seq"].apply(lambda x: len(x)==150).all()
		assert reads["methyl_seq"].apply(lambda x: len(x)==150).all()


def test_simulation_long_read(f_ref, output_dir, n_region, read_len):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, len_read=read_len)
	assert len(reads["dmr_label"].unique()) == n_region
	assert reads["methyl_seq"].apply(lambda x: len(x) == read_len).all()

def test_simulation_from_regions(f_ref, output_dir, f_region):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, f_region=f_region, a=1)
	df_region = pd.read_csv(f_region, sep="\t")
	assert len(reads["dmr_label"].unique()) == df_region.shape[0]

def test_simulation_plot(f_ref, output_dir, n_region):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, save_img=True)
	assert len(reads["dmr_label"].unique()) == n_region
	assert os.path.exists(output_dir)
	assert os.path.exists(os.path.join(output_dir, "region_methyl_level_sampling.png"))
	assert os.path.exists(os.path.join(output_dir, "regions/"))
	assert os.path.exists(os.path.join(output_dir, f"regions/region_{n_region-1}.png"))

def test_simulation_bulk(f_ref, output_dir, n_region, n_bulks):
	reads = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, save_img=True)
	bulk_simulation(reads=reads, n_bulks=n_bulks, output_dir=output_dir, std=0)
	assert os.path.exists(os.path.join(output_dir, f"bulk_{n_bulks-1}.txt"))
	assert os.path.exists(os.path.join(output_dir, "bulk_cell_type_proportions.csv"))

def test_random_seed(f_ref, output_dir, n_region):
	reads_1 = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, save_img=True, seed=42)
	reads_2 = read_simulation(output_dir=output_dir, f_ref=f_ref, n_region=n_region, a=1, save_img=True, seed=42)

	assert reads_1.equals(reads_2)
	


if __name__=="__main__":
	f_ref="../genome/hg19.fa"
	f_region="data/regions.csv"
	output_dir="data/output/"
	n_region=20

	test_simulation(f_ref, output_dir, n_region=n_region)
	for k in range(1,5):
		test_k_mers(f_ref, output_dir, n_region=n_region, k=k)
	test_simulation_long_read(f_ref, output_dir, n_region=n_region, read_len=500)
	test_random_seed(f_ref, output_dir, n_region=n_region)
	test_simulation_from_regions(f_ref, output_dir, f_region=f_region)
	test_simulation_plot(f_ref, output_dir, n_region=n_region)
	test_simulation_bulk(f_ref, output_dir, n_region=n_region, n_bulks=1)
	print("Test all passed!")