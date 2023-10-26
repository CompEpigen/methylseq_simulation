# Tumour-normal read-level methylation pattern and pseudo-bulk simulator

Read-level methylome simulator using a beta-binomial distribution. 

It currently supports only two cell-type simulations (tumour and normal).

Pseudo-bulk samples with random cell-type compositions can be also simulated with the read-level methylomes. 

## Set-up 
_methylseq_simulation_ requires Python version > 3.6 (So far, it is tested under Python version 3.7 and 3.8).
The dependencies are clarified in `requirements.txt`.

#### pip installation
If your environment (e.g. conda, Python venv and so on) is already activated, you can simply install the dependencies by `pip install -r requirements.txt`

#### Conda environment
[Conda](https://conda.io/projects/conda/en/latest/index.html) is a package and environment manager. 
It helps you with managing software dependencies independently from your local system. 
If you want to start conda, please find their installation guidance [here](https://docs.conda.io/projects/conda/en/23.1.x/user-guide/install/index.html).

1. Create a new conda enviroment with a specific python version `conda creat -n $your_env_name python=$python_version`
2. Install the dependencies by `pip install -r requirements.txt` 

#### Python venv
Python also supports a virtual environment. 

1. Create a Python virtual environment by `python -m venv $your_directory_path`
2. Activate the virtual environment by `source $your_directory_path/bin/activate`
3. Install dependencies by `pip install -r requirements.txt`

#### pip error handling 
If you get a version-related error message as below:
```
Could not find a version that satisfies the requirement numpy==1.21.4
(from versions: 1.14.5, 1.14.6, 1.15.0, 1.15.1, 1.15.2, 1.15.3, 1.15.4,
1.16.0, 1.16.1, 1.16.2, 1.16.3, 1.16.4, 1.16.5, 1.16.6, 1.17.0, 1.17.1,
1.17.2, 1.17.3, 1.17.4, 1.17.5, 1.18.0, 1.18.1, 1.18.2, 1.18.3, 1.18.4,
1.18.5, 1.19.0, 1.19.1, 1.19.2, 1.19.3, 1.19.4, 1.19.5, 1.20.0, 1.20.1,
1.20.2, 1.20.3, 1.21.0, 1.21.1, 1.21.2, 1.21.3)
```
You may want to upgrade your `pip` by `pip install --upgrade pip`

## Quick start
You can simulate read-level methylomes and pseudo-bulk samples (with `--bulk` option) by running `main.py` as below:

````
python main.py --help
usage: main.py [-h] -r F_REF [-d F_REGION] [-o OUTPUT_DIR] [--save_img]
               [-ng N_REGIONS] [-a A] [-b B] [-nr N_READS] [-k K_MER]
               [-l LEN_READ] [--seed SEED] [--bulk] [-nb N_BULKS] [-s STD]

optional arguments:
  -h, --help            show this help message and exit
  -r F_REF, --f_ref F_REF
                        .fasta file for the reference genome
  -d F_REGION, --f_region F_REGION
                        tab-separated .csv file, DMRs should be given with
                        mean methylation level of each cell type, chr, start
                        and end. n_regions will be selected from hg19 CpG
                        islands if the file is not given
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        a directory where all generated results will be saved
                        (default: ./)
  --save_img            Save simulated methylation patterns as a .png file
                        (default: False)
  -ng N_REGIONS, --n_regions N_REGIONS
                        Number of regions to select from CGIs when the region
                        file is not given (default: 100)
  -a A, --a A           alpha parameter of the beta-binomial distribution to
                        simulate read-level methylomes (default: 1e-6)
  -b B, --b B           beta parameter of the beta-binomial distribution to
                        simulate read-level methylomes (default: 5)
  -nr N_READS, --n_reads N_READS
                        Read coverage to simulate in each DMR (default: 120)
  -k K_MER, --k_mer K_MER
                        K to process the simulated read sequences into K-mer
                        sequence (default: 1)
  -l LEN_READ, --len_read LEN_READ
                        Read length to simulate (default: 150)
  --seed SEED           seed number (default: 950410)
  --bulk                Whether you want to generate pseudo-bulks or the
                        entire dataset (default: False)
  -nb N_BULKS, --n_bulks N_BULKS
                        Number of bulks, Applicable only when --bulk is True
                        (default: 1)
  -s STD, --std STD     Standard deviation to sample local proportions. The
                        larger value is given, the more varying local
                        proportions are sampled from a Gaussian distribution
                        centred at the global proportion (default: 0.0)

````
## Output example 
If `--save_img` option is on, you can get the summary of results as figures in the designated `output_dir`. 

#### Sampled region-wise methylation level
<img src="data/output/region_methyl_level_sampling.png" width="500">

#### Simulated read-level methylomes
<img src="data/output/regions/region_1.png" width="500">
