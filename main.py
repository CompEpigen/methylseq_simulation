from simulation import read_simulation, bulk_simulation

import pandas as pd
import os, argparse, logging


def arg_parser():
    parser = argparse.ArgumentParser()
    # Required hyperparameters 
    parser.add_argument("-d", "--f_region", required=False, default=None, type=str, help="tab-separated .csv file, DMRs should be given with mean methylation level of each cell type, chr, start and end")
    parser.add_argument("-o", "--output_dir", required=False, default="./", type=str, help="a directory where all generated results will be saved")
    parser.add_argument("-r", "--f_ref", required=True, type=str, help=".fasta file for the reference genome")
    
    parser.add_argument("-a", "--a", required=False, type=float, default=1e-6, help="alpha parameter of the beta-binomial distribution to simulate read-level methylomes")
    parser.add_argument("-b", "--b", required=False, type=float, default=5, help="beta parameter of the beta-binomial distribution to simulate read-level methylomes")

    parser.add_argument("--save_img", required=False, action="store_true", default=False, help="Save simulated methylation patterns as a .png file")
    parser.add_argument("-ng", "--n_regions", type=int, default=100, help="Number of regions to select from CGIs when the region file is not given")

    #Hyperparameters for read-level methylome simulation 
    parser.add_argument("-nr", "--n_reads", type=int, default=120, help="Read coverage to simulate in each DMR")
    parser.add_argument("-l", "--len_read", type=int, default=150, help="Read length to simulate")
    parser.add_argument("--seed", type=int, default=950410, help="seed number")

    # Hyperparameters for pseudo-bulks 
    parser.add_argument("--bulk", default=False, action="store_true", help="Whether you want to generate pseudo-bulks or the entire dataset")
    parser.add_argument("-nb", "--n_bulks", type=int, default=None, help="Number of bulks, Applicable only when --bulk is True")
    parser.add_argument("-s", "--std", type=float, default=0.0, help="Standard deviation to sample local proportions. The larger value is given, the more varying local proportions are sampled from a Gaussian distribution centred at the global proportion")
    
    return parser.parse_args()


if __name__=="__main__":
    args = arg_parser()
    

    # Save hyperparameters - TODO: make it config
    pd.DataFrame({"param":vars(args).keys(),
                  "values":vars(args).values()}).to_csv(args.output_dir+"/parameters.txt", sep="\t")

    reads = read_simulation(output_dir=args.output_dir, 
               f_ref=args.f_ref, 
               f_region=args.f_region, 
               n_region=args.n_regions,
               a=args.a,
               b=args.b, 
               n_reads=args.n_reads,
               len_read=args.len_read,
               save_img=args.save_img,
               seed=args.seed)

    if args.bulk:
        # pseudo_bulk simulation
        bulk_simulation(reads=reads, n_bulks=args.n_bulks, output_dir=args.output_dir, std=args.std)