import argparse
import subprocess
import os

def index_vcf(vcf_prefix, block_length):
    """Index the VCF file using the index_vcf.py script."""
    if not os.path.exists(f"{vcf_prefix}.index"):
        print(f"Indexing VCF file: {vcf_prefix}.vcf with block length: {block_length}")
        subprocess.run(["python", "index_vcf.py", vcf_prefix, str(block_length)])
    else:
        print(f"Index file {vcf_prefix}.index already exists. Skipping indexing.")


def run_singer_in_parallel(vcf_prefix, mutation_rate, recombination_to_mutation_ratio, block_length, num_iters, thinning_interval, Ne):
    """Run singer in parallel using the specified parameters."""
    
    # Read the breakpoints from the index file
    with open(f"{vcf_prefix}.index", 'r') as f:
        breakpoints = [int(line.split()[0]) for line in f.readlines()]

    # Prepare the singer commands
    cmd_list = []
    for i in range(len(breakpoints)):
        start = breakpoints[i]
        
        # Base command
        cmd = f"singer -Ne {Ne} -m {mutation_rate} -r {mutation_rate * recombination_to_mutation_ratio} -input {vcf_prefix} -start {start} -end {start + block_length} -n {num_iters} -thinning {thinning_interval}"
        
        cmd_list.append(cmd)

    # Execute the commands in parallel
    subprocess.run(["parallel", "{} :::"] + cmd_list)


def main():
    parser = argparse.ArgumentParser(description="Parallelize singer runs by cutting the chromosome")
    parser.add_argument('-Ne', type=float, required=True, help='Effective population size. Default: 1e4.')
    parser.add_argument("-m", type=float, required=True, help="Mutation rate.")
    parser.add_argument("-ratio", type=float, default=1, help="Recombination to mutation ratio. Default: 1.")
    parser.add_argument("-L", type=int, default=int(1e6), help="Block length. Default: 1e6.")
    parser.add_argument("-vcf", type=str, required=True, help="VCF file prefix (without .vcf or .vcf.gz extension).")
    parser.add_argument("-n", type=int, required=True, help="Number of MCMC samples.")
    parser.add_argument("-thinning", type=int, required=True, help="Thinning interval length.")
    parser.add_argument("-polar", type=float, default=0.5, required=False, help="Site flip probability. Default: 0.5.")
    args = parser.parse_args()

    print("Parameters:")
    print(f"Effective population size: {args.Ne}")
    print(f"Mutation rate: {args.mutation_rate}")
    print(f"Recombination to mutation ratio: {args.recombination_to_mutation_ratio}")
    print(f"Block length: {args.block_length}")
    print(f"VCF file prefix: {args.vcf_file_prefix}")
    print(f"Number of MCMC samples: {args.num_iters}")
    print(f"Thinning interval length: {args.thinning_interval}")
    print(f"Site flip probability: {args.polar}")

    index_vcf(args.vcf_file_prefix, args.block_length)
    run_singer_in_parallel(args.vcf_file_prefix, args.mutation_rate, args.recombination_to_mutation_ratio, args.block_length, args.num_iters, args.thinning_interval)

if __name__ == "__main__":
    main()
