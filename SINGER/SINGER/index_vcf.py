import argparse
import sys
from pathlib import Path
import gzip

# from stackoverflow: https://stackoverflow.com/questions/41525690/open-file-depending-on-whether-its-gz-or-not
def openfile(filename, mode='r'):
    if filename.endswith('gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def index_vcf(input_prefix, segment_length):
    if Path(input_prefix + ".vcf.gz").exists():
        input_file = f"{input_prefix}.vcf.gz"
        gz = True
    else:
        input_file = f"{input_prefix}.vcf"
        gz = False
    index_file = f"{input_prefix}.index"
    
    current_segment_start = -1
    byte_offset = 0
    last_variant_segment = -1
    
    # Clear existing index file
    with open(index_file, 'w') as f:
        f.write("")
    
    # Read VCF file line-by-line
    with openfile(input_file, 'r') as f:
        for line in f:
            # Skip header lines
            if not gz and line.startswith("#"):
                byte_offset += len(line.encode('utf-8'))
                continue
            if gz and line.decode().startswith("#"):
                byte_offset += len(line)
                continue
            
            # Extract position
            pos = int(line.decode().split("\t")[1]) if gz else int(line.split("\t")[1])
            
            # Calculate the start coordinate of the segment this variant belongs to
            segment_start = (pos // segment_length) * segment_length
            
            # If this is the first variant in a new segment, record its byte offset
            if segment_start != current_segment_start:
                # Record this segment start and byte offset
                with open(index_file, 'a') as f_idx:
                    f_idx.write(f"{segment_start}\t{byte_offset}\n")
                
                last_variant_segment = segment_start
                current_segment_start = segment_start
            
            byte_offset += len(line) if gz else len(line.encode('utf-8'))

def main():
    parser = argparse.ArgumentParser(description="Index a VCF file by block length.")
    parser.add_argument("vcf_file_prefix", type=str, help="VCF file prefix without .vcf extension")
    parser.add_argument("segment_length", type=int, help="Length of segments to index")
    args = parser.parse_args()

    print("Index vcf files")
    print("------------------")
    print(f"VCF file prefix: {args.vcf_file_prefix}")
    print(f"Block length: {args.segment_length}")
    print("------------------")

    index_vcf(args.vcf_file_prefix, args.segment_length)

if __name__ == "__main__":
    main()

