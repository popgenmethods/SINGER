import sys

def index_vcf(input_prefix, segment_length):
    input_file = f"{input_prefix}.vcf"
    index_file = f"{input_prefix}.index"
    
    current_segment_start = -1
    byte_offset = 0
    last_variant_segment = -1
    
    # Clear existing index file
    with open(index_file, 'w') as f:
        f.write("")
    
    # Read VCF file line-by-line
    with open(input_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                byte_offset += len(line.encode('utf-8'))
                continue
            
            # Extract position
            pos = int(line.split("\t")[1])
            
            # Calculate the start coordinate of the segment this variant belongs to
            segment_start = (pos // segment_length) * segment_length
            
            # If this is the first variant in a new segment, record its byte offset
            if segment_start != current_segment_start:
                # Record this segment start and byte offset
                with open(index_file, 'a') as f_idx:
                    f_idx.write(f"{segment_start}\t{byte_offset}\n")
                
                last_variant_segment = segment_start
                current_segment_start = segment_start
            
            byte_offset += len(line.encode('utf-8'))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <vcf_file_prefix> <segment_length>")
        sys.exit(1)
    
    input_prefix = sys.argv[1]
    segment_length = int(sys.argv[2])
    
    index_vcf(input_prefix, segment_length)
