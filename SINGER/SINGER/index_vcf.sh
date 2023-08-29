#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <vcf_file> <segment_length>"
  exit 1
fi

# Input VCF file and segment length
input_file="$1"
segment_length="$2"
index_file="${input_file}.byte_index"

# Initialize variables
current_segment_start=-1
byte_offset=0
last_recorded_segment_start=-1

# Clear existing index file
> "$index_file"

# Skip header lines and read VCF file line-by-line
grep -v '^#' "$input_file" | while read -r line; do
  # Extract position
  pos=$(echo "$line" | cut -f2)

  # Calculate the start coordinate of the segment this variant belongs to
  segment_start=$(( (pos / segment_length) * segment_length ))

  # If this is the first variant in a new segment, record its byte offset
  if [[ $segment_start -ne $current_segment_start ]]; then
    # If the previous segment had no variants, mark it
    if [[ $segment_start -gt $((last_recorded_segment_start + segment_length)) ]]; then
      echo -e "$((last_recorded_segment_start + segment_length))\t-1" >> "$index_file"
    fi

    echo -e "$segment_start\t$byte_offset" >> "$index_file"
    last_recorded_segment_start="$current_segment_start"
    current_segment_start="$segment_start"
  fi

  # Update byte offset (add the length of the line + 1 for the newline)
  byte_offset=$((byte_offset + ${#line} + 1))
done

# If the last segment had no variants, mark it
if [[ $current_segment_start -eq $last_recorded_segment_start ]]; then
  echo -e "$current_segment_start\t-1" >> "$index_file"
fi
