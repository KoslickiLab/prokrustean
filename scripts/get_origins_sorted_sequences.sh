#!/bin/bash

# Default settings
directory="./"
file_pattern="*"
file_pattern_full=""
output="output.txt"
output_origin="output.txt.origin.txt"
fasta=false
fastq=false

print_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -p <pattern>    Specify the file pattern to search for (default: '*')"
    echo "  -d <directory>  Specify the directory for input and output files (default: './')"
    echo "  -o <output>     Specify the output file name (default: 'output.txt')"
    echo "  -a              Process files as FASTA format"
    echo "  -q              Process files as FASTQ format"
    echo "  -h              Display this help and exit"
    echo ""
    echo "This script processes input FASTQ/FASTA files matching a specified file pattern. It performs the following operations: 
    - Merging: Combines sequences from multiple files into a single dataset.
    - Sorting: Sorts the combined sequences lexicographically based on their content.
    - Indexing: Outputs a list of origin indices corresponding to each sequence in the sorted output. 
    The origin index indicates the file from which each sequence originated, and is listed in a separate informative output file."
    echo "Example:"
    echo "  $0 -p 'S*' -q -o output.txt"
}


# Parse command-line options
while getopts "p:d:o:fq" opt; do
  case $opt in
    p) file_pattern="$OPTARG"  # Pattern for file selection
    ;;
    d) directory="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    f) fasta=true
       fastq=false
    ;;
    q) fastq=true
       fasta=false
    ;;
    h) print_help
    exit 0
    ;;
    \?) print_help
        exit 1
    ;;
  esac
done

# Set file pattern based on fasta or fastq selection
if [ "$fasta" = true ]; then
  file_pattern_full="${directory}${file_pattern}.fasta.gz"
elif [ "$fastq" = true ]; then
  file_pattern_full="${directory}${file_pattern}.fastq.gz"
fi
output="$directory$output"
output_origin="$output.origin.txt"

# Check if either fasta.gz or fastq.gz is selected
if [ "$fasta" = false ] && [ "$fastq" = false ]; then
  echo "Please specify file type: -f for FASTA or -q for FASTQ"
  exit 1
fi

# input parameters
echo "Directory: $directory"
echo "File pattern: $file_pattern"
if [ "$fasta" = true ]; then
  echo "File type: FASTA"
else
  echo "File type: FASTQ"
fi
echo "Input: $file_pattern_full"
echo "Output: $output"
echo "Output origin information: $output_origin"

i=1
> $output_origin  # Clear the file before writing
echo "# Pairs of the origin index and the file" >> "$output_origin"
for f in $file_pattern_full; do
    echo "$i $f" >> "$output_origin"
    ((i++))
done

# Conditional processing based on file type
if [ "$fasta" = true ]; then
    echo "Processing for files $file_pattern_full..."
    i=1
    for f in $file_pattern_full; do
        # echo "Processing "${f} as index ${i}"
        gunzip -c "$f" | awk -v idx="$i" '/^>/ {header=$0; getline seq; print seq "\t" header "\t:" idx;}' &
        i=$((i + 1))
    done | sort -k1,1 | awk -F"\t" '{print $NF}' > "$output"
elif [ "$fastq" = true ]; then
    echo "Processing for files $file_pattern_full..."
    # 
    i=1
    for f in $file_pattern_full; do
        gunzip -c "$f" | awk -v idx="$i" 'NR % 4 == 1 {header=$0} NR % 4 == 2 {seq=$0} NR % 4 == 3 {plus=$0} NR % 4 == 0 {print seq "\t" header "\t" plus "\t" $0 "\t" idx;}' &
        i=$((i + 1))
    done | sort -k1,1 | awk -F"\t" '{print $NF}' > "$output"
fi
