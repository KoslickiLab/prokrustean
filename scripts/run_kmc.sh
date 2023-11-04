#!/bin/bash

# constant
ABUNDANCE=1
THREADS=8

#  variables
file=""
left=""
right=""
logfile="timings.log"

# Parse command-line options
while getopts "i:l:r:" opt; do
  case $opt in
    i) file="$OPTARG"
    ;;
    l) left="$OPTARG"
    ;;
    r) right="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
  esac
done

# Check if the necessary parameters are provided
if [ -z "$file" ] || [ -z "$left" ] || [ -z "$right" ]; then
  echo "Usage: $0 -i <file> -l <left> -r <right>"
  exit 1
fi

# ../KMC/bin/kmc -k24 -ci1 -t8 ../data/downloads/SRR20044276.part_002.fastq.gz 24mers ../KMC/output
cd "$(dirname "$0")/../../KMC/bin"
ulimit -n 2048

rm -rf ../log 
mkdir ../log
> ../log/"$logfile"

for ((k=left; k<=right; k++)); do
    rm -rf ../output
    mkdir ../output
    echo "Running ./kmc -i $file and value $k"
  
    # ./kmc  -k"$k" -ci"$ABUNDANCE" -t"$THREADS" -hp "$file" ../output/"$k"mers ../output
    ( time ./kmc  -k"$k" -ci"$ABUNDANCE" -t"$THREADS" -hp "$file" ../output/"$k"mers ../output > /dev/null 2>&1 ) 2>> temp_time.log
    real_time=$(grep real temp_time.log)
    echo "$real_time" >> ../log/"$logfile"
    rm temp_time.log 
done

# Execute the kmc file
# ./kmc -k24 -ci1 -t8 ../../data/downloads/SRR20044276.part_002.fastq.gz 24mers ../output
