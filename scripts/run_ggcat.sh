#!/bin/bash

# constant
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



cd "$(dirname "$0")/../../ggcat"
ulimit -n 10000

rm -rf ./log 
mkdir ./log
> ./log/"$logfile"

for ((k=left; k<=right; k++)); do
    rm -rf ./output
    mkdir ./output
    echo "Running ./ggcat -i $file and value $k"

    # ggcat build -k "$k" -j "$THREADS" "$file" -o ./output/out"$k".ggcat
    ( time ggcat build -k "$k" -j "$THREADS" "$file" -o ./output/out"$k".ggcat > /dev/null 2>&1 ) 2>> temp_time.log
    real_time=$(grep real temp_time.log)
    echo "$real_time" >> ./log/"$logfile"
    rm temp_time.log
done
