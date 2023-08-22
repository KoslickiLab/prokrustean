SEQ_FASTQ="../data/seq/---.fastq.gz"
BWT="../data/bwt/---.fastq.gz"
mkdir -p "../data"
mkdir -p "../data/seq"
mkdir -p "../data/bwt"

# ropebwt2 project
OUTPUT="../data/bwt/sample.ropebwt2.bwt"
git clone git@github.com:lh3/ropebwt2.git
cd ropebwt2
make
./ropebwt2 -o $OUTPUT
# remark: line break exists at the last line.