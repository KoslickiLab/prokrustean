
# Counting k-mers
We get one of the results introduced in [the paper](https://www.biorxiv.org/content/10.1101/2023.11.21.568151). First, download the eBWT of [ERR3450203](https://www.ncbi.nlm.nih.gov/sra/?term=ERR3450203), a metagenomic dataset (short reads).

```
wget -O ERR3450203.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EcSCz4NjH1ZPjRzuhxgrk50BCHziK9JLz3eyLK1LJPL0GA?e=0mIaNV&download=1"
```

Compute Prokrustean graph. 
```
./prokrustean -i ./ERR3450203.bwt
```

Count k-mers of a range of k
```
./prokrustean_kmer_count -p ./ERR3450203.bwt.prokrustean
```

### Printing distinct k-mers
As a proof of concept, we support printing k-mers. We showcase the process with a smaller dataset. First, download the eBWT of [SRR20044276](https://www.ncbi.nlm.nih.gov/sra/?term=SRR20044276).
```
wget -O SRR20044276.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EakiUQImGQhLuMI8n7PAPIIBda3Qje88lVxqcy5-BeVQIA?e=skQrOA&download=1"
```
Recover sequences of the bwt.
```
./prokrustean_access -i ./SRR20044276.bwt
```
Generate the Prokrustean graph.
```
./prokrustean -i ./SRR20044276.bwt
```
Print the distinct k-mers. Below code assumes ./SRR20044276.bwt.prokrustean.access exists. 
```
./prokrustean_kmer -p ./SRR20044276.bwt.prokrustean -k 20
```

# Counting unitigs
Refer to the previous section to get the Prokrustean graph of [ERR3450203](https://www.ncbi.nlm.nih.gov/sra/?term=ERR3450203). 

Then below code counts maximal unitigs of a range of k.
```
./prokrustean_unitig_count -p ./ERR3450203.bwt.prokrustean
```

### Printing the compacted de Bruijn graph
The command ./prokrustean_unitig_count counts the maximal unitigs of the given sequence set. Similar to printing distinct k-mers, we confirm that de Bruijn graphs of any k can be extracted.

Again, below code assumes ./SRR20044276.bwt.prokrustean.access exists. 
```
./prokrustean_cdbg -p ./SRR20044276.bwt.prokrustean -k 20
```

# Computing Bray-Curtis Dissimilarity
We compute the Bray-Curtis dissimilarity between two samples for all k values. We refer to the datasets used in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10696585/). First, download two datasets [SAMEA14382676](https://www.ebi.ac.uk/ena/browser/view/SAMEA14382676) and [SAMEA14382411](https://www.ebi.ac.uk/ena/browser/view/SAMEA14382411), which are both metagenomic sequencing data from an infant of 12 months and 3 weeks, respectively. 
```
# Work with new directory.  
mkdir braycurtis
cd braycurtis

wget -O sample0.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR975/001/ERR9752001/ERR9752001_1.fastq.gz"

wget -O sample1.fastq.gz "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR975/008/ERR9751998/ERR9751998_1.fastq.gz"
```
Merge two sequences.
```
zcat sample0.fastq.gz sample1.fastq.gz | gzip > merged_samples.fastq.gz
```
Use BCR to generate the eBWT. (Note, using BCR is beneficial for this task because the order of the reads in the bwt is preserved, and hence the sample ids can be assigned in the form of $0,...,0,1,...,1$.)
```
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA
cd BCR_LCP_GSA
make FASTQ=1
./BCR_LCP_GSA ../braycurtis/merged_samples.fastq.gz ../braycurtis/merged_samples.bwt
```
Compute the Prokrustean graph.
```
./prokrustean -i ../braycurtis/merged_samples.bwt
```
Generate the sample ids. Count sequences in each file and generate rows annotating each sequence with the sample id 0 or 1. The generated file sample_ids.txt contains 11678075 rows of 0s and 21575210 rows of 1s.
```
num0=$(zcat sample0.fastq.gz 2>/dev/null | grep -c "^@" || echo "0")
num1=$(zcat sample1.fastq.gz 2>/dev/null | grep -c "^@" || echo "0")
awk -v num0="$num0" -v num1="$num1" 'BEGIN { for (i=1; i<=num0; i++) print 0; for (i=1; i<=num1; i++) print 1; }' > sample_ids.txt
``` 
Below code outputs Bray-Curtis dissimilarities between sample 0 and 1 for all k values. 
```
./prokrustean_braycurtis -p ../braycurtis/merged_samples.bwt.prokrustean -s ../braycurtis/samples_ids.txt
```
