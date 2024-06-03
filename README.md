The Procrustean graph is a space-efficient substring index designed for rapid enumeration of k-mer contexts across varying k values, with a time complexity that remains constant regardless of the k range. The Procrustean graph is designed to effectively describe the structural transitions of k-mers with respect to k values. We utilize this capability to compute metrics relevant to k-mer-based applications in pangenomics and metagenomics.

This repository covers (1) Prokrustean graph construction from an ebwt of a sequence set, and (2) applications using Prokrustean graphs. Please cite the paper below, when using our code for academic purposes:
> **Adam Park and David Koslicki.** "Prokrustean Graph: A substring index supporting rapid enumeration across a range of k-mer sizes." [![dx.doi.org/10.1101/2023.11.21.568151](https://img.shields.io/badge/doi-10.1101%2F2023.11.21.568151-blue.svg)](http://dx.doi.org/10.1101/2023.11.21.568151)

# Quick start
### Install
Prokrustean is a header-only stand-alone c++ project. The graph construction utilizes the wavelet tree of the input ebwt in [SDSL project](https://github.com/simongog/sdsl-lite), but we copied the relevant code from the project to avoid dependency issues.
```
git clone https://github.com/KoslickiLab/prokrustean.git
cd prokrustean
cmake -B build .
cd build
make
```

### Make a prokrustean graph
This approach requires a eBWT file that includes a single string with a seperator (normally '$'). Either download [our example data](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/EpeyylRQoyhAmi60bt8ne3IBaTDVXzsdVVB8ODAKZ0CPRw?e=dCw1Oi), or refer to the [BWT Section](#BWT) below to make your own. Below code generates a binary file (./SRR20044276.bwt.prokrustean) representing the prokrustean graph of a short read dataset of SRR20044276.

```
# kmin is 20 by default.
./prokrustean -i ./SRR20044276.bwt -l 20
```
### Application: Count distinct k-mers
The command below counts the number of k-mers for a range k, [$kmin$, the maximum length possible]. Note that the output may not match that of other libraries like KMC. This is due to Prokrustean not using canonical k-mers, restricting to an `A,C,T,G` alphabet only, etc. But to check correctness, we have utilized brute-force implementations for testing purposes. ([/tests/naive_impl.cpp](/tests/naive_impl.cpp)).
```
./prokrustean_kmer_count -p ./SRR20044276.bwt.prokrustean
```

### Application: Count maximal unitigs of de Bruijn graphs
Below command counts the number of maximal unitigs of de Bruijn graphs for a range of k, [$kmin$, the maximum length possible]. 
```
./prokrustean_unitig_count -p ./SRR20044276.bwt.prokrustean
```

### Application: Compute Bray-Curtis
The command below computes the Bray-Curtis dissimilarities between samples across a range of k, [$kmin$, the maximum length possible]. For comparing two read sets, `A.fastq.gz` and `B.fastq.gz`, first construct `merged.bwt.prokrustean` from the ebwt of the merged files. Then, the second parameter annotes each sample ids for each sequence, `merged.bwt.samples_ids.txt`. This file contains N rows of 0s or 1s, where N is the total number of sequences. Each integer in the ith row indicates whether the ith sequence in the Procrustean graph is from A (0) or B (1).
```
./prokrustean_braycurtis -p ./merged.bwt.prokrustean -s ./merged.bwt.samples_ids.txt
```

### Application: Count vertex degrees of the overlap graph
The prokrustean graph contains an hierarchical representation of the overlap graph of the original sequence set. This compact representation makes computing vertex degrees of the overlap graph beneficial. Below command outputs the incoming and outgoing degrees of each sequence, where $kmin$ is the overlap threshold (overlap at least $kmin$ are counted)
```
./prokrustean_overlap -p ./merged.bwt.prokrustean
```

<!-- # Examples -->
<!-- More detailed usages are introduced in the `/examples` folder.  

[/examples/datasets](examples/datasets.md) includes datasets used to generate figures in the paper. 

[/examples/experiments](examples/experiments.md) details how the actual experiments in the paper are conducted.
-->

<div id="BWT"></div>

## How to get eBWTs
There are many libraries computing eBWTs from fastq or fasta files. Their outputs may differ based on the employed strategies, but all are valid for generating the Procrustean graph. Here, we provide some examples that worked in our tasks. 
### optimalBWT
We used [optimalBWT](https://github.com/davidecenzato/optimalBWT) for normal cases on linux environments.
```
git clone https://github.com/davidecenzato/optimalBWT.git
cd optimalBWT
git submodule update --init --recursive
make
```
Run optimalBWT with a sample data, e.g. download this fastq file [SRR20044276.fastq](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR20044276&display=download) and unzip it.
```
python3 optimalBWT.py SRR20044276.fastq SRR20044276.fastq --algorithm sais --fasta --verbose

# Change the file name for simplicity of subsequent processes.
mv SRR20044276.fastq.optbwt SRR20044276.bwt
``` 

### BCR
[BCR](https://github.com/giovannarosone/BCR_LCP_GSA) was specifically used for correctness tests and computing Bray-Curtis dissimilarity because it preserves the order of sequences from the original dataset. This preservation is crucial for annotating sequences; for example, in computing Bray-Curtis dissimilarity, each sequence should be annotated with thier sample Id.
```
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA
cd BCR_LCP_GSA
make FASTQ=1
```
Then generate eBWT.
```
./BCR_LCP_GSA ./SRR20044276.fastq.gz ./SRR20044276.bwt
```

### grlBWT
[grlbwt](https://github.com/ddiazdom/grlBWT) is another recently published library. We used it specifically when other compilers were incompatible with our Mac-based environments. This library requires a concatenated sequence string as input, thus some preprocessing is essential. 

To install the prerequisites for creating bwts, SDSL and grlBWT, please do the following:
```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cd ..
git clone https://github.com/ddiazdom/grlBWT.git
cd grlBWT
mkdir build
cd build
cmake ..
make
```
You will then need to make sure that grlBWT is on your `$PATH`. One way to do this is:
```
mkdir -p ~/bin
cp grl2plain ~/bin/
cp bwt_stats ~/bin/
cp grlbwt2rle ~/bin/
cp grlbwt-cli ~/bin/
cp reverse_bwt ~/bin/
cp split_runs ~/bin/
export PATH="$HOME/bin:$PATH"
cd ../..
```

There is a helper script `scrpts/get_bwt_grl.py` that converts fastq/a files into a bwt. If you wish to use it, you will need to install BioPython:
```
conda install -c bioconda biopython
```

#### Generate ebwt
Below script uses [grl bwt](https://github.com/ddiazdom/grlBWT) to compute the ebwt of a sequence file (fastq.gz). Any ebwt using a seperator (normally '$') is accepted. 
```
# concatenate sequences, run grlbwt, convert to a text file.
python3 scripts/get_bwt_grl.py -i some_sequences.fastq.gz
```

### ropebwt2

[ropebwt2](https://github.com/lh3/ropebwt2/blob/master/main.c) may have been familiarized by some users. However, it is not recommended to use the library to generate the input of this project as it generates an incorrect data for fairly large datasets. For example, get the bwt of SRR20044276.fasq.gz with ropebwt2 and recover the sequences with ./prokrustean_access. It fails at recovering the original sequences. The same task consistently works with both grlBWT and optimalBWT, and our naive implementation of bwt developed for testing purposes.

# Supplementary Details
## Structure of the prokrustean graph
* Vertices represent sequences and maximal repeats, labeled with their sizes.
* Edges represent stratifying regions (defined in [the paper](https://www.biorxiv.org/content/10.1101/2023.11.21.568151)) in a sequence or maximal repeat. For each vertex, a list of position:repeatId pairs are stored. Therefore, intervals of a stratifying region is extracted by (position, position + size(repeatId)).
* For each maximal repeat, "left/right counts" are also collected as default. Those mean the number of characters that can extend left/right each maximal repeat. For example, for a maximal repeat R, if cR exists in a sequence, then the left count increases by 1. 

Use -c to obtain a readable prokrustean graph to see stored information. This result is not reusable in applications.
```
./prokrustean -i ./SRR20044276.bwt -l 20 -c
```

## Accessing texts
Prokrustean graph is not a full text index, so sequences are identified by their order which depends on how the bwt library works. Since ebwt lacks any standard,  [as discussed in this paper](https://arxiv.org/abs/2202.13235), sequence indices in ebwts can either follow the original sequence order or lexicographical order or some other information maximizing the compression of the index. Below code simply recovers sequences from the given ebwt to support printing output - kmers or de Bruijn graphs.
```
# generates ./SRR20044276.bwt.prokrustean.sequences
./prokrustean_access -i ./SRR20044276.bwt
```

## Correctness of Applications
All applications were tested alongside their corresponding brute-force approaches. The results were not directly compared bit by bit with those of other tools mentioned in the paper, such as KMC and GGCAT. This is because these mature tools process the dataset in a specific manner to better support practitioners. For instance, [KMC](https://github.com/refresh-bio/KMC) offers numerous options for filtering k-mers by default, including parameters like abundances, canonical k-mers, and a limited alphabet consisting of A, G, C, and T. 
* k-mer extracting/counting applications: k-mers from prokrustean were compared with those collected directly from the original file by a naive implementation: chop the sequences by k and gather the unique ones.
* de bruijn graph related applications: compacted de bruijn graph built from prokrustean was compared with the graph naively built from the original file.




