# ProKrustean

The Prokrustean graph is a space-efficient data structure capable of quickly querying to extract a wide range of information about kmers for any given value of k, or multiple k-values, with a time complexity that is independent of the k size range. A theoretical contribution of the Prokrustean graph is at elucidating the relationship between substring similarities in genomic sequences and the set of kmers. Its expected utility lies in defining quantitative metrics that describe the behaviors of kmer-based applications, which can be rapidly computed using a Prokrustean graph for a range of k values. Additionally, it offers potential benefits in building competitive multi-k (variable-order) representations of sequences or reads, which are particularly promising for genome assembly and alignment.

This repository supports (1) Prokrustean graph construction from an ebwt of a sequence set, and (2) applications using Prokrustean graphs. 

# Quick start
#### Install
Prokrustean is a header-only stand-alone c++ project. The graph construction utilizes the wavelet tree of the input ebwt in [SDSL project](https://github.com/simongog/sdsl-lite), but we copied the relevant code from the project to avoid dependency issues.
```
git clone git@github.com:KoslickiLab/prokrustean.git
cd prokrustean
cmake -B build .
cd build
make
```

#### Make a prokrustean graph
This approach requires a eBWT file that includes a single string with a seperator (normally '$'). Either download [our example data](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/EpeyylRQoyhAmi60bt8ne3IBaTDVXzsdVVB8ODAKZ0CPRw?e=dCw1Oi), or refer to the [BWT Section](#BWT) below to make your own. Below code generates a binary file (./SRR20044276.bwt.prokrustean) representing the prokrustean graph.

```
# -l is the smallest length of maximal repats to be computed. (Lmin in the paper)
./prokrustean -i ./SRR20044276.bwt -l 20
```
#### Application1: Count distinct k-mers
```
# default k=[L..R] (L=Lmin, R=maximum sequence length in the input).
./prokrustean_kmer_count -p ./SRR20044276.bwt.prokrustean
```

#### Application2: Count (maximal) unitigs of de Bruijn graphs
```
# default k=[L..R] (L=Lmin, R=maximum sequence length in the input)
./prokrustean_unitig_count -p ./SRR20044276.bwt.prokrustean
```

#### Make an access to sequences for prokrustean graph
Prokrustean graph is a location-based model that does not directly store sequences. Instead, it refers to the locations of substring similarity information of a given sequence set. Application 3&4 requires the string sequences. Below code simply recover sequences from the given bwt. Since ebwt lacks any standard,  [as discussed in this paper](https://arxiv.org/abs/2202.13235), sequence indices in ebwts may or may not match with the order of original sequences, so the clearest way to secure the correctness is to recover the sequences from ebwt. 
```
# generates ./SRR20044276.bwt.prokrustean.sequences 
./prokrustean_access -i ./SRR20044276.bwt
```  
#### Application3: Print distinct k-mers
```
# ./SRR20044276.bwt.prokrustean.access is assumed to be there.
./prokrustean_kmer -p ./SRR20044276.bwt.prokrustean -k 20
```
#### Application4: Print compacted de Bruijn graph
```
# ./SRR20044276.bwt.prokrustean.access is assumed to be there.
./prokrustean_cdbg -p ./SRR20044276.bwt.prokrustean -k 20
```

<div id="BWT"></div>

## BWT Construction
### optimalBWT
We recommend [optimalBWT](https://github.com/davidecenzato/optimalBWT) for generating an input file. Note, some mac-based compiler fails at installing optimalBWT. Then using the grlBWT at the next section is recommended.
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
### grlBWT
For users who fail to utilize optimalBWT, we recommend [grlbwt](https://github.com/ddiazdom/grlBWT). But this library requires a concatenated sequence string as an input, so the initial sequence files have to be converted.

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

#### Get input (ebwt)
Below script uses [grl bwt](https://github.com/ddiazdom/grlBWT) to compute the ebwt of a sequence file (fastq.gz). Any ebwt using a seperator (normally '$') is accepted. 
```
# concatenate sequences, run grlbwt, convert to a text file.
python3 scripts/get_bwt_grl.py -i some_sequences.fastq.gz
```

Or here are some preprocessed [ebwt files](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/EpeyylRQoyhAmi60bt8ne3IBaTDVXzsdVVB8ODAKZ0CPRw?e=dCw1Oi) 
```
wget -O SRR20044276.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EakiUQImGQhLuMI8n7PAPIIBda3Qje88lVxqcy5-BeVQIA?e=skQrOA"
```
### ropebwt2

[ropebwt2](https://github.com/lh3/ropebwt2/blob/master/main.c) may have been familiarized by some users. However, it is not recommended to use the library to generate the input of this project as it generates an incorrect data for fairly large datasets. For example, get the bwt of SRR20044276.fasq.gz with ropebwt2 and recover the sequences with ./prokrustean_access. It fails at recovering the original sequences. The same task consistently works with both grlBWT and optimalBWT, and our naive implementation of bwt developed for testing purposes.

# Discussions
## Structure of prokrustean
The data structure represents [the Prokrustean Graph](https://www.biorxiv.org/content/10.1101/2023.11.21.568151). 
* Vertices represent sequences and strata. The size of each is stored only. 
* Edges represent stratifying regions in a sequence or stratum. For each vertex, a list of position:stratumId pairs are stored. Therefore, intervals of a stratifying region is extracted by (position, position + size(stratumId)).
* For each stratum, "left/right counts" are also collected as default. Those mean the number of characters that can extend left/right each stratum. For example, for a stratum R, if cR exists in a sequence, then the left count increases by 1. This information is used in de Bruijn graph related applications. 

For checking the structure with a sample dataset, use -c to get a readable prokrustean to see stored information. This result is not reusable in applications.
```
./prokrustean -i ./SRR20044276.bwt -l 20 -c
```
## Correctness of Applications
All applications were tested alongside their corresponding brute-force approaches. The results were not directly compared bit by bit with those of other tools mentioned in the paper, such as KMC and GGCAT. This is because these mature tools process the dataset in a specific manner to better support practitioners. For instance, [KMC](https://github.com/refresh-bio/KMC) offers numerous options for filtering k-mers by default, including parameters like abundances, canonical k-mers, and a limited alphabet consisting of A, G, C, and T. 
* k-mer extracting/counting applications: k-mers from prokrustean were compared with those collected directly from the original file by a naive implementation: chop the sequences by k and gather the unique ones.
* de bruijn graph related applications: compacted de bruijn graph built from prokrustean was compared with the graph naively built from the original file.




