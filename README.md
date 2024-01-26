# ProKrustean

Prokrustean graph is a space efficient data structure that can quickly be queried to extract essentially arbitrary information regarding kmers for any value of k, or multiple k-values with time complexity essentially independent of k size range. The contribution of Prokrustean graph is to reveal the general structure of similarities of genetic sequences and their connections to kmers. It is expected to be useful at constructing mathematical models that can more precisely anticipate the behaviors of kmer-based applications, and building practical applications that can advance the existing kmer-based methods.

This repository supports (1) Prokrustean graph construction from an ebwt of a sequence set, and (2) applications using Prokrustean graphs. 
Prokrustean graph is a location-based model that does not directly store sequences. Instead, it contains every similarity information of a given sequence set, which is evident in the k-independent applications.

# Quick start
#### Install
Prokrustean is a header-only stand-alone c++ project. The graph construction requires the wavelet tree of [SDSL project](https://github.com/simongog/sdsl-lite), but I copied the required part from the large project to avoid annoying dependency collsions. Credit goes to their contribution.

To install SDSL, please do the following:
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

You are then ready to install Prokrustean:

```
git clone git@github.com:KoslickiLab/prokrustean.git
cd prokrustean
cmake -B build .
cd build
make
```
There is a helper script `get_bwt_grl.py` that converts fastq/a files into a bwt. If you wish to use it, you will need to install BioPython:
```
conda install -c bioconda biopython
```

#### Get input (ebwt)
Below script uses [grl bwt](https://github.com/ddiazdom/grlBWT) to compute the ebwt of a sequence file (fastq.gz). Any ebwt using a seperator (normally '$') is accepted. (Having a simpler way would be great. Any recommendation is appreciated.)
```
# concatenate sequences, run grlbwt, convert to a text file.
python3 scripts/get_bwt_grl.py -i some_sequences.fastq.gz
```

Or here are some preprocessed [ebwt files](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/EpeyylRQoyhAmi60bt8ne3IBaTDVXzsdVVB8ODAKZ0CPRw?e=dCw1Oi) 
```
wget -O SRR20044276.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EakiUQImGQhLuMI8n7PAPIIBda3Qje88lVxqcy5-BeVQIA?e=skQrOA"
```

#### Make a prokrustean graph
```
# -l is the smallest length of maximal repats to be computed. (Lmin in the paper)
./prokrustean -i ./SRR20044276.bwt -l 20
```
Above command generates a binary file. Use -c to get a readable (but not reusable) prokrustean if you want to see information stored.
```
./prokrustean -i ./SRR20044276.bwt -l 20 -c
```
# Applications

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
Application 3&4 requires the string sequences. Below code simply recover sequences from the given bwt and store them so that prokrustean can use. This is such a dumb way making a permutated copy of the sequences. However, since ebwt lacks any standard,  [as discussed in this paper](https://arxiv.org/abs/2202.13235), sequence indices in ebwts may or may not match with the order of sequences in fastq.gz, the most clear way to secure the correctness is to recover the sequences from ebwt.
```
# generates ./SRR20044276.bwt.prokrustean.sequences 
./prokrustean_access -i ./SRR20044276.bwt
```  
Or it can be built during prokrustean construction (-r).
```
./prokrustean -i ./SRR20044276.bwt -l 20 -r
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
# Issues
I will list noteworthy general technical/theoretical issues if there is any.

# Discussions
## Structure of prokrustean
Prokrustean simply represents the Prokrustean Graph in the paper. 
* Sequence or stratum: index and size, no string is stored.
* Regions of a sequence or stratum: a list of (position, stratumId) pairs. Intervals of a region is extracted by (position, position + size(stratumId)).
* For each stratum, "left/right counts" are also collected as default. Those mean the number of characters that can extend left/right each stratum. For example, for a stratum R, if cR exists in a sequence, then the left count increases by 1. This information is used in de Bruijn graph related applications. It can also be directly calculated from the graph, but then the graph has to be scanned one more time, so I decided to use a bit more space (1 byte per stratum) for faster computation.  
* src/prokrustean.hpp contains the code. 

## Credibility of prokrustean and applications
The results of the applications can be very different from other correspondences. For example, the results of k-mer counting of [KMC](https://github.com/refresh-bio/KMC) is different from that of prokrustean. KMC has a lot of options to filter k-mers by default, i.e. abundances, canonical k-mers, and limited alphabet of `A`,`G`, `C`, and `T`. So, the applications were tested by brute-force approaches on a sample of 70 million total length.
* collect k-mers application: k-mers from prokrustean were compared with those collected directly from the original file by a naive implementation: chop the sequences by k and gather the unique ones.
* count k-mers application: k-mer counting application was tested with the result above.
* build de bruijn graphs application:  compacted de bruijn graph built from prokrustean was also compared with the graph built independently from the original file, vertex by vertex and edge by edge.
* count unitigs application:  maximal unitig counting application was tested with the result above.

## Simpler ways to get an ebwt?
ebwt is actively advancing, but not many available softwares are there and some projects have compiler issues. It would be nice to receive any recommendation.
* [ropebwt](https://github.com/lh3/ropebwt2/blob/master/main.c) is currently not used because I experienced it does not correctly recover the sequences, eg. after making SRR20044276.bwt from SRR20044276.fasq.gz file with ropebwt, and recover the sequences with the ./prokrustean_access, the result is different from that of other bwt techniques. There can be many reasons, but anyway it is not compactible with prokrustean for now.
* [grlbwt](https://github.com/ddiazdom/grlBWT) used in this project requires a concatenated string as an input, which is a long journey for generating the bwt. 

