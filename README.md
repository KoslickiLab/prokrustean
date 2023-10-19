# ProKrustean

```terminal
cmake -B build .
cd build
make
```

Construct prokrustean
```terminal
./prokrustean -i ../../prokrustean_data/SRR20044276.bwt -l 10 -e -f
```

Construct prokrustean (text)
```terminal
./prokrustean -i ../../prokrustean_data/SRR20044276.bwt -l 10 -c
```

Distinct kmers (requires both bwt, prokrustean)
```terminal
./prokrustean -i ../../prokrustean_data/SRR20044276.bwt -l 10
./prokrustean_kmer -i ../../prokrustean_data/SRR20044276.bwt -g ../../prokrustean_data/SRR20044276.bwt.prokrustean -k 15
```

Maximal unitigs (-e required: character extensions)
```terminal
./prokrustean -i ../../prokrustean_data/SRR20044276.bwt -l 10 -e
./prokrustean_unitig_count -i ../../prokrustean_data/SRR20044276.bwt.prokrustean 
```