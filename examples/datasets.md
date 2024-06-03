## Metagenomic sequencing data

Below datasets were mainly used to generate figures in the paper. 

[SRR20044276](https://www.ncbi.nlm.nih.gov/sra/?term=SRR20044276) (bacteria short reads < 100mb)

[ERR3450203](https://www.ncbi.nlm.nih.gov/sra/?term=ERR3450203) (human gut metagenomic short reads, 8.3G bases)

[SRR18495451](https://www.ncbi.nlm.nih.gov/sra/?term=SRR18495451) (Long reads, metagenomic sequencing of three ripening washed-rind cheeses, 11.3G bases)

Or here are preprocessed [ebwt files](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/EpeyylRQoyhAmi60bt8ne3IBaTDVXzsdVVB8ODAKZ0CPRw?e=dCw1Oi) that can be downloaded as below.
```
wget -O SRR20044276.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EakiUQImGQhLuMI8n7PAPIIBda3Qje88lVxqcy5-BeVQIA?e=skQrOA&download=1"

wget -O ERR3450203.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EcSCz4NjH1ZPjRzuhxgrk50BCHziK9JLz3eyLK1LJPL0GA?e=0mIaNV&download=1"

wget -O SRR18495451.bwt "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EdxDjR__9RVCkbOu4zWsMPcBu1859XAmWMrc7wPkFqY_dg?e=49340S&download=1"
```

## Pangenomic references

We tested the growth of the Prokrustean graph accumulating chromosome 1 contigs/sequences included in human references as listed below. Note that it is extremely important to remove letters `N` in the references not to contain unnecessarily many maximal repeats of `NNNNN`, which can generate an error during the prokrustean graph construction.

* [AP023461.1](https://www.ncbi.nlm.nih.gov/nuccore/AP023461.1/)
* [CH003448.1](https://www.ncbi.nlm.nih.gov/nuccore/CH003448.1/)
* [CH003496.1](https://www.ncbi.nlm.nih.gov/nuccore/CH003496.1/)
* [CM000462.1](https://www.ncbi.nlm.nih.gov/nuccore/CM000462.1/)
* [CM001609.2](https://www.ncbi.nlm.nih.gov/nuccore/CM001609.2/)
* [CM003683.2](https://www.ncbi.nlm.nih.gov/nuccore/CM003683.2/)
* [CM009447.1](https://www.ncbi.nlm.nih.gov/nuccore/CM009447.1/)
* [CM009872.1](https://www.ncbi.nlm.nih.gov/nuccore/CM009872.1/)
* [CM010808.1](https://www.ncbi.nlm.nih.gov/nuccore/CM010808.1/)
* [CM021568.2](https://www.ncbi.nlm.nih.gov/nuccore/CM021568.2/)
* [CM034951.1](https://www.ncbi.nlm.nih.gov/nuccore/CM034951.1/)
* [CM035659.1](https://www.ncbi.nlm.nih.gov/nuccore/CM035659.1/)
* [CM039011.1](https://www.ncbi.nlm.nih.gov/nuccore/CM039011.1/)
* [CM045155.1](https://www.ncbi.nlm.nih.gov/nuccore/CM045155.1/)
* [CM073952.1](https://www.ncbi.nlm.nih.gov/nuccore/CM073952.1/)
* [CM074009.1](https://www.ncbi.nlm.nih.gov/nuccore/CM074009.1/)
* [CP139523.1](https://www.ncbi.nlm.nih.gov/nuccore/CP139523.1/)
* [NC_000001.11](https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11/)
* [NC_060925.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_060925.1/)

Below code downloads the eBWT of the merged sequences of all 19 references above.
```
wget -O pangenome-merged19.bwt.zip "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EX9_6KkOD9RGumH7gRpIcekBH8vPYu2-fOn_Xd_9c9K3Gg?e=xe5RL3&download=1"

unzip pangenome-merged.bwt.zip
```

## Metagenomic references
We used e.coli references to check the growth of the Prokrustean graph. 

[3682 E. coli assemblies in NCBI circa 2020](https://zenodo.org/records/6577997)

Below code downloads the eBWT of the merged sequences of 3500 references.
```
wget -O metagenome-merged3500.bwt.zip "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/Edc0BvCVE29Clx7k4wa_HjEBAnQ5NI7ne3X48wpIu7nT1g?e=20kbcf&download=1"

unzip metagenome-merged3500.bwt.zip
```
