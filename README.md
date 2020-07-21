# peptides

Currently the only script in this repo (download_domain.py) is to download sequences from the NCBI refseq that will be compatible with Kraken2. 

This is to get around the issues with using rsync, wget and ftp (that seem to be related to the newer NCBI file structure) that are mentioned on the Kraken2 issues [here](https://github.com/DerrickWood/kraken/issues/114) and [here](https://github.com/DerrickWood/kraken2/issues/272) as none of the fixes listed helped me. 

I have based this upon the scripts at [this repo](https://github.com/fischuu/Kraken_db_install_scripts/) (which crashed before downloading all bacteria for me) to make one script that will download any domain you choose, has the option of downloading only complete genomes or all genomes, can download either protein sequences or nucleotide (DNA) sequences and can download only the human genome from the vertebrate_mammalian section (when given the appropriate options).
It will also not download something that has already been downloaded, and adds the NCBI taxid to each sequence ID so it is compatible with Kraken2. 

It will also give a log file that tells you each of the sequences it has downloaded and added kraken taxids to, but also any files that it had a problem at some point (search for "Didn't" in the text file to see this).

Please feel free to [email me](mailto:robyn.wright@dal.ca) with any questions.

**It requires the additional packages** (as well as the standard os and argparse packages):
- pandas 
```{bash, eval=FALSE}
conda install pandas
```
- biopython
```{bash, eval=FALSE}
conda install -c conda-forge biopython
```

## Example usage

Example 1 - this will download the protein sequences of complete bacterial reference genomes:
```{bash, eval=FALSE}
python download_domain.py --domain bacteria --complete True --ext protein
```


Example 2 - this will download the nucleotide sequences of the reference human genome:
```{bash, eval=FALSE}
python download_domain.py --domain vertebrate_mammalian --complete True --ext dna --human False
```


Example 3 - this will download the nucleotide sequences of all fungal reference genomes (including scaffolds, contigs, etc.):
```{bash, eval=FALSE}
python download_domain.py --domain fungi --complete False --ext dna
```

Just a further note that the human genome is listed as Chromosome rather than Complete. I have made it so that if adding the complete flag will still download the genome, but this is just something to be aware of.
