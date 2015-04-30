

Dependencies
```
pip install bx-python
```

Doing a transformation

```
wget https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF4.0/geneSet.v4_0.gaf
./gaf2gtf.py geneSet.v4_0.gaf > geneSet.v4_0.gtf
synapse get syn3243594
./interval2matrix.py -gtf geneSet.v4_0.gtf --fix-chrome --log2-ave -bed broad.mit.edu_TGCT_Genome_Wide_SNP_6.hg19.cna.bed > broad.mit.edu_TGCT_Genome_Wide_SNP_6.hg19.cna.matrix
```
