# the motif enrichment functions for hammerhead output


### get example data and run hammerhead to get the bed file indicating the modifided sites
```bash 
wget https://figshare.com/ndownloader/files/46437190 -O ecoli.fa

# small test data, only have 642 sites
# The data will give 642 modified sequences from BED, which only contains the GATC sites, not the CCWGG sites.

wget https://figshare.com/ndownloader/files/46437193 -O test.fastq.gz

# full test data, have 30000 sites for GATC and 1600 sites for CCWGG
wget https://figshare.com/ndownloader/files/42651082 -O ecoli_re.fastq.gz

# run the full set
#hammerhead --ref ecoli.fa --read test.fastq.gz --min_depth 5 --min_depth_strand 3
hammerhead --ref ecoli.fa --read  ecoli_re.fastq.gz --min_depth 5 --min_depth_strand 3


```
The code will generate the files blow:
```bash
ls
ecoli.fa      enrichment.bed  mapping.mpileup.txt  potential_modification_site.bed  test.fastq.gz
ecoli.fa.fai  enrichment.fa   mapping.sort.bam     potential_modification_site.txt
```


### start the motif enrichment analysis using two files:
-  potential_modification_site.bed 
-  ecoli.fa
```bash