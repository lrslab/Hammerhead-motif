# the motif enrichment functions for hammerhead output


### get example data and run hammerhead to get the bed file indicating the modifided sites
```bash 
wget https://figshare.com/ndownloader/files/46437190 -O ecoli.fa
wget https://figshare.com/ndownloader/files/46437193 -O test.fastq.gz

hammerhead --ref ecoli.fa --read test.fastq.gz --min_depth 5 --min_depth_strand 3
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