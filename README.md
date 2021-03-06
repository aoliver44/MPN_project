# MPN_project
 Analysis of MPN microbiomes

Data availibility (16S sequences) deposed on the SRA under the bioproject: 

This project was created and maintained using Renv. It is relatively easy to pull down a docker image of R (version 3.4.4 was used for this project) and run most of the statistical analysis done. Included above is a Dockerfile (you need Docker!) if you wish to build the image yourself (R is not light weight...the image is about 2GB).

Prior to statistical analysis, the sequences were processed using qiime2-2018.4 as below:

```
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path raw_data \
   --output-path imported_data.qza

qiime demux emp-paired \
  --m-barcodes-file metadata.tsv \
  --m-barcodes-column BarcodeSequence \
  --i-seqs imported_data.qza \
  --o-per-sample-sequences demux.qza

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 5 \
  --p-trim-left-r 5 \
  --p-trunc-len-f 285 \
  --p-trunc-len-r 238 \
  --o-table table.qza \
  --p-n-threads 32 \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier path/to/515_926classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

unzip table.qza 
# cd into that directory, then the data folder. Should be feature-table.biom in there.
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
biom head -i feature-table.tsv

unzip taxonomy.qza
# cd into that directory, then the data folder. Should be taxonomy.tsv in there.
# Merge by featureID in Bash or R to get final OTU table
```
