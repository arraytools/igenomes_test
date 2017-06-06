
# UCSC/hg19

## Create dbsnp vcf (& tbi index) files from chromosome 1 only

### http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz # 75MB
### http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz # 70MB, 

cd ~/GSE48215/hg19
gunzip chr1.fa.gz
/opt/SeqTools/bin/bwa-0.7.12/bwa index chr1.fa

awk '$1=="1"'  '/home/brb/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf' > genes_chr1.gtf

## Create .dict file
java -Xmx10g -jar /opt/SeqTools/bin/picard-tools-1.141/picard.jar CreateSequenceDictionary R=chr1.fa O=chr1.dict
## Create .fai index file
/opt/SeqTools/bin/samtools-1.3/samtools faidx chr1.fa

## Create dbsnp vcf (& tbi index) files from chromosome 1 only

### Since the chromosome in chr1.gtf file uses 'chr', we download the dbsnp file from the gatk version.

cd ~/GSE48215/hg19
cp ~/SeqTestdata/usefulvcf/hg19/gatk/common_all_20160601.vcf.gz* .
# guznip common_all_20160601.vcf.gz # Don't do this!

# subset
export PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3:$PATH  # tabix and bgzip
tabix -h common_all_20160601.vcf.gz chr1: > common_all_20160601_1.vcf
# zip

bgzip -c common_all_20160601_1.vcf > common_all_20160601_1.vcf.gz
# create index file
tabix -f -p vcf common_all_20160601_1.vcf.gz

rm common_all_20160601.vcf*
