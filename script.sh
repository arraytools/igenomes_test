export PATH=/opt/SeqTools/bin/samtools-1.3:$PATH # samtools
export PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3:$PATH  # tabix and bgzip
export PATH=/opt/SeqTools/bin/bwa-0.7.15:$PATH  # bwa
export PATH=/opt/SeqTools/bin/bowtie2-2.2.9/:$PATH # bowtie2
export PICARD_PATH=/opt/SeqTools/bin/picard-tools-1.141

################
# Ensembl 37
################
mkdir -p ~/igenomes_test/GRCh37
cd ~/igenomes_test/GRCh37
URL_DBSNP=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/common_all_20160601.vcf.gz
CHR=chr1
## Assume the full igenomes has been downloaded to ~/igenomes/Homo_sapiens directory
## bwa index
cp ~/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes/1.fa .
bwa index 1.fa # 3-4 minutes
## gtf file
awk '$1=="1"'  ~/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf > genes_1.gtf
## dict file
java -Xmx10g -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R=1.fa O=1.dict
## Create .fai index file
samtools faidx 1.fa
## Create dbsnp vcf (& tbi index) files from chromosome 1 only
curl -L -O $URL_DBSNP # 1GB, slow, same as UCSC hg38
curl -L -O $URL_DBSNP.tbi
## subset
tabix -h common_all_20160601.vcf.gz chr1: > common_all_20160601_1.vcf
## zip
bgzip -c common_all_20160601_1.vcf > common_all_20160601_1.vcf.gz
## create index file
tabix -f -p vcf common_all_20160601_1.vcf.gz
## delete the original vcf.gz file
## rm common_all_20160601.vcf*

################
# NCBI GRCh38
################
mkdir -p ~/igenomes_test/GRCh38
cd ~/igenomes_test/GRCh38
URL_DBSNP=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/VCF/GATK/common_all_20160527.vcf.gz
CHR=chr1
## Assume the full igenomes has been downloaded to ~/igenomes/Homo_sapiens directory
## bwa index
cp ~/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/chr1.fa .
bwa index chr1.fa # 3-4 minutes
## gtf file
awk '$1=="chr1"'  ~/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf > genes_chr1.gtf
## dict file
java -Xmx10g -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R=chr1.fa O=chr1.dict
## Create .fai index file
samtools faidx chr1.fa
## Create dbsnp vcf (& tbi index) files from chromosome 1 only
curl -L -O $URL_DBSNP # 1GB, slow, same as UCSC hg38
curl -L -O $URL_DBSNP.tbi
## subset
tabix -h common_all_20160527.vcf.gz chr1: > common_all_20160527_1.vcf
## zip
bgzip -c common_all_20160527_1.vcf > common_all_20160527_1.vcf.gz
## create index file
tabix -f -p vcf common_all_20160527_1.vcf.gz
## delete the original vcf.gz file
## rm common_all_20160527.vcf*

#############
# UCSC hg38
#############
mkdir -p ~/igenomes_test/hg38
cd ~/igenomes_test/hg38
URL_DBSNP=ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh38p2/VCF/GATK/common_all_20160527.vcf.gz
CHR=chr1
## Assume the full igenomes has been downloaded to ~/igenomes/Homo_sapiens directory
## bwa index
cp ~/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes/chr1.fa .
bwa index chr1.fa # 3-4 minutes
## gtf file
awk '$1=="chr1"'  ~/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf > genes_chr1.gtf
## dict file
java -Xmx10g -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R=chr1.fa O=chr1.dict
## Create .fai index file
samtools faidx chr1.fa

## bowtie2 index
bowtie2-build chr1.fa genome

## Create dbsnp vcf (& tbi index) files from chromosome 1 only
curl -L -O $URL_DBSNP # 1GB, slow
curl -L -O $URL_DBSNP.tbi
## subset
tabix -h common_all_20160527.vcf.gz chr1: > common_all_20160527_1.vcf
## zip
bgzip -c common_all_20160527_1.vcf > common_all_20160527_1.vcf.gz
## create index file
tabix -f -p vcf common_all_20160527_1.vcf.gz
## delete the original vcf.gz file
## rm common_all_20160527.vcf*

#############
# UCSC/hg19
#############

cd ~/GSE48215/hg19
### http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz # 70MB,
gunzip chr1.fa.gz
/opt/SeqTools/bin/bwa-0.7.12/bwa index chr1.fa

awk '$1=="1"'  ~/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf > genes_chr1.gtf

## Create .dict file
java -Xmx10g -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R=chr1.fa O=chr1.dict
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
