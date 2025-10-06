# Creating a conda enviroment for fastp

conda create -n fastp -c bioconda fastp
conda activate fastp

# Install the required tools within the enviroment
conda install -c bioconda fastp -y

# to deactivate the current conda enviroment use
conda deactivate

# Download the .fastq files from SRA explorer 
# LINK: https://sra-explorer.info/?utm_source=chatgpt.com
# In our case the project number (PRJNA1051785)
# execute the bash downloading script on frontend to
mkdir data
cd data

mkdir raw
cd raw

cat sra_explorer_sra_download.sh  | parallel -j 8 
# downloads 8 fastq files in parallel

# Download the reference genome from NCBI
mkdir genomes
cd genomes

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz

# Download the annotation file from NCBI
mkdir annotation
cd annotation

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz

# Now we have all the data to start the RNA-seq analysis
# Now we start with the QC using fastqc and multiqc to aggregate the results

# we already have the directory for our raw fastq files in data/raw
# Now we start directly with the trimming using fastp (it also does QC and trimming in the same step
# The trimmed directory will be created in the data/raw directory
# so we have data/raw/trimmed
# GNU parallel is used to parallelize the trimming of multiple samples. Here we are using 8 cores (-j 8) to trim 8 samples simultaneously.

mkdir -p trimmed

cat samples.txt | parallel -j 8 '
fastp \
  -i {}_sample_1.fastq.gz \
  -I {}_sample_2.fastq.gz \
  -o trimmed/{}_R1.trimmed.fastq.gz \
  -O trimmed/{}_R2.trimmed.fastq.gz \
  -w 8 \
  -h trimmed/{}_fastp.html \
  -j trimmed/{}_fastp.json
'


# After trimming we will do the hisat2 alignment
# First we need to build the hisat2 index from the reference genome
# This step is done only once
# The index files will be created in the genomes directory

module load gcc/14.2.0
module load hisat2

cd ../genomes

hisat2-build /path/to/genomes/ref_genome.fna.gz ref_genome_index

# This will create multiple index files with the prefix ref_genome_index
# Now we can align the trimmed reads to the reference genome using hisat2
# The output will be in the form of .sam files
# We will create a directory for the sam files
mkdir -p ../data/raw/trimmed/sam
cd ../data/raw/trimmed

# we will create a samples.txt files to store all the base names of the samples
# This file will be used to loop through all the samples for alignment

ls *_R1.trimmed.fastq.gz | sed 's/_R1.trimmed.fastq.gz//' > samples.txt

# Now for the hisat2 alignment again we will use GNU parallel to align multiple samples simultaneously

cat samples.txt | parallel -j 8 ' hisat2 -p 8 \ -x /path/to/genomes/ref_genome.fna.gz ref_genome_index \ -1 {}_R1.trimmed.fastq.gz \ -2 {}_R2.trimmed.fastq.gz \ -S aligned/{}.sam '

# Note: upto trimming step we are working with paired end reads, meaning 12 samples will have 24 fastq files (12 R1 and 12 R2).
# After alignment we will have 12 .sam files in the alignments/sam directory

# After alignment we will convert the .sam files to .bam files and sort them using samtools
# We will create a directory for the bam files which will need to be sorted and indexed
mkdir -p ../alignments/bam

module load samtools

mkdir -p bam_sorted

ls aligned/*.sam | sed 's/.sam//' | parallel -j 8 '
samtools sort -@ 8 -o bam_sorted/{/}.sorted.bam {}.sam
'
 # followed by indexing the sorted bam files

 ls bam_sorted/*.sorted.bam | parallel -j 8 'samtools index {}'
# This will create .bai index files for each sorted bam file

# Flagstat summary and mapping summary



# Now we will use featureCounts to count the reads mapped to each gene
# The output will be a matrix of counts with genes as rows and samples as columns
# We will create a directory for the counts files
mkdir -p ../counts

module load subread

featureCounts -T 8 -p -t exon -g gene_id \
  -a /path/to/annotation.gtf \
  -o counts/gene_counts.txt \
  bam_sorted/*.sorted.bam
# The output file gene_counts.txt will have the counts for each gene in each sample

# Now we can use the counts file for downstream analysis like differential expression analysis using DESeq2 or edgeR in R
# This part will be done in R and is not included in this script
# End of the RNA-seq pipeline