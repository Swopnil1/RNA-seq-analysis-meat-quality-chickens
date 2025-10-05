# RNA-sequencing pipeline and differential gene expression anaylsis

Swopnil Pradhan

Octotber 1, 2025

## Introduction
RNA sequencing is a high throughput technique used to study the transcriptome (the complete set of RNA molecules within a cell at a given time). 

RNA-seq analysis applications includes:
  1. Gene expression profiling
    Identifying genes which are upregulated or downregulated in response to particular treatment or condition.
  2. Disease research
    Identfying biomarkers or molecular mechanisms in cancer, neurological disorders and infections. Understanding the effects of treating organisms with therapeutics. 
  3.   Alternative splicing analysis
     Understanding different isoforms produced from the same gene.

The analysis explained below is mainly for identifying differentially expressed genes and discovering underlying pathways and mechanisms which can affect meat quality of two different muscle types breast (B) and Leg (L) in two different strains of chicken namely Chinese Dagu chicken (DG) and AA+ broiler roosters (AA). 

## Obtaining raw data from SRA

The dataset which we will be working with comes from (Zhu et al., 2024). To find the raw sequencing data, we can navigate through the SRA explorer https://sra-explorer.info/.  The sequencing data can be accessed through entering the Bioproject accession number **PRJNA1051785** in the SRA explorer. 

<img width="1440" height="813" alt="image" src="https://github.com/user-attachments/assets/fe04fb89-99e9-4124-b43f-5c705b03ef7d" />

Steps to follow after entering PRJNA105178. 
Click "Title" checkbox to select all the samples -> Add to collection -> 12 saved datasets -> SRA downloads -> Download "Bash script for downloading SRA files (nice filenames)" -> Run the 
script sra_explorer_sra_download.sh 

> (Optional)
You can also utilize GNU parallel to simultaneously download multiple files at once.
```
cat sra_explorer_sra_download.sh  | parallel -j 8 
```
The raw data is in FASTQ format which is required for bulk RNA-sequencing analysis. 
The FASTQ file format stores nucleotide sequences with corresponding quality scores primarily used for raw sequencing reads from NGS platforms where quality of each base matters. 

## Download reference genome and annotation files
Now, we have the raw sequencing data in FASTQ format. However, we require additional reference genome file and annotation file as well before we can start our analysis. 

To download the reference genome, we can either do it through NCBI or ENSEMBL. 

```
mkdir genomes
cd genomes

# When downloading from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
```

```
# When downloading from ENSEMBL
wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
gunzip Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
````

Similarly, download the annotation file (in our case in .GTF file format) from

```
mkdir annotation
cd annotation

# When downloading from NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz
gunzip GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gtf.gz
```

```
# When downloading from ENSEMBL
wget https://ftp.ensembl.org/pub/release-115/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf.gz
gunzip Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf.gz
```

Now, we have all the files which we require to start the RNA-seq analysis process. 

## RNA-sequencing pipeline followed in the anaylsis:
1. Trimming and QC using fastp
2. Alignment using Hisat2
3. Sorting and indexing using samtools
4. Read summarization using featureCounts from Subread
5. Differntial expression analysis in RStudio


### Trimming and QC using fastp
fastp is a fast, all in one FASTQ preprocessor which is primarily designed for quality control (QC), filtering and trimming of next-generation sequencing (NGS) reads. Few other major functions of fastp includes adapter trimming, length filtering and quality control reporting by generating HTML and JSON reports. 

> create conda enviroment to install the required dependencies which we use
```
conda create -n fastp -c bioconda fastp
conda activate fastp

conda install -c bioconda fastp -y
```

> We create a file samples.txt which contains only the sample names of each samples and not include each prefix.
```
ls *_sample_1.fastq.gz | sed 's/_sample_1.fastq.gz//' > samples.txt
cat samples.txt

```

```
# Making use of GNU parallel to trim multiple samples simultaneously.
mkdir -p trimmed
cd trimmed

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

conda deactivate fastp
```

<img width="634" height="802" alt="image" src="https://github.com/user-attachments/assets/b078715d-47a0-417a-9856-2ff4dcb4c552" />


### Hisat2 alignment
HISAT2 (Hierarchial Indexing for Spliced Alignment of Transcripts 2) is a fast, memory efficient tool used for aligning RNA-seq reads to a reference genome. Hisat2 is capable of handling the complexities of eukaryotic transcriptomes such as splicing events, alternative isoforms and intron exon boundries meaning Hisat2 is  a splice site aware aligner which is beneficial as compared to standard alingers. 

> Create conda enviroment for Hisat2
```
conda create -n hisat2
conda activate hisat2

conda install -c bioconda hisat2 -y

# to check if hisat2 is installed properly
hisat2 --version
```

> Hisat2 requires building a reference index
```
# genome_index is the index prefix name
hisat2-build /path/to/reference/genome.fa genome_index
```

```
mkdir -p ../data/raw/trimmed/sam
cd ../data/raw/trimmed

# again we can use the samples.txt containing the samples name without the prefix.

cat samples.txt | parallel -j 8 ' hisat2 -p 8 \ -x /path/to/reference/genome.fa genome_index \ -1 {}_R1.trimmed.fastq.gz \ -2 {}_R2.trimmed.fastq.gz \ -S aligned/{}.sam '
```
> During the trimming step we still have paired ends reads for example sample1_R1 and sample1_R2. However, after the alignment step we would have only one read for each sample, meaning in our analysis we go from 24 samples after the trimming step to 12 samples after the alignment.

### Sorting and indexing SAM files 
For the next step we require to convert SAM -> BAM to make the downstream analysis process faster and compatible with the modern bioinfomratics tools. 

Sorting by using genomic coordinates is crucial for downstream process as it orders the alignments based on their positions within the reference genome. 

> Create samtools enviroment using conda
```
conda create -n samtools -c bioconda samtools
conda activate samtools

conda install -c bioconda samtools -y
```

```
mkdir -p ../alignments/bam
ls aligned/*.sam | sed 's/.sam//' | parallel -j 8 '
samtools sort -@ 8 -o bam{/}.sorted.bam {}.sam
'
 # followed by indexing the sorted bam files

 ls bam_sorted/*.sorted.bam | parallel -j 8 'samtools index {}'
# This will create .bai index files for each sorted bam file
```

> Since BAM files are binary files simply using 'less' or 'head' to view the files would not work. We would need to use samtools to view the BAM files. 
```
samtools view -h sample_name.bam | less -S
```
<img width="292" height="239" alt="image" src="https://github.com/user-attachments/assets/a3f8b398-4e01-48fa-abc5-322708bc44ce" />

### Flagstat summary 
In the next step, we will do a flagstat summary using samtools which focuses more on the detailed technical QC snapshot of each sample. It reports counts and percentages for the various alignment catergories - mapped properly, paired, duplicates and singeltons. 

> It is a crucial step, as it assess the quality and integrity of each of the individual sample's alignemnt.

```
# Run flagstat
mkdir -p flagstat_reports
ls bam_sorted/*.sorted.bam | parallel -j 8 '
samtools flagstat {} > flagstat_reports/{/.}.flagstat.txt
'
```
<img width="612" height="244" alt="image" src="https://github.com/user-attachments/assets/31d61802-45d3-40f9-ab1d-53e844f4a604" />

### Mapping summary
Mapping summary is simply a summarized table which contains all the key statistics from multiple flagsta files in the analysis. It helps compare alignment quality scores across all samples at a glance and identify failed runs. 

```
echo -e "Sample\tTotal\tMapped\tPercentMapped" > mapping_summary.tsv

for f in flagstat_reports/*.flagstat.txt; do
    sample=$(basename "$f" .sorted.flagstat.txt)
    total=$(grep "in total" "$f" | awk '{print $1}')
    mapped=$(grep " mapped (" "$f" | head -n1 | awk '{print $1}')
    percent=$(grep " mapped (" "$f" | head -n1 | awk -F'[()%]' '{print $2}')
    echo -e "${sample}\t${total}\t${mapped}\t${percent}" >> mapping_summary.tsv
done
```
<img width="722" height="194" alt="image" src="https://github.com/user-attachments/assets/0b5f5bed-362f-43b3-b20c-16b749fe3b9d" />






