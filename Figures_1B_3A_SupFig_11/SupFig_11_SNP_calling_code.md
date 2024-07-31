---
output:
  html_document: default
  pdf_document: default
---
# SNPs Calling  (Sup. Fig. 11) 

The code used to access the file `officinalis.KSL4-2.merged.snps.snpeff.haplotypecaller.tsv` which is utilized in Sup. Fig. 11, is provided below. Please note that all this code was executed in a different directory, and I am providing only the relative path to each file

## Download the sequencing data from the following accessions: 

```
# accessions_list.txt
W1830	DRR057221
W1270	DRR057213
W1291	DRR057214
W0566	DRR057210
W1252	DRR057212
W1301	DRR057215
W1361	DRR057219
W1200	DRR057211
W1302	DRR057216
W1308	DRR057217
W1315	DRR057218
W1814	DRR057220
W0002	DRR014203
W0065	DRR057209
W1930	DRR057222
```

```{bash}
# sra-tools 2.11.0 
prefetch --option-file accessions_list.txt
for f in DRR*; do fasterq-dump $f; done
```

## Run nextflow pipeline

O. officinalis genome assembly (Ooffic_genome.fa):

```
https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_008326285.1/ 
```

samplesheet.csv:
```
patient,sample,lane,fastq_1,fastq_2
W1200,DRR057211,lane_1,DRR057211_1.fastq.gz,DRR057211_2.fastq.gz
W1930,DRR057222,lane_1,DRR057222_1.fastq.gz,DRR057222_2.fastq.gz
W0002,DRR014203,lane_1,DRR014203_1.fastq.gz,DRR014203_2.fastq.gz
W0065,DRR057209,lane_1,DRR057209_1.fastq.gz,DRR057209_2.fastq.gz
W0566,DRR057210,lane_1,DRR057210_1.fastq.gz,DRR057210_2.fastq.gz
W1252,DRR057212,lane_1,DRR057212_1.fastq.gz,DRR057212_2.fastq.gz
W1270,DRR057213,lane_1,DRR057213_1.fastq.gz,DRR057213_2.fastq.gz
W1291,DRR057214,lane_1,DRR057214_1.fastq.gz,DRR057214_2.fastq.gz
W1301,DRR057215,lane_1,DRR057215_1.fastq.gz,DRR057215_2.fastq.gz
W1302,DRR057216,lane_1,DRR057216_1.fastq.gz,DRR057216_2.fastq.gz
W1308,DRR057217,lane_1,DRR057217_1.fastq.gz,DRR057217_2.fastq.gz
W1315,DRR057218,lane_1,DRR057218_1.fastq.gz,DRR057218_2.fastq.gz
W1361,DRR057219,lane_1,DRR057219_1.fastq.gz,DRR057219_2.fastq.gz
W1814,DRR057220,lane_1,DRR057220_1.fastq.gz,DRR057220_2.fastq.gz
W1830,DRR057221,lane_1,DRR057221_1.fastq.gz,DRR057221_2.fastq.gz
```

```{bash}

# nextflow/22.10.4

nextflow \
  run \
    nf-core/sarek \
    --input samplesheet.csv \
    --fasta Ooffic_genome.fa \
    --tools haplotypecaller \
    --skip_tools baserecalibrator \
    --igenomes_ignore \
    --save_bam_mapped \
    --outdir snp_calling_results \
    -profile biohpc_gen \
    -r 3.0.1
```

## Merging `.vcf` files of individual genotypes

```{bash}
# Create a list of paths of each vcf
find snp_calling_results/variant_calling/haplotypecaller -type f -name "*haplotypecaller.vcf.gz" > officinalis.merged.snps.haplotypecaller.list
```

Merging with bcftools and indexing the merged `.vcf` 

```{bash}
# bcftools 1.13
bcftools merge -l officinalis.merged.snps.haplotypecaller.list -o officinalis.merged.snps.haplotypecaller.vcf.gz
bcftools index officinalis.merged.snps.haplotypecaller.vcf.gz
```

## Annotate snp effects 

1. Build snpeff database

I installed snpeff in a conda environment that I called `snpeff`
```{bash}
# Create directory for this new genome (O. officinalis)
mkdir ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/data/Oryza_officinalis

# Get annotation files
cp Ooffi_maker_gene_annotation.gff ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/data/Oryza_officinalis

# Get the genome reference sequence file
cp Ooffic_genome.fa ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/data/Oryza_officinalis

# change the name of the annotation and fasta files 
cd ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/data/Oryza_officinalis/
mv Ooffi_maker_gene_annotation.gff genes.gff
gzip genes.gff

cd ..
mkdir genomes
mv Oryza_officinalis/Ooffic_genome.fa genomes/Oryza_officinalis.fa
gzip genomes/Oryza_officinalis.fa


cd ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/
```

Open the `snpEff.config` file 


Add the following below `Databases & Genomes`:
It is very important that version matches the names of the folder and genome.
```
# Oryza officinalis genome, version Oryza_officinalis
Oryza_officinalis.genome : Oryza_officinalis
Oryza_officinalis.codonTable : Standard
```

This part of the file will look like:
```
#-------------------------------------------------------------------------------
# Databases & Genomes
#
# One entry per genome version.
#
# For genome version 'ZZZ' the entries look like
#       ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
#       ZZZ.reference           : [Optional] Comma separated list of URL to site/s where information for building ZZZ database was extracted.
#       ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
#
#-------------------------------------------------------------------------------

# Oryza officinalis genome, version Oryza_officinalis
Oryza_officinalis.genome : Oryza_officinalis
Oryza_officinalis.codonTable : Standard
```
2. Add proteins as check point when building the database:

Using the `.fasta` and the `.gff3` annotation files, one can generate the `Ooffi_prot.fa` file, which contains all O. officinalis proteins with `gffread` package (it's anyway provided in this repository)

The >ID of each chromosome in the assembly `.fasta` was adjusted to align with the chromosome names in the annotation file.

```bash
# gffread 0.12.7
gffread -S -y Ooffi_prot.fa -g Ooffic_genome.fa Ooffi_maker_gene_annotation.gff
```

Add the proteome of O. officinalis to O. officinalis database:
```bash
cd ~/miniconda3/envs/snpeff/share/snpeff-5.0-1/data/Oryza_officinalis
cp input/Ooffi_prot.fa protein.fa
gzip protein.fa
```

3. Create the snpEff database for O. officinalis: 
```bash
java -jar snpEff.jar build -gff3 -v Oryza_officinalis
```

4. Run snp eff annotation
```{bash}
SNPEFF=~/miniconda3/envs/snpeff/share/snpeff-5.0-1/snpEff.jar

java -Xmx8G -jar $SNPEFF eff Oryza_officinalis \
  officinalis.merged.snps.haplotypecaller.vcf.gz \
    > officinalis.merged.snps.snpeff.haplotypecaller.vcf
```


## Extract SNPs from the KSL4-2 region 

```{bash}
# Compress
bgzip officinalis.merged.snps.snpeff.haplotypecaller.vcf
# Index
bcftools index officinalis.merged.snps.snpeff.haplotypecaller.vcf.gz
# Extract KSL4-2's region 
bcftools \
  view \
    -r OoffiChr04:7878000-7887000 \
    officinalis.merged.snps.snpeff.haplotypecaller.vcf.gz \
    -o officinalis.KSL4-2.merged.snps.snpeff.haplotypecaller.vcf.gz
# Index   
bcftools index officinalis.KSL4-2.merged.snps.snpeff.haplotypecaller.vcf.gz

# Extract vcf file into a table to further process the snps

bcftools \
  query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%AF\t%INFO\n' \
    officinalis.KSL4-2.merged.snps.snpeff.haplotypecaller.vcf.gz > officinalis.KSL4-2.merged.snps.snpeff.haplotypecaller.tsv
```
