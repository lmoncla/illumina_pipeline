# README #

This repository contains the TCF Illumina sequence data analysis pipeline. This pipeline handles going from raw fastq files that came off the MiSeq to generating variant calls and performing diversity analyses. 

* TCF Illumina pipeline
* Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

## How do I get set up? ##

This pipeline is written in Python 2.7. To use, download the python file and config file, fill out the config file, and add to your path. 

## Basic usage: ##
illumina_pipeline_fastq_to_snps.py config

## Dependencies ##

This pipeline is a wrapper that runs a series of publicly available and free softwares. All of these softwares are already installed in /usr/local/bin on the Coltrane server and are ready to use if the analyses are run on any part of your account on Coltrane. It relies on the following dependencies: 

#### Trimmomatic
http://www.usadellab.org/cms/?page=trimmomatic
Trimmomatic performs fastq file trimming based on quality scores and removes short reads. 

#### bowtie2
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Bowtie2 maps trimmed reads to a reference sequence, producing an output sam file. 

#### lofreq  
http://csb5.github.io/lofreq/
Lofreq calls variants and filters them based on quality scores, forward/reverse read balance, coverage, frequency, quality, etc.. 

#### varscan
http://varscan.sourceforge.net/
Varscan is another option for a variant caller. Allows variants to be called/filtered based on the same criteria as Lofreq (forward/reverse read balance, coverage, frequency, quality)

#### samtools
http://samtools.sourceforge.net/
Necessary for converting sam/bam files, sorting bam files and generating pileup files

#### Trinity
https://github.com/trinityrnaseq/trinityrnaseq/wiki
A de novo assembler written for RNA seq data. Used in this pipeline for de novo assembling trimmed reads to identify any unexpected contaminant contigs. 

#### bbmap  
https://sourceforge.net/projects/bbmap/
Used for extracting mapped fastq files from sam/bam files for input into Trinity de novo assembly for quality control. 

#### BLAST command line tools
https://www.ncbi.nlm.nih.gov/books/NBK279690/
Contigs generated from Trinity assembly are piped to BLAST to return the top 10 hits that most closely match the contig sequences. 

#### snpEff 
http://snpeff.sourceforge.net/
snpEff is used for annotating output vcf files with coding region coordinates and inferring amino acid changes from nucleotide variant calls. 

#### popoolation
https://sourceforge.net/p/popoolation/wiki/Main/
Popoolation is used for calculating pi, piN and piS across full genes, genomes or in sliding windows. 


### Contribution guidelines ###


### Who do I talk to? ###

* Louise Moncla lhmoncla@gmail.com