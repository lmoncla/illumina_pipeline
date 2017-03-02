# README #

This repository contains the TCF Illumina sequence data analysis pipeline. This pipeline handles going from raw fastq files that came off the MiSeq to generating variant calls and performing diversity analyses. 

* TCF Illumina pipeline
* Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
This pipeline is written in Python 2.7. Simply download the python file and configuration file and run them. 

* Dependencies

This pipeline is a wrapper that runs a series of publicly available and free softwares. It relies on the following dependencies: 

1. Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic
2. bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
3. lofreq - http://csb5.github.io/lofreq/
4. varscan - http://varscan.sourceforge.net/
5. samtools - http://samtools.sourceforge.net/
6. Trinity - https://github.com/trinityrnaseq/trinityrnaseq/wiki
7. BLAST command line tools - https://www.ncbi.nlm.nih.gov/books/NBK279690/
8. snpEff - http://snpeff.sourceforge.net/
9. popoolation - https://sourceforge.net/p/popoolation/wiki/Main/
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact