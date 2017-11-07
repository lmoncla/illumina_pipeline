# TCF lab pipeline for Illumina sequence data analysis

# README #

This repository contains the TCF Illumina sequence data analysis pipeline. This pipeline handles going from raw fastq files that came off the MiSeq to generating variant calls and performing diversity analyses.

* TCF Illumina pipeline
* Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

## How do I get set up? ##

This pipeline is written in Python 2.7. To use, download the python file and config file, fill out the config file, and add to your path.

## Dependencies ##

This pipeline is a wrapper that runs a series of publicly available and free softwares. All of these softwares are already installed in /usr/local/bin on the Coltrane server and are ready to use if the analyses are run on any part of your account on Coltrane. It relies on the following dependencies:

#### Trimmomatic
http://www.usadellab.org/cms/?page=trimmomatic
Trimmomatic performs fastq file trimming based on quality scores and removes short reads.

#### bowtie2
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Bowtie2 maps trimmed reads to a reference sequence, producing an output sam file.

#### picard
http://broadinstitute.github.io/picard/
Picard is a commonly used software for removing duplicate reads with the MarkDuplicates tool.

#### lofreq  
http://csb5.github.io/lofreq/
Lofreq calls variants and filters them based on quality scores, forward/reverse read balance, coverage, frequency, quality, etc..

#### varscan
http://varscan.sourceforge.net/
Varscan is another option for a variant caller. Allows variants to be called/filtered based on the same criteria as Lofreq (forward/reverse read balance, coverage, frequency, quality)

#### samtools
http://samtools.sourceforge.net/
Necessary for converting sam/bam files, sorting bam files, removing low mapping quality reads from sam files, and generating pileup files

#### Trinity
https://github.com/trinityrnaseq/trinityrnaseq/wiki
A de novo assembler written for RNA seq data. Used in this pipeline for de novo assembling trimmed reads to identify any unexpected contaminant contigs.

#### bbmap  
https://sourceforge.net/projects/bbmap/
Used for extracting mapped fastq files from sam/bam files for input into Trinity de novo assembly for quality control. Also used for coverage normalization, with the tool bbnorm.sh.  

#### BLAST command line tools
https://www.ncbi.nlm.nih.gov/books/NBK279690/
Contigs generated from Trinity assembly are piped to BLAST to return the top 10 hits that most closely match the contig sequences.

#### snpEff
http://snpeff.sourceforge.net/
snpEff is used for annotating output vcf files with coding region coordinates and inferring amino acid changes from nucleotide variant calls.

#### popoolation
https://sourceforge.net/p/popoolation/wiki/Main/
Popoolation is used for calculating pi, piN and piS across full genes, genomes or in sliding windows.


## Basic usage: ##
`illumina_pipeline_fastq_to_snps.py`

#### input files
As this is currently written, the script will run on all of the fastq files that are in your current directory. These files should all be standard fastq files and need to end in .fastq. The program will automatically combine forward and reverse reads that derive from the same sample together into the same folder. The pipeline has been written to allow you to either map all samples to the same reference sequence (as you would want to do for any sort of experimental evolution/infection study) or to map each sample to a unique reference (as you might want to do for clinical samples). The main difference that you need to worry for specifying between these 2 options is in the reference sequence section of the config file (see: Filling in the config file below).

#### the config file
All parameters including values you alter in the various programs, which analyses you wish to perform, and your reference sequence are specified in the config file. The config.py file contains annotated notes about which fields to fill in and what they are, although you should always consult the manuals for the specific programs you are calling for more specific details about what they do. These arguments will then be passed to illumina_pipeline_fastq_to_snps.py, which will incorporate those arguments at the appropriate time in the pipeline. After running, the program will produce a .params file which will contain all of the arguments that you used for the analyses. The intention here is that you shouldn't have to alter the illumina_pipeline_fastq_to_snps.py file, and should be able to just alter the config.py file. Because you will have an output params file, it will let you keep a record of what you ran so that you can directly edit the config.py file and repeat the analyses later if you want.

#### a note on de novo assembly
This pipeline performs de novo assembly with Trinity, which is a de novo assembler originally developed for RNAseq data. It does not perform well for de novo assembling RNA genomes and results in multiple contigs per genome. However, unlike IVA, it will detect even low-frequency contigs that are contaminant sequences (like bacterial, human, dog) that are inevitably always in a sequencing run. The point of including it in this pipeline is for quality control: if you perform a de novo assembly, and BLAST the result, and the BLAST result shows that you have an extra viral sequence that shouldn't be there, then you should be concerned. More specifically, if you do find that you are picking up on a strange flu contig that should not be there, you should extract that contig and map your reads back to it to determine how many contaminant reads there are. I have written in this de novo assembly step so that you can de novo assembly all reads or only the ones that were mapped to the reference. I personally find this second option a little better because it takes less time, and my main concern if for contaminant reads in the mapped file that is used for downstream analysis.

## Filling in the config file ##

### SECTION 1: SPECIFY WHICH TASKS YOU WANT TO DO HERE
You may elect to perform trimming, mapping, SNP calling, duplicate read removal, coverage normalization, de novo assembly and run popoolation using this pipeline. To enable these analyses, simply type "True" (make sure to use a capital T) after the = each option. The config file is divided into 3 basic sections: basic tasks (trimming, mapping, calling and annotating SNPs and de novo assembly), cleaning tasks, which are performed after mapping on sam/bam files (coverage normalization and duplicate read removal), and running popoolation. Specifics are below:

#### Basic tasks
#### self.trim = `True` or `False`
Use Trimmomatic to trim the ends of your reads. This must be done for all raw fastq files. If set to True, fill in the parameters under the "SET TRIMMING PARAMETERS" section.

#### self.map = `True` or `False`
Use bowtie2 to map to a reference sequence. If set to True, fill in the "SPECIFY REFERENCE SEQUENCE" section.

#### self.call_snps = `True` or `False`
Use either Varscan or Lofreq to call SNPs after mapping to a reference sequence. You must map to a reference before calling SNPs, as the input file for SNP calling is output file for mapping. If set to True, fill in the "SET SNP CALLING PARAMETERS" section of the config file.

#### self.annotate_aa_changes = `True` or `False`
Use SNPEff to annotate coding region changes after variant calling. This will take as input the vcf file that has been produced from self.call_snps and use it as input into SNPEff for annotation. This requires that your genomes are indexed with SNPEff. There are instructions for this on the SNPEff home page as well as some notes on the bottom of this page.

#### self.de_novo_assembly = `True` or `False`
Use Trinity to de novo assemble all trimmed fastq files. All results from Trinity will be output to a folder called "trinity_output". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt".

#### self.de_novo_assemble_mapped_reads = `True` or `False`
Extract all mapped reads from the sam file, and then use Trinity to de novo assemble all of those extracted, trimmed, mapped reads. All results from Trinity will be output to a folder called "trinity_de_novo_assembly_mapped_reads_only". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt".

#### Data cleaning tasks
#### self.remove_duplicate_reads = `True` or `False`
Use picard's MarkDuplicates tool to remove duplicate reads. Picard considers reads to be duplicates if they have the exact same 5' start site. In some comparisons with samblaster, picard and dedupe, picard performed the best. Picard improved data reproducibility for within-host influenza data, while samblaster and dedupe resulted in an increase in the number of low-frequency variants called, but did result in consistent calls. Therefore, I have elected to use picard for duplicate read removal in this pipeline. Duplicate read removal will occur after mapping, but before coverage normalization.

#### self.normalize_coverage = `True` or `False`
Use bbnorm.sh from the bbmap software package to normalize coverage across the genome. This is performed after mapping. If duplicate read removal is also to be performed, normalization will occur after duplicate read removal. BBnorm.sh normalizes coverage using fastq files as input, so when self.normalize_coverage is set to `True`, reads will be extracted from the sam file, used as input for for bbnorm, and then remapped. The resulting output file will end in `.normalized.coverage.sam` where `coverage` will be the desired coverage depth that you set.

#### self.coverage_normalization_depth = `integer`
Set the desired depth of coverage you wish to achieve after coverage normalization to `integer`.
=======

### SECTION 2 : SET/ALTER PARAMETERS

#### SET TRIMMING PARAMETERS:
#### self.minlength = `integer`
After read trimming has been performed, discard reads that are shorter than integer length. Value must be an integer. I would recommend 100.

#### self.window_size = `integer`
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by `integer`.

#### self.trim_qscore = `integer`
Phred-based quality score threshold to use during trimming. If you would like to use a Q30 threshold, you would specify 30. 30 is recommended.

#### self.minlength = `integer`
After read trimming has been performed, discard reads that are shorter than `integer` length. Value must be an integer. I would recommend 100.

#### self.window_size = `integer`
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by `integer`.

#### self.trim_qscore = `integer`
Phred-based quality score threshold to use during trimming. If you would like to use a Q30 threshold, you would specify 30. 30 is recommended.

=======

#### SPECIFY REFERENCE SEQEUNCE AND MAPPING QUALITY:
One important note here is that this pipeline is meant to run with a single reference sequence file. If you want to specify multiple gene segments, simply put all of them into the same fasta file. The fasta file must end in .fasta or .fa.

#### self.use_different_reference_for_each_sample = `True` or `False`
Specify True to map all of the samples to the same reference sequence or False to map each sample to it's own reference. If specifying False, then you need to put the fasta reference file into the same folder as the trimmed fastqs. The easiest way to do this is to run the pipeline and do only Trimming, which will combine the forward and reverse fastq files and make folders with their specific names. Then just move the fasta reference files into the appropriate folder.

#### self.reference_sequence = `path to reference sequence`
If you are mapping everything to the same reference sequence, then you have to specify the full path to the reference sequence you wish to use. The reference sequence should be in fasta format and can end in .fasta or .fa. Ex: User/Documents/CA04_HA.fa

#### self.reference_sequence_name = `name`
Bowtie2 requires the user to input a "base name" for the reference sequence. For this, specify the actual name of the reference sequence, NOT the path. Ex: CA04_HA.fasta

#### self.mapping_quality_threshold = `30`
After mapping with bowtie, it is a good idea to remove reads from the sam file that have a low mapping quality score. Reads with mapping quality scores less than this value will be removed from the sam file. If you do not wish to use this option, set to 0. These values are specified as Phred scores. 

=======
#### SET SNP CALLING PARAMETERS:
One important note here is that if you would like SNPs to be annotated as to whether they cause a coding region change, then you need to put together gtf files and configure new genomes in snpEff. Instructions for how to do that are at the end of this document.

#### self.use_lofreq = `True` or `False`
#### self.use_varscan = `True` or `False`
For each, set to True to call SNPs with that program. The pipeline can be run using either, neither or both. All output files will be written to a sub-folder called "snp_calls".

#### self.min_coverage = `integer`
This will set the minimum coverage required at a base in order to perform variant calling at that base. Variants at positions with coverage less than integer will not be called.

#### self.snp_qual_threshold = `integer`
This will set the minimum quality score required at a base in order to perform variant calling at that base. Variants with quality scores lower than integer will not be called.

#### self.snp_frequency = `float`
Variants that are present at a frequency less than that set by decimal will not be called. decimal values should range from 0 to 1, with a value of 0.01 specifying that SNPs should be called at a 1% frequency cutoff.

=======

### Output
After this has been run, a folder will be made for each sample, which will contain the original fastq files, trimmed fastq files, mapping files in sam and bam format, variant calls in vcf format, and a parameters file. Details are below:

#### sample.fastq
the original fastq files

#### sample.trimmed.fastq
trimmed fastq files

#### sample.sam
trimmed fastq files that were mapped to the reference sequence

#### sample.bam
trimmed fastq files that were mapped to the reference sequence, but in bam format

#### sample.sorted.bam
trimmed fastq files that were mapped to the reference sequence, but in sorted bam format (necessary for SNP calling)

#### sample.lofreq.`snp_frequency`.vcf
unfiltered lofreq variant calls. These have NOT been filtered to account for minimum quality, coverage, or SNP frequency. Here, `snp_frequency` specifies the minimum frequency a variant had to be present to be called that was applied to filtering.

#### sample.lofreq.filtered.`snp_frequency`.vcf
filtered variant calls, filtered with the parameters you specified in the config file.

#### sample.lofreq.annotated.`snp_frequency`.vcf
filtered lofreq variant calls.filtered with the parameters you specified in the config file and annotated by snpEff to include information about amino acid changes to coding regions.

#### sample.varscan.snps.`snp_frequency`.vcf
variants called by varscan according to the parameters you specified in the config file. Not yet annotated.

#### sample.varscan.annotated.snps.`snp_frequency`.vcf
variants called by varscan according to the parameters you specified in the config file and annotated by snpEff to include information about amino acid changes to coding regions.

#### trinity_output
A folder containing the output files for Trinity when run using all reads as input. The output contigs are in a file called Trinity.fasta and the restuls from sending those contigs to BLAST are contained in Trinity_BLAST_result.txt.

#### trinity_de_novo_assembly_mapped_reads_only
A folder that contains the output files for Trinity when run using only mapped reads as input. The output contigs are in a file called Trinity.fasta and the restuls from sending those contigs to BLAST are contained in Trinity_BLAST_result.txt.

#### sample.params
a parameters file containing a summary of the commands you specified


=======
## Configuring new genomes in snpEff ##

Unfortunately, putting together coding regions for each genome is a little tedious, especially if you are mapping each sample to its own reference and need to annotate a new set of coordinates for each sample. There is no very good way that I have thought of to automate this, so we are unfortunately left with these fairly specific and tedious instructions for building new genome databases with snpEff for each new genome you are using.

snpEff is a program that annotates vcf files. It functions by using genome databases that are publicly available to determine coding region coordinate. Unfortunately, this is strongly skewed towards mammalian genomes and there are only a handful of viral sequences that are already prepared. It is therefore very likely that you will need to put together your own genomic coordinates and configure new genomes for each sample. Here is how to do it.

Detailed instructions are provided on snpEff's website: http://snpeff.sourceforge.net/SnpEff_manual.html#databases

However, there are a few important quirks that are not specified in this document, so here are my notes for how to do this.

#### A. Configure a new genome by adding to the the snpEff.config file.
1. Open this file, which is in the snpEff folder in any text editor.
2. Scroll to the very bottom of the file
3. Add your new reference sequence to this file using the following format: `reference_sequence_name.genome : nickname for that reference sequence`. Here, `reference_sequence_name` is the full name of the fasta file you are using and the nickname is whatever you want it to be. The .genome is necessary. snpEff will not correctly make the reference genome if this does not end in .genome.
4. save the altered configuration file


#### B. Build the database using a gtf file. This is recommended.
1. cd into snpEff/data
2. make a new directory and name is `reference_sequence_name`. The name of the folder should be exactly identical to `reference_sequence_name`. It should NOT include the .genome.
3. place your genome fasta file and the associated gtf file into the reference_sequence_name folder. The reference sequence fasta file MUST be named sequences.fa. The gtf MUST be named genes.gtf.
4. run the following command: `java -jar snpEff.jar build -gtf22 -v reference_sequence_name`

* Note: when you run the build command, reference_sequence_name does NOT include the .genome. If you have multiple gene segments whose coordinates are all specified in 1 gtf, that is fine. They just need to all have the same base names in the gtf and sequences.fa file (so CA04_HA and CA04_NA is fine).

* Once these steps are complete, your genome should be built in a zipped file called snpEffectPredictor.bin

#### C. Annotate your vcfs
Now, simply run the command: `java -jar snpEff.jar reference_sequence_name input.vcf > output.vcf`



### Questions and comments ###

Louise Moncla
lhmoncla@gmail.com
