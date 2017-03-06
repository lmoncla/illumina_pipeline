 HEAD
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


## Basic usage: ##
`illumina_pipeline_fastq_to_snps.py config`

#### input files
As this is currently written, the script will run on all of the fastq files that are in your current directory. These files should all be standard fastq files and need to end in .fastq. The program will automatically combine forward and reverse reads that derive from the same sample together into the same folder. The pipeline has been written to allow you to either map all samples to the same reference sequence (as you would want to do for any sort of experimental evolution/infection study) or to map each sample to a unique reference (as you might want to do for clinical samples). The main difference that you need to worry for specifying between these 2 options is in the reference sequence section of the config file (see: Filling in the config file below). 

#### the config file
All parameters including values you alter in the various programs, which analyses you wish to perform, and your reference sequence are specified in the config file. The config.py file contains annotated notes about which fields to fill in and what they are, although you should always consult the manuals for the specific programs you are calling for more specific details about what they do. These arguments will then be passed to illumina_pipeline_fastq_to_snps.py, which will incorporate those arguments at the appropriate time in the pipeline. After running, the program will produce a .params file which will contain all of the arguments that you used for the analyses. The intention here is that you shouldn't have to alter the illumina_pipeline_fastq_to_snps.py file, and should be able to just alter the config.py file. Because you will have an output params file, it will let you keep a record of what you ran so that you can directly edit the config.py file and repeat the analyses later if you want. 

#### a note on de novo assembly
This pipeline performs de novo assembly with Trinity, which is a de novo assembler originally developed for RNAseq data. It does not perform well for de novo assembling RNA genomes and results in multiple contigs per genome. However, unlike IVA, it will detect even low-frequency contigs that are contaminant sequences (like bacterial, human, dog) that are inevitably always in a sequencing run. The point of including it in this pipeline is for quality control: if you perform a de novo assembly, and BLAST the result, and the BLAST result shows that you have an extra viral sequence that shouldn't be there, then you should be concerned. More specifically, if you do find that you are picking up on a strange flu contig that should not be there, you should extract that contig and map your reads back to it to determine how many contaminant reads there are. I have written in this de novo assembly step so that you can de novo assembly all reads or only the ones that were mapped to the reference. I personally find this second option a little better because it takes less time, and my main concern if for contaminant reads in the mapped file that is used for downstream analysis. 

## Filling in the config file ##

### SECTION 1: SPECIFY WHICH TASKS YOU WANT TO DO HERE
You may elect to perform trimming, mapping, SNP calling, de novo assembly and run popoolation using this pipeline. To enable these analyses, simply type "True" (make sure to use a capital T) after the = each option. Specifics are specified below: 

#### self.trim = `True` or `False`
Use Trimmomatic to trim the ends of your reads. This must be done for all raw fastq files. If set to True, fill in the parameters under the "SET TRIMMING PARAMETERS" section. 

#### self.map = `True` or `False`
Use bowtie2 to map to a reference sequence. If set to True, fill in the "SPECIFY REFERENCE SEQUENCE" section. 

#### self.call_snps = `True` or `False` 
Use either Varscan or Lofreq to call SNPs after mapping to a reference sequence. You must map to a reference before calling SNPs, as the input file for SNP calling is output file for mapping. If set to True, fill in the "SET SNP CALLING PARAMETERS" section of the config file. All SNP calls will be output to a folder called "snp_calls".

#### self.de_novo_assembly = `True` or `False`
Use Trinity to de novo assemble all trimmed fastq files. All results from Trinity will be output to a folder called "trinity_output". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt". 

#### self.de_novo_assemble_mapped_reads = `True` or `False`
=======
#### self.trim = `True` or `False`
Use Trimmomatic to trim the ends of your reads. This must be done for all raw fastq files. If set to True, fill in the parameters under the "SET TRIMMING PARAMETERS" section. 

#### self.map =`True` or `False`
Use bowtie2 to map to a reference sequence. If set to True, fill in the "SPECIFY REFERENCE SEQUENCE" section. 

#### self.call_snps = `True` or `False` 
Use either Varscan or Lofreq to call SNPs after mapping to a reference sequence. You must map to a reference before calling SNPs, as the input file for SNP calling is output file for mapping. If set to True, fill in the "SET SNP CALLING PARAMETERS" section of the config file. All SNP calls will be output to a folder called "snp_calls".

#### self.de_novo_assembly = `True` or `False`
Use Trinity to de novo assemble all trimmed fastq files. All results from Trinity will be output to a folder called "trinity_output". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt". 

#### self.de_novo_assemble_mapped_reads = `True` or `False`
Extract all reads that were mapped to the reference and use Trinity to perform a de novo assembly on those reads. All results from Trinity will be output to a folder called "trinity_de_novo_assembly_mapped_reads_only". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt." 



### SECTION 2 : SET/ALTER PARAMETERS 

#### SET TRIMMING PARAMETERS:
#### self.minlength = `integer` 
After read trimming has been performed, discard reads that are shorter than integer length. Value must be an integer. I would recommend 100. 

#### self.window_size = `integer`
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by integer. 

#### self.trim_qscore = `integer`
=======
#### self.minlength = `integer` 
After read trimming has been performed, discard reads that are shorter than `integer` length. Value must be an integer. I would recommend 100. 

#### self.window_size = `integer`
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by `integer`. 

#### self.trim_qscore = `integer`
Phred-based quality score threshold to use during trimming. If you would like to use a Q30 threshold, you would specify 30. 30 is recommended. 




#### SPECIFY REFERENCE SEQEUNCE:
One important note here is that this pipeline is meant to run with a single reference sequence file. If you want to specify multiple gene segments, simply put all of them into the same fasta file. The fasta file must end in .fasta or .fa. 

#### self.use_different_reference_for_each_sample = True or False
Specify True to map all of the samples to the same reference sequence or False to map each sample to it's own reference. If specifying False, then you need to put the fasta reference file into the same folder as the trimmed fastqs. The easiest way to do this is to run the pipeline and do only Trimming, which will combine the forward and reverse fastq files and make folders with their specific names. Then just move the fasta reference files into the appropriate folder. 

#### self.reference_sequence = path to reference sequence
If you are mapping everything to the same reference sequence, then you have to specify the full path to the reference sequence you wish to use. The reference sequence should be in fasta format and can end in .fasta or .fa. Ex: User/Documents/CA04_HA.fa

#### self.reference_sequence_name = name
=======
#### self.use_different_reference_for_each_sample = `True` or `False`
Specify True to map all of the samples to the same reference sequence or False to map each sample to it's own reference. If specifying False, then you need to put the fasta reference file into the same folder as the trimmed fastqs. The easiest way to do this is to run the pipeline and do only Trimming, which will combine the forward and reverse fastq files and make folders with their specific names. Then just move the fasta reference files into the appropriate folder. 

#### self.reference_sequence = `path to reference sequence`
If you are mapping everything to the same reference sequence, then you have to specify the full path to the reference sequence you wish to use. The reference sequence should be in fasta format and can end in .fasta or .fa. Ex: User/Documents/CA04_HA.fa

#### self.reference_sequence_name = `name`
Bowtie2 requires the user to input a "base name" for the reference sequence. For this, specify the actual name of the reference sequence, NOT the path. Ex: CA04_HA.fasta



#### SET SNP CALLING PARAMETERS:
One important note here is that if you would like SNPs to be annotated as to whether they cause a coding region change, then you need to put together gtf files and configure new genomes in snpEff. Instructions for how to do that are at the end of this document. 

#### self.use_lofreq = True or False
#### self.use_varscan = True or False
For each, set to True to call SNPs with that program. The pipeline can be run using either, neither or both. All output files will be written to a sub-folder called "snp_calls". 

#### self.min_coverage = integer 
This will set the minimum coverage required at a base in order to perform variant calling at that base. Variants at positions with coverage less than integer will not be called. 

#### self.snp_qual_threshold = integer
This will set the minimum quality score required at a base in order to perform variant calling at that base. Variants with quality scores lower than integer will not be called. 

#### self.snp_frequency = decimal
Variants that are present at a frequency less than that set by decimal will not be called. decimal values should range from 0 to 1, with a value of 0.01 specifying that SNPs should be called at a 1% frequency cutoff. 
=======
#### self.use_lofreq = `True` or `False`
#### self.use_varscan = `True` or `False`
For each, set to True to call SNPs with that program. The pipeline can be run using either, neither or both. All output files will be written to a sub-folder called "snp_calls". 

#### self.min_coverage = `integer`
This will set the minimum coverage required at a base in order to perform variant calling at that base. Variants at positions with coverage less than `integer` will not be called. 

#### self.snp_qual_threshold = `integer`
This will set the minimum quality score required at a base in order to perform variant calling at that base. Variants with quality scores lower than `integer` will not be called. 

#### self.snp_frequency = `integer`
Variants that are present at a frequency less than that set by `integer` will not be called. `integer` values should range from 0 to 1, with a value of 0.01 specifying that SNPs should be called at a 1% frequency cutoff. 



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

#### sample.lofreq.(frequency_cutoff).vcf
unfiltered lofreq variant calls. These have NOT been filtered to account for minimum quality, coverage, or SNP frequency. Here, (frequency_cutoff) specifies the minimum frequency a variant had to be present to be called that was applied to filtering. 

#### sample.lofreq.filtered.(frequency_cutoff).vcf
filtered variant calls, filtered with the parameters you specified in the config file.

#### sample.lofreq.annotated.(frequency_cutoff).vcf
filtered lofreq variant calls.filtered with the parameters you specified in the config file and annotated by snpEff to include information about amino acid changes to coding regions. 

#### sample.varscan.snps.(frequency_cutoff).vcf
variants called by varscan according to the parameters you specified in the config file. Not yet annotated. 

#### sample.varscan.annotated.snps.(frequency_cutoff).vcf
variants called by varscan according to the parameters you specified in the config file and annotated by snpEff to include information about amino acid changes to coding regions. 

#### trinity_output
A folder containing the output files for Trinity when run using all reads as input. The output contigs are in a file called Trinity.fasta and the restuls from sending those contigs to BLAST are contained in Trinity_BLAST_result.txt.

#### trinity_de_novo_assembly_mapped_reads_only
A folder that contains the output files for Trinity when run using only mapped reads as input. The output contigs are in a file called Trinity.fasta and the restuls from sending those contigs to BLAST are contained in Trinity_BLAST_result.txt.
 
#### sample.params
a parameters file containing a summary of the commands you specified




### Questions and comments ###

Louise Moncla lhmoncla@gmail.com
=======
# illumina_pipeline
TCF lab pipeline for Illumina sequence data analysis
 db1bde78c999a8649bc64adba3e3f241cdac1633
=======
