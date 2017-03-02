README for illumina_pipeline_fastq_to_snps.py

January 27, 2016

#########################################################################################
WHAT:
#########################################################################################
This README file describes how to use the illumina_piepline_fastq_to_snps.py, what it does, and which requirements it needs. This script is simply a wrapper for performing trimming, mapping to reference sequence, variant calling, and diversity analysis on pooled virus sequence data from Illumina sequencing experiments. I have elected to use Trimmomatic for trimming, bowtie2 for mapping, either lofreq or varscan for variant calling, and snpeff for annotation/coding region changes. I have also added in de novo assembly with Trinity, followed by a BLAST search for the contigs produced as a quality control step. 


#########################################################################################
REQUIREMENTS/DEPENDENCIES:
#########################################################################################
1. Trimmomatic - http://www.usadellab.org/cms/?page=trimmomatic
2. bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
3. lofreq - http://csb5.github.io/lofreq/
4. varscan - http://varscan.sourceforge.net/
5. samtools - http://samtools.sourceforge.net/
6. Trinity - https://github.com/trinityrnaseq/trinityrnaseq/wiki
7. BLAST command line tools - https://www.ncbi.nlm.nih.gov/books/NBK279690/
8. snpEff - http://snpeff.sourceforge.net/
9. popoolation - https://sourceforge.net/p/popoolation/wiki/Main/

If you are using the TCF lab Coltrane server, then I have already installed all of these programs in /usr/local/bin, which is a directory that is automatically in your path. I have also added the absolute path of "/usr/local/bin/" for each program that I call in the pipeline. If you are not using this on the TCF Coltrane server, all of these programs should be downloaded and added to your path or installed in a directory that is already in your path. 

When you download the illumina pipeline package, make sure that it includes the python file, "config.py". Place these 2 scripts in a folder that is already in your path or add them to your path. The following are instructions for adding these to your path if you have not done it before:

1. Open terminal
2. Type "nano .bash_profile". The .bash_profile file tells your computer which directories terminal should look through to find files/scripts. Nano is a text editor that will open and display the contents of the file you specify. So you will now see the contents of .bash_profile.
3. Add a new directory to your path. In the screen, type: "export PATH="$PATH:$HOME/directory" with directory being the name of the folder you want to add to your path. 
4. Exit by pressing control + X
5. Save changes by pressing "y". Press enter to keep the same name. Do not change the name. 

If you encounter a permission denied error when you try to run illumina_pipeline_fastq_to_snps.py, navigate to the directory it is housed in and give it the proper permissions by typing in terminal "chmod u+x illumina_pipeline_fastq_to_snps.py" and "chmod u+x config.py". Restart terminal. It should now work. 

Alternatively, if you do not want to add these to your path, you don't have to! But instead of running:
illumina_pipeline_fastq_to_snps.py config
you will need to run:
python illumina_pipeline_fastq_to_snps.py config


#########################################################################################
USAGE:
#########################################################################################
Navigate to a directory that houses your fastq files and a reference sequence. Then run:

illumina_pipeline_fastq_to_snps.py config

## Input files: 
As this is currently written, the script will run on all of the fastq files that are in your current directory. These files should all be standard fastq files and need to end in .fastq. The program will automatically combine forward and reverse reads that derive from the same sample together into the same folder. The pipeline has been written to allow you to either map all samples to the same reference sequence (as you would want to do for any sort of experimental evolution/infection study) or to map each sample to a unique reference (as you might want to do for clinical samples). The main difference that you need to worry for specifying between these 2 options is in the reference sequence section of the config file (see: Filling in the config file below). 

## The config file: 
All parameters including values you alter in the various programs, which analyses you wish to perform, and your reference sequence are specified in the config file. The config.py file contains annotated notes about which fields to fill in and what they are, although you should always consult the manuals for the specific programs you are calling for more specific details about what they do. These arguments will then be passed to illumina_pipeline_fastq_to_snps.py, which will incorporate those arguments at the appropriate time in the pipeline. After running, the program will produce a .params file which will contain all of the arguments that you used for the analyses. The intention here is that you shouldn't have to alter the illumina_pipeline_fastq_to_snps.py file, and should be able to just alter the config.py file. Because you will have an output params file, it will let you keep a record of what you ran so that you can directly edit the config.py file and repeat the analyses later if you want. 

## De novo assembly: 
This pipeline performs de novo assembly with Trinity, which is a de novo assembler originally developed for RNAseq data. It does not perform well for de novo assembling RNA genomes and results in multiple contigs per genome. However, unlike IVA, it will detect even low-frequency contigs that are contaminant sequences (like bacterial, human, dog) that are inevitably always in a sequencing run. The point of including it in this pipeline is for quality control: if you perform a de novo assembly, and BLAST the result, and the BLAST result shows that you have an extra viral sequence that shouldn't be there, then you should be concerned. More specifically, if you do find that you are picking up on a strange flu contig that should not be there, you should extract that contig and map your reads back to it to determine how many contaminant reads there are. I have written in this de novo assembly step so that you can de novo assembly all reads or only the ones that were mapped to the reference. I personally find this second option a little better because it takes less time, and my main concern if for contaminant reads in the mapped file that is used for downstream analysis. 


#########################################################################################
FILLING IN THE CONFIG FILE:
#########################################################################################
The config.py file consists of 2 major sections: specifying which analyses/tasks you want to do and specifying how you want them performed (set/alter paramters). 


### SECTION 1: SPECIFY WHICH TASKS YOU WANT TO DO HERE ##################################

You can elect to perform trimming, mapping, SNP calling, de novo assembly and run popoolation using this pipeline. To enable these analyses, simply type "True" (make sure to use a capital T) after the = each option. Specifics are specified below: 

self.trim = <True> or <False>
Use Trimmomatic to trim the ends of your reads. This must be done for all raw fastq files. If set to True, fill in the parameters under the "SET TRIMMING PARAMETERS" section. 

self.map = <True> or <False>
Use bowtie2 to map to a reference sequence. If set to True, fill in the "SPECIFY REFERENCE SEQUENCE" section. 

self.call_snps = <True> or <False> 
Use either Varscan or Lofreq to call SNPs after mapping to a reference sequence. You must map to a reference before calling SNPs, as the input file for SNP calling is output file for mapping. If set to True, fill in the "SET SNP CALLING PARAMETERS" section of the config file. All SNP calls will be output to a folder called "snp_calls".

self.de_novo_assembly = <True> or <False>
Use Trinity to de novo assemble all trimmed fastq files. All results from Trinity will be output to a folder called "trinity_output". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt". 

self.de_novo_assemble_mapped_reads = <True> or <False>
Extract all reads that were mapped to the reference and use Trinity to perform a de novo assembly on those reads. All results from Trinity will be output to a folder called "trinity_de_novo_assembly_mapped_reads_only". Within that folder, output contigs will be written to Trinity.fasta. Those contigs will be piped to the BLAST server and the top 10 BLAST hits for each contig are reported in the output file "Trinity_BLAST_result.txt." 



### SECTION 2 : SET/ALTER PARAMETERS ####################################################


# SET TRIMMING PARAMETERS:
self.minlength = <integer> 
After read trimming has been performed, discard reads that are shorter than <integer> length. Value must be an integer. I would recommend 100. 

self.window_size = <integer>
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by <integer>. 

self.trim_qscore = <integer>
Phred-based quality score threshold to use during trimming. If you would like to use a Q30 threshold, you would specify 30. 30 is recommended. 


# SPECIFY REFERENCE SEQEUNCE:

One important note here is that this pipeline is meant to run with a single reference sequence file. If you want to specify multiple gene segments, simply put all of them into the same fasta file. The fasta file must end in .fasta or .fa. 

self.use_different_reference_for_each_sample = <True> or <False>
Specify True to map all of the samples to the same reference sequence or False to map each sample to it's own reference. If specifying False, then you need to put the fasta reference file into the same folder as the trimmed fastqs. The easiest way to do this is to run the pipeline and do only Trimming, which will combine the forward and reverse fastq files and make folders with their specific names. Then just move the fasta reference files into the appropriate folder. 

self.reference_sequence = <path to reference sequence>
If you are mapping everything to the same reference sequence, then you have to specify the full path to the reference sequence you wish to use. The reference sequence should be in fasta format and can end in .fasta or .fa. Ex: User/Documents/CA04_HA.fa

self.reference_sequence_name = <name>
Bowtie2 requires the user to input a "base name" for the reference sequence. For this, specify the actual name of the reference sequence, NOT the path. Ex: CA04_HA.fasta


# SET SNP CALLING PARAMETERS:

One important note here is that if you would like SNPs to be annotated as to whether they cause a coding region change, then you need to put together gtf files and configure new genomes in snpEff. Instructions for how to do that are at the end of this document. 

self.use_lofreq = <True> or <False>
self.use_varscan = <True> or <False>
For each, set to True to call SNPs with that program. The pipeline can be run using either, neither or both. All output files will be written to a sub-folder called "snp_calls". 

self.min_coverage = <integer> 
This will set the minimum coverage required at a base in order to perform variant calling at that base. Variants at positions with coverage less than <integer> will not be called. 

self.snp_qual_threshold = <integer>
This will set the minimum quality score required at a base in order to perform variant calling at that base. Variants with quality scores lower than <integer> will not be called. 

self.snp_frequency = <decimal>
Variants that are present at a frequency less than that set by <decimal> will not be called. <decimal> values should range from 0 to 1, with a value of 0.01 specifying that SNPs should be called at a 1% frequency cutoff. 



#########################################################################################
OUTPUT: 
#########################################################################################
After this has been run, a folder will be made for each sample, which will contain the original fastq files, trimmed fastq files, mapping files in sam and bam format, variant calls in vcf format, and a parameters file. Details are below:

sample.fastq: the original fastq files
sample.trimmed.fastq: trimmed fastq files 
sample.sam: trimmed fastq files that were mapped to the reference sequence
sample.bam: trimmed fastq files that were mapped to the reference sequence, but in bam format
sample.sorted.bam: trimmed fastq files that were mapped to the reference sequence, but in sorted bam format (necessary for SNP calling)
sample.vcf: unfiltered lofreq variant calls. These have NOT been filtered to account for minimum quality, coverage, or SNP frequency. 
sample.filtered.vcf: filtered variant calls, filtered with the parameters you specified in the config file. 
sample.params: a parameters file containing a summary of the commands you specified




#########################################################################################
Configuring new genomes in snpEff
#########################################################################################

Unfortunately, putting together coding regions for each genome is a little tedious, especially if you are mapping each sample to its own reference and need to annotate a new set of coordinates for each sample. There is no very good way that I have thought of to automate this, so we are unfortunately left with these fairly specific and tedious instructions for building new genome databases with snpEff for each new genome you are using. 

snpEff is a program that annotates vcf files. It functions by using genome databases that are publicly available to determine coding region coordinate. Unfortunately, this is strongly skewed towards mammalian genomes and there are only a handful of viral sequences that are already prepared. It is therefore very likely that you will need to put together your own genomic coordinates and configure new genomes for each sample. Here is how to do it. 

Detailed instructions are provided on snpEff's website: http://snpeff.sourceforge.net/SnpEff_manual.html#databases

However, there are a few important quirks that are not specified in this document, so here are my notes for how to do this. 

1. Configure a new genome by adding to the the snpEff.config file. 
	a. Open this file, which is in the snpEff folder in any text editor. 
	b. Scroll to the very bottom of the file
	c. add your new reference sequence to this file using the following format: 

reference_sequence_name.genome : nickname for that reference sequence

so, for example, when I added the CA04 HA entry from GenBank GQ117044, my reference_sequence_name was the full name of the fasta file I was using, and my nickname was A_California_04_2009_HA. The .genome is absolutely necessary. It will NOT work if the .genome is not present in the name. 

	d. save the altered configuration file


2. Build the database using a gtf file. This is recommended. 
	a. cd into snpEff/data
	b. make a new directory with the name of the reference sequence. The name of the 			folder should be exactly identical to reference_sequence_name. It should NOT include the 	.genome. 
	c. place your genome fasta file and the associated gtf file into the reference_sequence_name folder. The reference sequence fasta file MUST be named sequences.fa. The gtf MUST be named genes.gtf. 
	d. run the following command: java -jar snpEff.jar build -gtf22 -v reference_sequence_name
	
	* Note: when you run the build command, reference_sequence_name does NOT include the .genome. If you have multiple gene segments, that is fine. They just need to all have the same base names in the gtf and sequences.fa file. (so CA04_HA and CA04_NA is fine). 

	* Once these steps are complete, your genome should be built in a zipped file called snpEffectPredictor.bin
	
3. Annotate your vcfs
Now, simply run the command: java -jar snpEff.jar reference_sequence_name input.vcf > output.vcf



#########################################################################################
Questions and contact information #########################################################################################

If you have questions, contact Louise Moncla at lhmoncla@gmail.com
