# Bioinformatics pipeline for Illumina sequence data analysis

This repository contains an Illumina sequence data analysis pipeline that I developed during PhD in Thomas Friedrich's lab. This pipeline handles going from raw fastq files that came off the MiSeq to generating variant calls. 

* Illumina pipeline
* Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

## How do I get set up? ##

This pipeline is written in Python 3.7, but should be compatible with python 2 as well. To use, download the python file and config file, fill out the config file, and add to your path. I have included a .yml file in this repository which will create an environment that contains all of the necessary dependencies installed. 

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) for your machine. 

2. Clone this repo: `git clone https://github.com/lmoncla/illumina_pipeline.git`

3. Navigate into the `illumina-pipeline` directory and create the conda environment:
` cd illumina-pipeline/`

`conda env create -f illumina-pipeline.yml`

4. Activate the environment: `conda activate illumina-pipeline`

5. `unzip example-data.zip` or double-click to unzip it. 

6. `cd example-data`

7. Fill out config file, following instructions below, and run pipeline on example data. 

`python ../illumina_pipeline_fastq_to_snps.py config`

As specified when downloaded, the config file will trim, map, and call variants on the example files. These example files are influenza HA sequences from the strain A/California/04/2009, that were grown in cell culture. The reference sequence `CA04_HA_GQ117044.fa` provided in the `example-files` folder can be found on GenBank [here](https://www.ncbi.nlm.nih.gov/nuccore/GQ117044).

## Basic usage: ##

From within the directory containing your list of fastq files, run: 

`python ../illumina_pipeline_fastq_to_snps.py config`

Replacing `../illumina_pipeline_fastq_to_snps.py` with the path to that file in your directory system. 

### Input files
As this is currently written, the script will run on all of the fastq files that are in your current directory. These files should all be standard fastq files and need to end in .fastq or .fastq.gz. The program will automatically combine forward and reverse reads that derive from the same sample together into the same folder. This is determined based on those samples having the same name, but being differentiated by having the `R1` or `R2` designation. The pipeline has been written to allow you to either map all samples to the same reference sequence (as you would want to do for any sort of experimental evolution/infection study) or to map each sample to a unique reference (as you might want to do for clinical samples). The main difference that you need to worry for specifying between these 2 options is in the reference sequence section of the config file (see: Filling in the config file below).

### The config file
All parameters including values you alter in the various programs, which analyses you wish to perform, and your reference sequence are specified in the config file. The config.py file contains annotated notes about which fields to fill in and what they are, although you should always consult the manuals for the specific programs you are calling for more specific details about what they do. These arguments will then be passed to `illumina_pipeline_fastq_to_snps.py`, which will incorporate those arguments at the appropriate time in the pipeline. After running, the pipeline will output a log file for each sample, which will contain all of the arguments that you used for the analyses as well as any error messages. The intention for using a config file here is that you shouldn't have to alter the `illumina_pipeline_fastq_to_snps.py` file, and should be able to just alter the `config.py` file. Because you will have an log file, it will let you keep a record of what you ran so that you can directly edit the config.py file and repeat the analyses later if you want.


## Filling in the config file ##

### SECTION 1: SPECIFY WHICH TASKS YOU WANT TO DO HERE
You may elect to perform trimming, mapping, SNP calling, duplicate read removal, coverage normalization, and human read removal. To enable these analyses, simply type "True" (make sure to use a capital T) after the = each option. This pipeline does not automatically call consensus genomes. I usually really like to look at the raw, mapped data, so I usually import these into the trial version of Geneious and call a consensus genome from there (you can do this with the free, trial version). The config file is divided into 3 basic sections: basic tasks (trimming, mapping, calling and annotating SNPs and de novo assembly), cleaning tasks, which are performed after mapping on sam/bam files (coverage normalization and duplicate read removal). Specifics are below:

#### self.remove_human_reads = `True` or `False`
Sometimes I need to remove human reads from my fastq files as a first step to protect patient privacy. If you select True, this will map to a human reference genome with bowtie 2, and output all of the reads that did not map. It will retain the raw fastq files with human reads in a newly created folder. If you want to run this, you need to set the path to the human reference genome in the config file. Bowtie2 provides pre-compiled human reference genomes that you can download. I was not able to upload mine here, due to github's 100 MB size constraint. I use GR38, which can be downloaded [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). 

#### self.trim = `True` or `False`
Use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to trim the ends of your reads. This must be done for all raw fastq files. If set to True, fill in the parameters under the "SET TRIMMING PARAMETERS" section. Trimmomatic trims in sliding windows, for which you can specify the q-score cutoff and the window size. You can also elect to remove index adapter sequences. I've included the Nextera adapter sequences in this report as `Nextera_XT_adapter.fa`, but you will need to specify your own adapter sequence as a fasta file if you are using some other library prep method. 

#### self.map = `True` or `False`
Use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to map to a reference sequence. If set to True, fill in the "SPECIFY REFERENCE SEQUENCE" section.

#### self.call_snps = `True` or `False`
Use [Varscan](http://varscan.sourceforge.net/) to call SNPs after mapping to a reference sequence. You must map to a reference before calling SNPs, as the input file for SNP calling is output file for mapping. If set to True, fill in the "SET SNP CALLING PARAMETERS" section of the config file.

#### self.remove_duplicate_reads = `True` or `False`
Use [picard's MarkDuplicates tool](http://broadinstitute.github.io/picard/) to remove duplicate reads. Picard considers reads to be duplicates if they have the exact same 5' start site. In some comparisons with samblaster, picard and dedupe, picard performed the best. Picard improved data reproducibility for within-host influenza data, while samblaster and dedupe resulted in an increase in the number of low-frequency variants called, but did result in consistent calls. Therefore, I have elected to use picard for duplicate read removal in this pipeline. Duplicate read removal will occur after mapping, but before coverage normalization.

#### self.normalize_coverage = `True` or `False`
Use bbnorm.sh from the [bbmap software package](https://sourceforge.net/projects/bbmap/) to normalize coverage across the genome. This is performed after mapping. If duplicate read removal is also to be performed, normalization will occur after duplicate read removal. BBnorm.sh normalizes coverage using fastq files as input, so when self.normalize_coverage is set to `True`, reads will be extracted from the sam file, used as input for for bbnorm, and then remapped. The resulting output file will end in `.normalized.coverage.sam` where `coverage` will be the desired coverage depth that you set.

#### self.coverage_normalization_depth = `integer`
Set the desired depth of coverage you wish to achieve after coverage normalization to `integer`.


=======

### SECTION 2 : SET/ALTER PARAMETERS

#### SET TRIMMING PARAMETERS:
#### self.remove_adapters = `True` or `False`
This specifies whether you would like Trimmomatic to remove index adapter sequences or not. This depends on whether this has already been done by your sequencing core or not.

#### self.adapters_fasta = `path_to_index_adapters_fasta`
If you select `True` for `self.remove_adapters`, then you need to provide a path to the index adapters fasta file here.

#### self.paired_trim = `True` or `False`
Trimmomatic can be run in either paired mode or single end mode. I usually opt to run single end mode, but you should read about the differences between them and decide for yourself. 

#### self.window_size = `integer`
Trimmomatic performs read end trimming by sliding along the read and calculating a running quality score in sliding windows. The width of those windows is specified by `integer`. I usually set this to 5, which I believe is the default value. 

#### self.trim_qscore = `integer`
Phred-based quality score threshold to use during trimming. If you would like to use a Q30 threshold, you would specify 30. 30 is recommended.

#### self.minlength = `integer`
After read trimming has been performed, discard reads that are shorter than `integer` length. Value must be an integer. I would recommend 100.

=======

#### SPECIFY REFERENCE SEQEUNCE AND MAPPING QUALITY:
One important note here is that this pipeline is meant to run with a single reference sequence file. If you want to specify multiple gene segments, simply put all of them into the same fasta file. The fasta file must end in .fasta or .fa.

#### self.use_different_reference_for_each_sample = `True` or `False`
Specify True to map all of the samples to the same reference sequence or False to map each sample to it's own reference. If specifying False, then you need to put the fasta reference file into the same folder as the trimmed fastqs. The easiest way to do this is to run the pipeline and do only Trimming, which will combine the forward and reverse fastq files and make folders with their specific names. Then just move the fasta reference files into the appropriate folder.

#### self.reference_sequence = `path to reference sequence`
If you are mapping everything to the same reference sequence, then you have to specify the full path to the reference sequence you wish to use. The reference sequence should be in fasta format and can end in .fasta or .fa. Ex: `example-data/CA04_HA_GQ117044.fa`

#### self.mapping_quality_threshold = `30`
After mapping with bowtie, it is a good idea to remove reads from the sam file that have a low mapping quality score. Reads with mapping quality scores less than this value will be removed from the sam file. If you do not wish to use this option, set to 0. These values are specified as Phred scores.

=======
#### SET SNV CALLING PARAMETERS:
One important note here is that if you would like SNVs to be annotated as to whether they cause a coding region change, then you need to put together gtf files. This functionality is coming shortly, but is not yet available. 

#### self.min_coverage = `integer`
This will set the minimum coverage required at a base in order to perform variant calling at that base. Variants at positions with coverage less than integer will not be called.

#### self.snp_qual_threshold = `integer`
This will set the minimum quality score required at a base in order to perform variant calling at that base. Variants with quality scores lower than integer will not be called.

#### self.snp_frequency = `float`
Variants that are present at a frequency less than that set by decimal will not be called. decimal values should range from 0 to 1, with a value of 0.01 specifying that SNPs should be called at a 1% frequency cutoff.

=======

### Output
After this has been run, a folder will be made for each sample, which will contain the original fastq files, trimmed fastq files, mapping files in sam and bam format, variant calls in vcf format, and a log file. The log file contains the standard out for all of the shell commands being run within this pipeline and will contain any errors that stem from the programs being called. Details are below:

#### sample.fastq
the original fastq files. If human reads were removed, then these files only contain reads that did not map to the human reference sequence. 

#### sample.trimmed.fastq
trimmed fastq files

#### sample.sam
trimmed fastq files that were mapped to the reference sequence

#### sample.bam
trimmed fastq files that were mapped to the reference sequence, but in bam format

#### sample.sorted.bam
trimmed fastq files that were mapped to the reference sequence, but in sorted bam format (necessary for SNP calling)

#### sample.varscan.snps.`snp_frequency`.vcf
variants called by varscan according to the parameters you specified in the config file. Not yet annotated.

#### raw-fastqs-with-human-reads
This folder contains the original, raw fastq files that still contain human reads. This folder is only produced if `self.remove_human_reads = True`

#### coverage_norm_and_duplicate_read_removal
This folder will contain all the results of analyses performed on data where duplicate reads have been removed or the coverage has been normalized. This folder will contain the mapping files (sam/bam) and variant call files (vcfs). Files that contain ".nodups" are ones in which duplicate reads were removed. Files that contain ".normalized.`coverage_value`x" are ones in which coverage was normalized to an average coverage of `coverage_value`. Output files that contain ".nodups.normalized.`coverage_value`x" are ones that had both duplicates removed, and then coverage normalized afterwards. Doing both is usually unnecessary, as often removing duplicate reads results in much less and more even coverage.

#### log_file.txt
a log file containing the output that would have been printed to the screen from the programs that have been called. Will contain any errors, so this is a good place to start if you notice that you are missing a file that should have been generated, or there seem to be problems with your output.


## Dependencies/software used in this pipeline ##

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

#### varscan
http://varscan.sourceforge.net/
Varscan is another option for a variant caller. Allows variants to be called/filtered based on the same criteria as Lofreq (forward/reverse read balance, coverage, frequency, quality)

#### samtools
http://samtools.sourceforge.net/
Necessary for converting sam/bam files, sorting bam files, removing low mapping quality reads from sam files, and generating pileup files

#### bbmap  
https://sourceforge.net/projects/bbmap/
Used for extracting mapped fastq files from sam/bam files for input into Trinity de novo assembly for quality control. Also used for coverage normalization, with the tool bbnorm.sh.  



### Questions and comments ###

Louise Moncla
lhmoncla@gmail.com
