#!/usr/bin/env python


class configuration(object):
	def __init__(self):
				

####### SPECIFY WHICH TASKS YOU WANT TO DO HERE #########################################
		
		# for each of the below (trimming, mapping and calling SNPs), set to true if you want to do it, set to false if you do not. If you set something to false then you do not need to change any of the parameters for the associated analysis
		self.trim = False
		self.map = True
		self.call_snps = True
		self.annotate_aa_changes = True
		self.de_novo_assembly = False
		self.de_novo_assemble_mapped_reads = False
		
		# this section is not yet ready for prime time
		self.calculate_genewise_pi = False
		self.calculate_genewise_piNpiS = False
		self.calculate_sliding_window_piNpiS = False
		
		
		
####### SET/ALTER PARAMETERS ############################################################
		
####### SET TRIMMING PARAMETERS #########################################################
### Fill these out if self.trim = True		
		
		# after trimming, discard reads below this length
		self.minlength = 100
		
		# during trimming, slide along reads in windows of this size (in base pairs)
		self.window_size = 5
		
		# trim reads using a quality score threshold of this (for Q30, set to 30, etc...)
		self.trim_qscore = 30
		
		
###### SPECIFY REFERENCE SEQUENCE ########################################################
### Fill these out if self.map = True		
		
		# here, put the full path for the reference sequence you wish to use for mapping. You can figure out the full path by dragging and dropping the file into the terminal, and then copying that file path into here
		self.reference_sequence = "/Volumes/LaCie/Users/lhmoncla/Documents/New_Illumina_Pipeline_Development/pipeline_testing/reference_based_mapping_pipeline/CA04_HA_GQ117044.fa"
		
		# specify what you would like your reference sequence to be named
		self.reference_sequence_name = "CA04_HA_GQ117044.fa"
		
		# If, instead of mapping everything to the same reference you would like to map sample to a different reference, then specify True here. This will also require that the references you wish to use have been placed in the same folder as the trimmed fastq files. 
		self.use_different_reference_for_each_sample = True		

####### SET SNP CALLING PARAMETERS ######################################################
### Fill these out if self.call_snps = True		
		
		# pick whether you would like to use LoFreq or Varscan to call variants
		self.use_lofreq = True
		self.use_varscan = True
				
		# minimum coverage, i.e., SNPs will not be called at sites that have coverage less than this value
		self.min_cov = 100
		
		# set base quality threshold, i.e., SNPs will not be called for bases that have a Qscore below this value
		self.snp_qual_threshold = 30
		
		# set SNP frequency cutoff (1% would be specified as 0.01). SNPs present below this frequency will not be reported
		self.snp_frequency = 0.01		
		

####### SPECIFY POPOOLATION PARAMETERS ##################################################
### Fill these out if self.calculate_genewise_pi, self.calculate_genewise_piNpiS, or self.calculate_sliding_window_pi = True		
		
		# Specify whether each sample should have its own gtf (set to True) or if the same one should be used for all (set to False)
		self.use_different_gtf_for_each_sample = True
		
		# if using a shared gtf for all files, specify the path to that gtf here; if using a different gtf for each sample, this isn't necessary, and the program will just assume the gtf is in the folder with the fastq files
		self.gtf_location = "/Directory/of/gtf"

		# set to True if you would like to subsample the pileup file. Then set desired coverage level with subsample_level = desired depth
		self.subsample = False
		self.subsample_level = 1000
		
		# this sets the --min-count parameter in popoolation, which is the minimum count of the minor allele for it to be counted
		self.min_count = 1
		
		# this sets the --min-coverage parameter in popoolation; sites with coverage less than this value will not be used for calculating pi
		self.min_coverage = 100 
		
		# this sets the --max-coverage parameter in popoolation; sites with coverage exceeding this value will not be used for calculating pi
		self.max_coverage = 1000000
		
		# this sets the --nonsyn-length-table in popoolation, which specifies which codon table you want to use (you can choose nsl-P1.txt or nsl-P2.txt)
		self.nonsyn_length_table = "nsl_p1.txt"
		
		

### Fill these out if self.calculate_sliding_window_pi = True		

		# for sliding window pi analyses, this specifies the --window-size and --step-size parameters; popoolation will slide along the gene in windows of size pi_step_size, calculating pi in windows of size pi_window_size
		self.pi_window_size = 9
		self.pi_step_size = 3
		




### NOTES ###
# To run the pipeline, navigate to the directory that contains your fastq files and run illumina_pipeline.py and the configuration file you wish to use. 
# This configuration file will be used to inform the commands run by illumina_pipeline.py. The pipeline uses Trimmomatic for trimming, bowtie2 for mapping and lofreq for variant calling. I will assume that if you are using this pipeline that you have read the documentation for these programs and understand how they work and what impact these parameters have on the analysis. For specific information about these programs, please consult their online documentation atTrimmomatic: http://www.usadellab.org/cms/?page=trimmomatic, and bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml. 