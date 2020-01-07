#!/usr/bin/env python


class configuration(object):
	def __init__(self):


####### SPECIFY WHICH TASKS YOU WANT TO DO HERE #########################################

		# for each of the below (trimming, mapping and calling SNPs), set to true if you want to do it, set to false if you do not. If you set something to false then you do not need to change any of the parameters for the associated analysis
		self.combine_fastqs = False  # calls combine fastq files into a folder and put results of analyses into that folder 
		self.remove_human_reads = False
		self.trim = False
		self.map = True
		self.call_snps = False
		self.annotate_aa_changes = False
		self.de_novo_assembly = False
		self.de_novo_assemble_mapped_reads = False

		# data cleaning tasks: coverage depth normalization with bbnorm and duplicate read removal with picard; these will be implemented upstream of variant calling, such that variant calling will use the de-duplicated or normalized sam/bam file
		self.remove_duplicate_reads = False
		
		# using bbnorm from the bbmap software package, normalize coverage across the sam or bam file to a set coverage depth, specified with self.coverage_normalization_depth = DEPTH
		self.normalize_coverage = False
		self.coverage_normalization_depth = 1000
	
		# popoolation analyses; set self.run_popoolation = True to run any of the subsequent popoolation analyses
		self.run_popoolation = False

		# run popoolation without subsampling
		self.calculate_genewise_pi = False
		self.calculate_genewise_piNpiS = False
		self.calculate_sliding_window_piNpiS = False

		# subsampled popoolation analyses
		self.perform_subsampling = False
		self.calculate_subsampled_pi = False
		self.calculate_subsampled_piNpiS = False
		self.calculate_subsampled_sliding_window_piNpiS = False



####### SET/ALTER PARAMETERS ############################################################

####### SET HUMAN REFERENCE GENOME TO USE TO REMOVE HUMAN READS ##########################
### Fill these out if self.remove_human_reads = True

		self.human_reference_sequence = "/Volumes/gradschool-and-postdoc-backups/post-doc/stored_files_too_big_for_laptop/Mumps/sequence_data/human-reference-genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"

####### SET TRIMMING PARAMETERS #########################################################
### Fill these out if self.trim = True

		# remove illumina adapters from sequence ends
		self.remove_adapters = True
		self.adapters_fasta = "/usr/local/bin/Trimmomatic-0.36/Nextera_XT_adapter.fa"

		# trim the reads as paired reads or as unpaired; set to True to run Trimmomatic in paired mode and False to run in unpaired mode
		self.paired_trim = False

		# after trimming, discard reads below this length
		self.minlength = 100

		# during trimming, slide along reads in windows of this size (in base pairs)
		self.window_size = 5

		# trim reads using a quality score threshold of this (for Q30, set to 30, etc...)
		self.trim_qscore = 30

###### SPECIFY REFERENCE SEQUENCE AND MAPPING QUALITY ########################################################
### Fill these out if self.map = True

		# here, put the full path for the reference sequence you wish to use for mapping. You can figure out the full path by dragging and dropping the file into the terminal, and then copying that file path into here
		self.reference_sequence = "/Volumes/gradschool-and-postdoc-backups/post-doc/stored_files_too_big_for_laptop/Mumps/sequence_data/mumps-seq-library-4/mapped-to-own-consensus/45/MuVs_Massachusetts.USA_42.16_1_G.fasta"

		# If, instead of mapping everything to the same reference you would like to map sample to a different reference, then specify True here. This will also require that the references you wish to use have been placed in the same folder as the trimmed fastq files.
		self.use_different_reference_for_each_sample = True

####### SET SNP CALLING PARAMETERS ######################################################
### Fill these out if self.call_snps = True

		# pick whether you would like to use LoFreq or Varscan to call variants
		self.use_lofreq = False
		self.use_varscan = True

		# minimum coverage, i.e., SNPs will not be called at sites that have coverage less than this value
		self.min_cov = 100

		# set base quality threshold, i.e., SNPs will not be called for bases that have a Qscore below this value
		self.snp_qual_threshold = 30

		# set SNP frequency cutoff (1% would be specified as 0.01). SNPs present below this frequency will not be reported
		self.snp_frequency = 0.01


####### SPECIFY POPOOLATION PARAMETERS ##################################################
### Fill these out if self.calculate_genewise_pi, self.calculate_genewise_piNpiS, or self.calculate_sliding_window_pi = True

		# set desired coverage level with subsample_level = desired depth
		self.subsample_level = 1000

		# this sets the --min-count parameter in popoolation, which is the minimum count of the minor allele for it to be counted
		self.min_count = 1

		# this sets the --min-coverage parameter in popoolation; sites with coverage less than this value will not be used for calculating pi
		self.min_coverage = 100

		# this sets the --max-coverage parameter in popoolation; sites with coverage exceeding this value will not be used for calculating pi
		self.max_coverage = 1000000

		# this sets the --min-quality parameter in popoolation; alleles with quality lower than this threshold will not be considered/used for diversity calculations
		self.min_quality = 0



### Fill these out if self.calculate_sliding_window_pi = True

		# for sliding window pi analyses, this specifies the --window-size and --step-size parameters; popoolation will slide along the gene in windows of size pi_step_size, calculating pi in windows of size pi_window_size
		self.pi_window_size = 9
		self.pi_step_size = 3





### NOTES ###
# To run the pipeline, navigate to the directory that contains your fastq files and run illumina_pipeline.py and the configuration file you wish to use.
# This configuration file will be used to inform the commands run by illumina_pipeline.py. The pipeline uses Trimmomatic for trimming, bowtie2 for mapping and lofreq for variant calling. I will assume that if you are using this pipeline that you have read the documentation for these programs and understand how they work and what impact these parameters have on the analysis. For specific information about these programs, please consult their online documentation atTrimmomatic: http://www.usadellab.org/cms/?page=trimmomatic, and bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml.
