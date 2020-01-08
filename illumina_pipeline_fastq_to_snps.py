#!/usr/bin/env python

import sys, subprocess, glob, os, shutil, re, importlib, datetime
from subprocess import call
from datetime import datetime
from time import gmtime, strftime
time_stamp = strftime("%Y-%m-%d", gmtime())

config_filename = sys.argv[1]
config = __import__(config_filename)						# imports config_filename as a module
cfg = config.configuration()								# make cfg an instance of class configuration accessed from the config module
date_time = str(datetime.now())
reference_sequence_name = cfg.reference_sequence.split("/")[-1]

file_list = []
sample_list = []
sample_dict = {}
reference_list = []


# make a list of all the fastq files in the directory and store their sample name as the fastq file name without the R1_001 or R2_001.
for file in glob.glob("*.fastq"):						# glob will find instances of all files ending in .fastq as shell would
	file_list.append(file)

for file in glob.glob("*.fastq.gz"):						# glob will find instances of all files ending in .fastq.gz as shell would
	file_list.append(file)

for file in file_list:									# find all reads that are pairs
	if "R1" in file:
		samplename = file.replace("_R1", "")		# rename R1 reads to take out the R1_001.fastq
		samplename = samplename.replace(".fastq.gz", "")
		samplename = samplename.replace(".fastq", "")

	elif "R2" in file:
		samplename = file.replace("_R2", "")		# rename R2 reads to take out the R2_001.fastq		CHANGE THIS BACK 
		samplename = samplename.replace(".fastq.gz", "")
		samplename = samplename.replace(".fastq", "")

	else:
		samplename=file									# if there are only nonpaired reads, just keep the same name
		samplename = samplename.replace(".fastq.gz", "")
		samplename = samplename.replace(".fastq", "")	

	sample_list.append(samplename)

	# set up sample_dictionary with the sample name as keys and filenames (read names) as values
	if samplename in sample_dict:
		sample_dict[samplename].append(file)

	else:
		sample_dict[samplename] = [file]

sample_list = list(set(sample_list))					# keep only unique list entries

def write_log_file():
	for s in sample_dict:
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "w") as log_file:
			log_file.write("")


# combine R1 and R2 fastq files and move into directories based on their sample name; if using different references for each mapping, this will also move each reference sequence into the sample folder
def combine_fastqfiles():

	# make folders for each sample name and combine the forward and reverse reads into each folder
	for s in sample_dict:
		call("mkdir {s}".format(s=s), shell=True)

		for file in file_list:
			if s in file.replace("_R2", "") or s in file.replace("_R1", ""):
				call("cp {file} {s}/{file}".format(s=s, file=file), shell=True)


# remove human reads by mapping to a human reference genome and then collecting the unmapped reads 
def remove_human_reads(sample_list):
	for s in sample_dict:
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		call("mkdir {s}/raw-fastqs-with-human-reads".format(s=s),shell=True)
		
		with open(log_filename, "w") as log_file:

			for fastq in sample_dict[s]:
				with_human_fastq_name = fastq.replace(".fastq",".with-human-reads.fastq")
				human_removed_fastq_name = fastq.replace(".gz","")

				print("removing human reads from %s" % fastq)
				# move raw fastq files into new directory and rename; remove human reads and place new, human-filtered reads into same directory
				call("mv {s}/{fastq} {s}/raw-fastqs-with-human-reads/{with_human_fastq_name}".format(fastq=fastq,s=s,with_human_fastq_name=with_human_fastq_name), shell=True)
				call("bowtie2 -x {human_reference_sequence} -U {s}/raw-fastqs-with-human-reads/{with_human_fastq_name} -S {s}/raw-fastqs-with-human-reads/{s}.human.sam --un {s}/{human_removed_fastq_name} --local".format(s=s,fastq=fastq, human_reference_sequence=cfg.human_reference_sequence,with_human_fastq_name=with_human_fastq_name,human_removed_fastq_name=human_removed_fastq_name), shell=True, stderr=log_file)
				call("rm {s}/raw-fastqs-with-human-reads/{s}.human.sam".format(s=s), shell=True)
				
				if ".gz" in fastq: 
					call("gzip {s}/{human_removed_fastq_name}".format(human_removed_fastq_name=human_removed_fastq_name, s=s), shell=True)
					
												
# perform trimming
def trim(sample_list):

	# set the ending for trimmed file names depending on whether a paired or unpaired trim was performed
	if cfg.paired_trim == True:
		trimmed_forward_file_ending = "_1P.fastq"
		trimmed_reverse_file_ending = "_2P.fastq"
	elif cfg.paired_trim == False:
		trimmed_forward_file_ending = ".trimmed.fastq"
		trimmed_reverse_file_ending = ".trimmed.fastq"

	for s in sample_dict:
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "w") as log_file:
		
			for value in sample_dict[s]:
				new_name = value.replace(".fastq.gz", "")
				new_name = new_name.replace(".fastq","")
				if cfg.paired_trim == False and cfg.remove_adapters == False:
					call("trimmomatic SE {s}/{value} {s}/{new_name}.trimmed.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, new_name=new_name, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True, stderr=log_file)
				if cfg.paired_trim == False and cfg.remove_adapters == True:
					call("trimmomatic SE {s}/{value} {s}/{new_name}.trimmed.fastq ILLUMINACLIP:{adapters_fasta}:1:30:10 SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, new_name=new_name, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore, adapters_fasta=cfg.adapters_fasta), shell=True, stderr=log_file)

				elif cfg.paired_trim == True and cfg.remove_adapters == False:
					call("trimmomatic PE {s}/{value} {s}/{value} -baseout {s}/{s}.trimmed.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True, stderr=log_file)

				elif cfg.paired_trim == True and cfg.remove_adapters == True:
					call("trimmomatic PE {s}/{value} {s}/{value} -baseout {s}/{s}.trimmed.fastq ILLUMINACLIP:{adapters_fasta}:1:30:10 SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, new_name=new_name, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore, adapters_fasta=cfg.adapters_fasta), shell=True, stderr=log_file)


# perform mapping
def map(sample_list):

	if cfg.use_different_reference_for_each_sample == False:
		print("generating reference sequence with bowtie2")
		call("bowtie2-build {reference_sequence} {reference_sequence_name}".format(reference_sequence=cfg.reference_sequence, reference_sequence_name=reference_sequence_name), shell=True)
		
		for s in sample_dict:
			log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
			with open(log_filename, "a") as log_file:

				# set the trimmed fastq files to use for mapping
				trimmed_reads = []
				for f in os.listdir(s):
					if ".trimmed" in f and f.endswith(".fastq"):
						trimmed = s + "/" + f
						trimmed_reads.append(trimmed)
				fastqs_to_map = ",".join(trimmed_reads)

				print("now mapping %s" % s)
				
				if  cfg.paired_trim == False: 
					call("bowtie2 -x {reference_sequence} -U {fastqs_to_map} -S {s}/{s}.sam --local".format(s=s,fastqs_to_map=fastqs_to_map, reference_sequence=cfg.reference_sequence), shell=True, stderr=log_file)
				
				elif  cfg.paired_trim == True:
					# define input fastqs for paired trimming
					paired1 = []
					paired2 = []
					unpaired1 = []
					unpaired2 = []
					
					for f in trimmed_reads:
						
						if f.endswith("_1P.fastq"):
							paired1.append(f)
						elif f.endswith("_2P.fastq"):
							paired2.append(f)
						elif f.endswith("_1U.fastq"):
							unpaired1.append(f)
						elif f.endswith("_2U.fastq"):
							unpaired2.append(f)
					
					paired1 = ",".join(paired1)
					paired2 = ",".join(paired2)
					unpaired = ",".join(unpaired1 + unpaired2)
						
					call("bowtie2 -x {reference_sequence} -1 {paired1} -2 {paired2} -U {unpaired} -S {s}/{s}.sam --local".format(s=s,paired1=paired1,paired2=paired2,unpaired=unpaired, reference_sequence=cfg.reference_sequence), shell=True, stderr=log_file)

				call("samtools view -bS {s}/{s}.sam|samtools sort|samtools view -h > {s}/{s}.sorted.sam".format(s=s), shell=True, stderr=log_file)
				call("java -Xmx4g -jar /usr/local/bin/BAMStats-1.25/BAMStats-1.25.jar -i {s}/{s}.sorted.sam > {s}/{s}.coverage_stats.txt".format(s=s), shell=True, stderr=log_file)

	elif cfg.use_different_reference_for_each_sample == True:
		for s in sample_dict:
			
			log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
			with open(log_filename, "a") as log_file:

				# set the trimmed fastq files to use for mapping
				trimmed_reads = []
				for f in os.listdir(s):
					if ".trimmed" in f and f.endswith(".fastq"):
						trimmed = s + "/" + f
						trimmed_reads.append(trimmed)
				fastqs_to_map = ",".join(trimmed_reads)

				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						reference_name = file
				
				print("now mapping %s" % s)
				
				call("rm {s}/bowtie_reference_files".format(s=s), shell=True)
				call("bowtie2-build {s}/{reference_name} {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True, stderr=log_file)
				call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
				
				if  cfg.paired_trim == False: 
					call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_to_map} -S {s}/{s}.sam --local".format(s=s, fastqs_to_map=fastqs_to_map, reference_name=reference_name), shell=True, stderr=log_file)

				elif  cfg.paired_trim == True:
					# define input fastqs for paired trimming
					paired1 = []
					paired2 = []
					unpaired1 = []
					unpaired2 = []
					
					for f in trimmed_reads:
						
						if f.endswith("_1P.fastq"):
							paired1.append(f)
						elif f.endswith("_2P.fastq"):
							paired2.append(f)
						elif f.endswith("_1U.fastq"):
							unpaired1.append(f)
						elif f.endswith("_2U.fastq"):
							unpaired2.append(f)
					
					paired1 = ",".join(paired1)
					paired2 = ",".join(paired2)
					unpaired = ",".join(unpaired1 + unpaired2)

					call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -1 {paired1} -2 {paired2} -U {unpaired} -S {s}/{s}.sam --local".format(s=s,paired1=paired1,paired2=paired2,unpaired=unpaired, reference_name=reference_name), shell=True, stderr=log_file)

				call("samtools view -bS {s}/{s}.sam|samtools sort|samtools view -h > {s}/{s}.sorted.sam".format(s=s), shell=True, stderr=log_file)
				call("bamstats -i {s}/{s}.sorted.sam > {s}/{s}.coverage_stats.txt".format(s=s), shell=True, stderr=log_file)
				

# perform duplicate read removal with picard
def remove_duplicate_reads():
	for s in sample_dict:
		print("removing duplicate reads for %s using picard" % s)
		
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "a") as log_file:
		
			#convert sam to bam, then sort bam file
			call("samtools view -bS {s}/{s}.sam|samtools sort|samtools view -h > {s}/{s}.sorted.sam".format(s=s), shell=True, stderr=log_file)
		
			# run picard
			call("mkdir {s}/coverage_norm_and_duplicate_read_removal".format(s=s), shell=True, stderr=log_file)
			call("picard MarkDuplicates I={s}/{s}.sorted.sam O={s}/coverage_norm_and_duplicate_read_removal/{s}.nodups.sam REMOVE_DUPLICATES=true M=file.params.txt".format(s=s), shell=True, stderr=log_file)
		
			# clean up names and remove the .sorted.sam file
			call("rm {s}/{s}.sorted.sam; rm {s}/{s}.sam.sorted.bam".format(s=s), shell=True, stderr=log_file)



def normalize_coverage():
	for s in sample_dict:
		print("normalizing coverage for %s using bbnorm" % s)
		
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "a") as log_file:

			call("mkdir {s}/coverage_norm_and_duplicate_read_removal".format(s=s), shell=True)
		
			deduped = ""
			original = ""
		
			# set which sam file to use for normalization, based on whether duplicate read removal is turned on or off
			for f in os.listdir(s):
				if f == s + ".sam":
					original = f
			for f in os.listdir(s+"/coverage_norm_and_duplicate_read_removal"):
				if f.endswith(".nodups.sam") == True:
					deduped = "coverage_norm_and_duplicate_read_removal/" + f

			if cfg.remove_duplicate_reads == True:
				sam_file = deduped
				output_sam_name = "coverage_norm_and_duplicate_read_removal/" + s + ".nodups.normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
				normalized_fastq_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".nodups.normalized.fastq"
			elif cfg.remove_duplicate_reads == False:
				sam_file = original
				output_sam_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
				normalized_fastq_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".normalized.fastq"

			call("reformat.sh in={s}/{sam_file} out={s}/{s}.reextracted_from_sam.fastq".format(s=s, sam_file=sam_file), shell=True, stderr=log_file)
	
			# run bbnorm on the combined file
			call("bbnorm.sh in={s}/{s}.reextracted_from_sam.fastq out={s}/{normalized_fastq_name} target={cov_depth}".format(s=s, normalized_fastq_name=normalized_fastq_name, cov_depth=cfg.coverage_normalization_depth), shell=True, stderr=log_file)
			call("rm {s}/{s}.reextracted_from_sam.fastq".format(s=s), shell=True, stderr=log_file)
		
		
			# remap with bowtie2 
			fastqs_to_map = s + "/" + normalized_fastq_name
			if cfg.use_different_reference_for_each_sample == False:
				call("bowtie2-build {reference_sequence} {reference_sequence_name}".format(reference_sequence=cfg.reference_sequence, reference_sequence_name=reference_sequence_name), shell=True, stderr=log_file)
				call("bowtie2 -x {reference_sequence} -U {fastqs_to_map} -S {s}/{output_sam_name} --local".format(s=s, output_sam_name=output_sam_name, fastqs_to_map=fastqs_to_map, reference_sequence=cfg.reference_sequence), shell=True, stderr=log_file)

			elif cfg.use_different_reference_for_each_sample == True:
				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						reference_name = file			
				call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_to_map} -S {s}/{output_sam_name} --local".format(s=s, output_sam_name=output_sam_name, fastqs_to_map=fastqs_to_map, reference_name=reference_name), shell=True, stderr=log_file)



def call_SNPs():
	for s in sample_dict:
		
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "a") as log_file:

			# define sam file based on whether duplicate read removal and coverage normalization are is turned on or off
			if cfg.remove_duplicate_reads == False and cfg.normalize_coverage == False:
				sam_file = s + ".sam"
			elif cfg.remove_duplicate_reads == True and cfg.normalize_coverage == False:
				sam_file = "coverage_norm_and_duplicate_read_removal/"+ s + ".nodups.sam"
			elif cfg.remove_duplicate_reads == True and cfg.normalize_coverage == True:
				sam_file = "coverage_norm_and_duplicate_read_removal/"+ s + ".nodups.normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
			elif cfg.remove_duplicate_reads == False and cfg.normalize_coverage == True:
				sam_file = "coverage_norm_and_duplicate_read_removal/"+ s + ".normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
		
			# define output vcf names
			vcf_name = sam_file.replace(".sam", ".varscan" + str(cfg.snp_frequency) + ".vcf")
			annotated_name = vcf_name.replace(".vcf", ".annotated.vcf")
				
			#convert sam to bam, then sort bam file
			call("samtools view -bS {s}/{sam_file}|samtools sort > {s}/{sam_file}.sorted.bam".format(s=s,sam_file=sam_file), shell=True, stderr=log_file)

			# assign names to snpEff_ref_name variable, which will be used for amino acid annotation with snpEff
			if cfg.use_different_reference_for_each_sample == False:
				snpEff_ref_name = reference_sequence_name

			elif cfg.use_different_reference_for_each_sample == True:
				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						snpEff_ref_name = file

			# fix names
			if ".fasta" in snpEff_ref_name:
				snpEff_ref_name = snpEff_ref_name.replace(".fasta", "")
			if ".fa" in snpEff_ref_name:
				snpEff_ref_name = snpEff_ref_name.replace(".fa", "")
			if "_H5_partial" in snpEff_ref_name:
				snpEff_ref_name = snpEff_ref_name.replace("_H5_partial","")
			if "full_genome" in snpEff_ref_name:
				snpEff_ref_name = snpEff_ref_name.replace("full_genome","")
			if "_mixed" in snpEff_ref_name:
				snpEff_ref_name = snpEff_ref_name.replace("_mixed","")


			# run SNP calling and amino acid change annotations if samples are mapped to their own references
			if cfg.use_different_reference_for_each_sample == True:
				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						reference_sequence = file

				print("now calling SNPs on %s using varscan" % s)
					
				call("samtools mpileup -A -d 1000000 {s}/{sam_file}.sorted.bam > {s}/{sam_file}.pileup -f {s}/{reference_sequence}".format(s=s, sam_file=sam_file, reference_sequence=reference_sequence), shell=True, stderr=log_file)
				call("varscan mpileup2snp {s}/{sam_file}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{vcf_name}".format(s=s, vcf_name=vcf_name, sam_file=sam_file, snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True, stderr=log_file)

				if cfg.annotate_aa_changes == True:
					call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name} > {s}/{annotated_name}".format(sam_file=sam_file, vcf_name=vcf_name, annotated_name=annotated_name, snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True, stderr=log_file)


			# run SNP calling and amino acid change annotations if samples are all mapped to the same reference
			elif cfg.use_different_reference_for_each_sample == False:
				print("now calling SNPs on %s using varscan" % s)
					
				call("samtools mpileup -A -d1000000 {s}/{sam_file}.sorted.bam > {s}/{sam_file}.pileup -f {reference_sequence}".format(s=s, sam_file=sam_file, reference_sequence=cfg.reference_sequence), shell=True, stderr=log_file)
				call("varscan mpileup2snp {s}/{sam_file}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{vcf_name}".format(s=s, vcf_name=vcf_name, sam_file=sam_file, snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True, stderr=log_file)

				if cfg.annotate_aa_changes == True:
					call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name} > {s}/{annotated_name}".format(snpEff_ref_name=snpEff_ref_name,s=s,vcf_name=vcf_name, annotated_name=annotated_name, snp_frequency=cfg.snp_frequency), shell=True, stderr=log_file)


# perform de novo assembly with Trinity
def de_novo_assembly():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		print("now performing de novo assembly on %s using Trinity" % s)
		
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "a") as log_file:

			# set the trimmed fastq files to use for de novo assembly
			trimmed_reads = []
			
			for f in os.listdir(s):
				if f.endswith(".trimmed.fastq") == True:
					trimmed = s + "/" + f
					trimmed_reads.append(trimmed)
			fastqs_to_assemble = ",".join(trimmed_reads)

			#call("mkdir {s}/trinity_output".format(s=s), shell=True)
			call("Trinity --seqType fq --single {fastqs_to_assemble} --output {s}/trinity_output".format(s=s,fastqs_to_assemble=fastqs_to_assemble), shell=True, stderr=log_file)

			# pipe the output to blastn to return the identity of the contigs
			call("blastn -db nt -query {s}/trinity_output/Trinity.fasta -out {s}/Trinity_BLAST_result.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length evalue stitle qcovhsp qstart qend sstart send' -remote".format(s=s), shell=True, stderr=log_file)


def de_novo_assemble_mapped_reads():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		print("now performing de novo assembly on mapped reads from %s using Trinity" % s)
		
		log_filename = s + "/"+ s + "_log_file_" + time_stamp + ".txt"
		with open(log_filename, "a") as log_file:

			# set the trimmed fastq files to use for de novo assembly
			trimmed_reads = []
			
			for f in os.listdir(s):
				if f.endswith(".trimmed.fastq") == True:
					trimmed = f
					trimmed_reads.append(trimmed)
			fastqs_to_assemble = ",".join(trimmed_reads)

			# cp the original and trimmed fastq files into the new directory, remove the white spaces and replace with _s
			print("removing white spaces from fastq files - a Trinity requirement")
			call("mkdir {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True, stderr=log_file)
		
			for i in range(0, len(trimmed_reads)):
				fastq = trimmed_reads[i]
				call("cp {s}/{fastq} {s}/trinity_de_novo_assembly_mapped_reads_only/{fastq}.nospaces.fastq".format(s=s, fastq=fastq), shell=True, stderr=log_file)
			call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.fastq; do sed -i '' $'s/\ /\_/g' $f; done".format(s=s), shell=True, stderr=log_file)

			# redo bowtie mapping and put the white spaces back into the sam file
			if cfg.use_different_reference_for_each_sample == True:
				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						reference_name = file
						reference_sequence = file
		
			elif cfg.use_different_reference_for_each_sample == False:
				reference_name = reference_sequence_name
				reference_sequence = cfg.reference_sequence
		
			trimmed_reads_no_spaces = []
			for i in trimmed_reads: 
				no_spaces_reads = s + "/trinity_de_novo_assembly_mapped_reads_only/" + i + ".nospaces.fastq"
				trimmed_reads_no_spaces.append(no_spaces_reads)
			fastqs_no_spaces = ",".join(trimmed_reads_no_spaces)
		
			print("rebuild bowtie reference sequence")
			call("rm {s}/bowtie_reference_files; rm {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True, stderr=log_file)
			call("mkdir {s}/bowtie_reference_files; cp {reference_sequence} {s}/bowtie_reference_files/{reference_name}".format(s=s,reference_sequence=reference_sequence, reference_name=reference_name), shell=True, stderr=log_file)
			call("bowtie2-build {s}/bowtie_reference_files/{reference_name} {s}/bowtie_reference_files/{reference_name}".format(s=s, reference_name=reference_name), shell=True, stderr=log_file)
			#call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
			call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_no_spaces} -S {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam --local".format(s=s, fastqs_no_spaces=fastqs_no_spaces,reference_name=reference_name), shell=True, stderr=log_file)

			# put back the spaces so that trinity can use it
			print("add spaces back into fastq files for bowtie")
			call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_1\:N/\ 1\:N/g' $f; done".format(s=s), shell=True, stderr=log_file)
			call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_2\:N/\ 2\:N/g' $f; done".format(s=s), shell=True, stderr=log_file)

			# split sam file into forward and reverse and then extract fastqs
			print("split into forward and reverse reads and pull out fastq files")
			call("splitsam.sh {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.unmapped.sam".format(s=s), shell=True, stderr=log_file)
			call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.fastq ow=t".format(s=s), shell=True, stderr=log_file)
			call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.fastq ow=t".format(s=s), shell=True, stderr=log_file)

			# perform de novo assembly
			print("use Trinity for de novo assembly of %s" % s)
			call("Trinity --seqType fq --single {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.fastq,{s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.fastq --max_memory 3G --output {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True, stderr=log_file)

			# pipe the output to blastn to return the identity of the contigs
			print("now BLAST searching contigs from Trinity assembly of mapped reads from %s" % s)
			call("blastn -db nt -query {s}/trinity_de_novo_assembly_mapped_reads_only/Trinity.fasta -out {s}/trinity_de_novo_assembly_mapped_reads_only/Trinity_BLAST_result.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length evalue stitle qcovhsp qstart qend sstart send' -remote".format(s=s), shell=True, stderr=log_file)



# create an output parameter file
def create_parameter_file():
	for s in sample_dict:
		outfile_name = s + ".params.txt"

		with open(outfile_name, "w") as outfile:
			outfile.write("")
		with open(outfile_name, "w") as outfile:
			outfile.write("TRIMMING: \nMinimum read length after trimming: %s \nTrim window size: %s \nTrim bases to quality threshold of: Q%s \n\n" % (cfg.minlength, cfg.window_size, cfg.trim_qscore))

			if cfg.use_different_reference_for_each_sample == True:
				for file in os.listdir(s):
					if file.endswith(".fasta") or file.endswith(".fa"):
						reference_name = file
				outfile.write("REFERENCE SEQUENCE:\n%s\n\n" % file)

			elif cfg.use_different_reference_for_each_sample == False:
				outfile.write("REFERENCE SEQUENCE:\n%s\n\n" % cfg.reference_sequence)

			outfile.write("SNP CALLING AND FILTERING: \nMinimum required coverage: %s\nMinimum quality score for the variant base: %s\nMinimum SNP frequency: %s" % (cfg.min_cov, cfg.snp_qual_threshold, cfg.snp_frequency))
		outfile.close

		call("mv {s}.params.txt {s}/{s}.params.txt".format(s=s), shell=True)



# RUN THE ANALYSES
#write_log_file()

if cfg.combine_fastqs == True: 
	combine_fastqfiles()
	
if cfg.remove_human_reads == True:
	remove_human_reads(sample_list)

if cfg.trim == True:
	trim(sample_list)

if cfg.map == True:
	map(sample_list)

if cfg.remove_duplicate_reads == True:
	remove_duplicate_reads()

if cfg.normalize_coverage == True:
	normalize_coverage()

if cfg.call_snps == True:
	call_SNPs()

if cfg.de_novo_assemble_mapped_reads == True:
	de_novo_assemble_mapped_reads()

create_parameter_file()

if cfg.de_novo_assembly == True:
	de_novo_assembly()

