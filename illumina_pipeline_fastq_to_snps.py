#!/usr/bin/env python2

import sys, subprocess, glob, os, shutil, re, importlib
from subprocess import call

config_filename = sys.argv[1]
config = __import__(config_filename)						# imports config_filename as a module
cfg = config.configuration()								# make cfg an instance of class configuration accessed from the config module
reference_sequence_name = cfg.reference_sequence.split("/")[-1]

file_list = []
sample_list = []
sample_dict = {}
reference_list = []


# make a list of all the fastq files in the directory and store their sample name as the fastq file name without the R1_001 or R2_001.
for file in glob.glob("*.fastq"):						# glob will find instances of all files ending in .fastq as shell would
	file_list.append(file)

for file in file_list:									# find all reads that are pairs
	if "R1" in file:
		samplename = file.replace("_R1", "")		# rename R1 reads to take out the R1_001.fastq
		samplename = samplename.replace(".fastq", "")

	elif "R2" in file:
		samplename = file.replace("_R2", "")		# rename R2 reads to take out the R2_001.fastq		CHANGE THIS BACK 
		samplename = samplename.replace(".fastq", "")

	else:
		samplename=file									# if there are only nonpaired reads, just keep the same name
		samplename = samplename.replace(".fastq", "")	

	sample_list.append(samplename)

	# set up sample_dictionary with the sample name as keys and filenames (read names) as values
	if samplename in sample_dict:
		sample_dict[samplename].append(file)

	else:
		sample_dict[samplename] = [file]

sample_list = list(set(sample_list))					# keep only unique list entries


# combine R1 and R2 fastq files and move into directories based on their sample name; if using different references for each mapping, this will also move each reference sequence into the sample folder
def combine_fastqfiles():

	# make folders for each sample name and combine the forward and reverse reads into each folder
	for s in sample_dict:
		call("mkdir {s}".format(s=s), shell=True)

		for file in file_list:
			if s in file:
				call("cp {file} {s}/{file}".format(s=s, file=file), shell=True)



# set the ending for trimmed file names depending on whether a paired or unpaired trim was performed
if cfg.paired_trim == True:
	trimmed_forward_file_ending = "_1P.fastq"
	trimmed_reverse_file_ending = "_2P.fastq"
elif cfg.paired_trim == False:
	trimmed_forward_file_ending = ".trimmed.fastq"
	trimmed_reverse_file_ending = ".trimmed.fastq"



# perform trimming
def trim(sample_list):
	for s in sample_dict:
		for value in sample_dict[s]:
			if cfg.paired_trim == False and cfg.remove_adapters == False:
				call("java -jar /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE {s}/{value} {s}/{value}.trimmed.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True)
			if cfg.paired_trim == False and cfg.remove_adapters == True:
				call("java -jar /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE {s}/{value} {s}/{value}.trimmed.fastq ILLUMINACLIP:{adapters_fasta}:1:30:10 SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore, adapters_fasta=cfg.adapters_fasta), shell=True)

			elif cfg.paired_trim == True:
				call("java -jar /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE {s}/{value} {s}/{value} -baseout {s}/{s}.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True)


# perform mapping
def map(sample_list):

	if cfg.use_different_reference_for_each_sample == False:
		call("bowtie2-build {reference_sequence} {reference_sequence_name}".format(reference_sequence=cfg.reference_sequence, reference_sequence_name=reference_sequence_name), shell=True)
		for s in sample_dict:
			
			# set the trimmed fastq files to use for mapping
			trimmed_reads = []
			for f in os.listdir(s):
				if f.endswith(".trimmed.fastq") == True:
					trimmed = s + "/" + f
					trimmed_reads.append(trimmed)
			fastqs_to_map = ",".join(trimmed_reads)

			call("bowtie2 -x {reference_sequence} -U {fastqs_to_map} -S {s}/{s}.sam --local".format(s=s,fastqs_to_map=fastqs_to_map, reference_sequence=cfg.reference_sequence), shell=True)
			#call("samtools view -h -q {mapping_quality} {s}/mapped.sam > {s}/{s}.sam; rm {s}/mapped.sam".format(s=s, mapping_quality=cfg.mapping_quality_threshold), shell=True)

	elif cfg.use_different_reference_for_each_sample == True:
		for s in sample_dict:
			
			# set the trimmed fastq files to use for mapping
			trimmed_reads = []
			
			for f in os.listdir(s):
				if f.endswith(".trimmed.fastq") == True:
					trimmed = s + "/" + f
					trimmed_reads.append(trimmed)
			fastqs_to_map = ",".join(trimmed_reads)

			for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					reference_name = file
			call("rm {s}/bowtie_reference_files".format(s=s), shell=True)
			call("bowtie2-build {s}/{reference_name} {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True)
			call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
			
			call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_to_map} -S {s}/{s}.sam --local".format(s=s, fastqs_to_map=fastqs_to_map, reference_name=reference_name), shell=True)
			#call("samtools view -h -q {mapping_quality} {s}/mapped.sam > {s}/{s}.sam; rm {s}/mapped.sam".format(s=s, mapping_quality=cfg.mapping_quality_threshold), shell=True)


# perform duplicate read removal with picard
def remove_duplicate_reads():
	for s in sample_dict:
		#convert sam to bam, then sort bam file
		call("samtools view -bS {s}/{s}.sam|samtools sort|samtools view -h > {s}/{s}.sorted.sam".format(s=s), shell=True)
		
		# run picard
		call("mkdir coverage_norm_and_duplicate_read_removal".format(), shell=True)
		call("java -jar /usr/local/bin/picard.jar MarkDuplicates I={s}/{s}.sorted.sam O={s}/coverage_norm_and_duplicate_read_removal/{s}.nodups.sam REMOVE_DUPLICATES=true M=file.params.txt".format(s=s), shell=True)
		
		# clean up names and remove the .sorted.sam file
		call("rm {s}/{s}.sorted.sam; rm {s}/{s}.sam.sorted.bam".format(s=s), shell=True)



def normalize_coverage():
	for s in sample_dict:
		print "normalizing coverage for %s using bbnorm" % s
		call("mkdir coverage_norm_and_duplicate_read_removal".format(), shell=True)
		
		deduped = ""
		original = ""
		
		# set which sam file to use for normalization, based on whether duplicate read removal is turned on or off
		for f in os.listdir(s):
			if f.endswith(".nodups.sam") == True:
				deduped = f
			elif f == s + ".sam":
				original = f
		if cfg.remove_duplicate_reads == True:
			sam_file = deduped
			output_sam_name = "coverage_norm_and_duplicate_read_removal/" + s + ".nodups.normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
			normalized_fastq_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".nodups.normalized.fastq"
		elif cfg.remove_duplicate_reads == False:
			sam_file = original
			output_sam_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
			normalized_fastq_name = "coverage_norm_and_duplicate_read_removal/"+ s + ".normalized.fastq"

		call("reformat.sh in={s}/{sam_file} out={s}/{s}.reextracted_from_sam.fastq".format(s=s, sam_file=sam_file), shell=True)	
	
		# run bbnorm on the combined file
		call("bbnorm.sh in={s}/{s}.reextracted_from_sam.fastq out={s}/{normalized_fastq_name} target={cov_depth}".format(s=s, normalized_fastq_name=normalized_fastq_name, cov_depth=cfg.coverage_normalization_depth), shell=True)
		call("rm {s}/{s}.reextracted_from_sam.fastq".format(s=s), shell=True)
		
		
	# remap with bowtie2 
		fastqs_to_map = s + "/" + normalized_fastq_name
		if cfg.use_different_reference_for_each_sample == False:
			call("bowtie2-build {reference_sequence} {reference_sequence_name}".format(reference_sequence=cfg.reference_sequence, reference_sequence_name=reference_sequence_name), shell=True)
			call("bowtie2 -x {reference_sequence} -U {fastqs_to_map} -S {s}/{output_sam_name} --local".format(s=s, output_sam_name=output_sam_name, fastqs_to_map=fastqs_to_map, reference_sequence=cfg.reference_sequence), shell=True)
			#call("samtools view -h -q {mapping_quality} {s}/mapped.sam > {s}/{output_sam_name}; rm {s}/mapped.sam".format(s=s, output_sam_name=output_sam_name, mapping_quality=cfg.mapping_quality_threshold), shell=True)

		elif cfg.use_different_reference_for_each_sample == True:
			for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					reference_name = file
				#call("rm {s}/bowtie_reference_files".format(s=s), shell=True)
				#call("bowtie2-build {s}/{reference_name} {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True)
				#call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
			
			call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_to_map} -S {s}/{output_sam_name} --local".format(s=s, output_sam_name=output_sam_name, fastqs_to_map=fastqs_to_map, reference_name=reference_name), shell=True)
			#call("samtools view -h -q {mapping_quality} {s}/mapped.sam > {s}/{output_sam_name}; rm {s}/mapped.sam".format(s=s, output_sam_name=output_sam_name, mapping_quality=cfg.mapping_quality_threshold), shell=True)



def call_SNPs():
	for s in sample_dict:
		
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
		if cfg.use_lofreq == True:
			vcf_name = sam_file.replace(".sam", ".lofreq" + str(cfg.snp_frequency) + ".vcf")
		elif cfg.use_varscan == True:
			vcf_name = sam_file.replace(".sam", ".varscan" + str(cfg.snp_frequency) + ".vcf")
		annotated_name = vcf_name.replace(".vcf", ".annotated.vcf")
				
		#convert sam to bam, then sort bam file
		call("samtools view -bS {s}/{sam_file}|samtools sort > {s}/{sam_file}.sorted.bam".format(s=s,sam_file=sam_file), shell=True)
		call("rm {s}/{sam_file}.lofreq.{freq}.vcf; rm {s}/{sam_file}.filtered.{freq}.vcf".format(s=s,sam_file=sam_file,freq=cfg.snp_frequency), shell=True)

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

			if cfg.use_lofreq == True:
				print "now calling SNPs on %s using lofreq" % s
				call("lofreq call -f {s}/{reference_sequence} -o {s}/{vcf_name} {s}/{sam_file}.sorted.bam".format(s=s, sam_file=sam_file, vcf_name=vcf_name, reference_sequence=reference_sequence), shell=True)
				call("lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {vcf_name} -o {vcf_name}.filtered.vcf".format(s=s, vcf_name=vcf_name, min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold, snp_frequency=cfg.snp_frequency), shell=True)

				if cfg.annotate_aa_changes == True:
					call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name}.filtered.vcf > {s}/{annotated_name}".format(snpEff_ref_name=snpEff_ref_name,s=s,vcf_name=vcf_name, annotated_name=annotated_name, snp_frequency=cfg.snp_frequency), shell=True)

			if cfg.use_varscan == True:
				print "now calling SNPs on %s using varscan" % s
				call("samtools mpileup -d 1000000 {s}/{sam_file}.sorted.bam > {s}/{sam_file}.pileup -f {s}/{reference_sequence}".format(s=s, sam_file=sam_file, reference_sequence=reference_sequence), shell=True)
				call("java -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp {s}/{sam_file}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{vcf_name}".format(s=s, vcf_name=vcf_name, sam_file=sam_file, snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True)

				if cfg.annotate_aa_changes == True:
					call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name} > {s}/{annotated_name}".format(sam_file=sam_file, vcf_name=vcf_name, annotated_name=annotated_name, snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True)


		# run SNP calling and amino acid change annotations if samples are all mapped to the same reference
		elif cfg.use_different_reference_for_each_sample == False:
			# call variants with lofreq and filter them
			if cfg.use_lofreq == True:
				print "now calling SNPs on %s using lofreq" % s
				call("lofreq call -f {reference_sequence} -o {s}/{vcf_name} {sam_file}.sorted.bam".format(s=s, vcf_name=vcf_name, sam_file=sam_file, snp_frequency=cfg.snp_frequency, reference_sequence=cfg.reference_sequence), shell=True)
				call("lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {s}/{vcf_name} -o {s}/{vcf_name}.filtered.vcf".format(s=s, vcf_name=vcf_name, min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold, snp_frequency=cfg.snp_frequency), shell=True)

				if cfg.annotate_aa_changes == True:
					call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name} > {s}/{annotated_name}".format(snpEff_ref_name=snpEff_ref_name,s=s,vcf_name=vcf_name, annotated_name=annotated_name, snp_frequency=cfg.snp_frequency), shell=True)

			if cfg.use_varscan == True:
				print "now calling SNPs on %s using varscan" % s
				call("samtools mpileup -d1000000 {s}/{sam_file}.sorted.bam > {s}/{sam_file}.pileup -f {reference_sequence}".format(s=s, sam_file=sam_file, reference_sequence=cfg.reference_sequence), shell=True)
				call("java -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp {s}/{sam_file}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{vcf_name}".format(s=s, vcf_name=vcf_name, sam_file=sam_file, snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True)

			if cfg.annotate_aa_changes == True:
				call("java -jar /usr/local/bin/snpEff_latest_core/snpEff/snpEff.jar {snpEff_ref_name} {s}/{vcf_name} > {s}/{annotated_name}".format(snpEff_ref_name=snpEff_ref_name,s=s,vcf_name=vcf_name, annotated_name=annotated_name, snp_frequency=cfg.snp_frequency), shell=True)


		#call("cd {s}; rm -rf snp_calls; mkdir snp_calls; for f in *.vcf; do mv $f snp_calls/$f; done".format(s=s), shell=True)

	#if cfg.use_lofreq == True and cfg.annotate_aa_changes == True:
		#call("rm combined.{snp_frequency}.snps.txt; grep -r --include='*.lofreq.annotated.{snp_frequency}.vcf' . * >> combined.lofreq.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'/#/d' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\;/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\:/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\t*\=/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\.fastq.*\.vcf//g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)


	#if cfg.use_varscan == True and cfg.annotate_aa_changes == True:
		#call("rm combined.{snp_frequency}.snps.txt; grep -r --include='*.varscan.snps.annotated.{snp_frequency}.vcf' . * >> varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'/Position/d' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\:/\t/g' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		#call("sed -i '' $'s/\.fastq.*\.txt//g' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)



# perform de novo assembly with Trinity
def de_novo_assembly():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		
		# set the trimmed fastq files to use for de novo assembly
		trimmed_reads = []
			
		for f in os.listdir(s):
			if f.endswith(".trimmed.fastq") == True:
				trimmed = s + "/" + f
				trimmed_reads.append(trimmed)
		fastqs_to_assemble = ",".join(trimmed_reads)

		call("/usr/local/bin/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --single {fastqs_to_assemble} --max_memory 3G --output {s}/trinity_output".format(s=s,fastqs_to_assemble=fastqs_to_assemble), shell=True)

		# pipe the output to blastn to return the identity of the contigs
		call("/usr/local/bin/ncbi-blast-2.6.0+/bin/blastn -db nt -query {s}/trinity_output/Trinity.fasta -out {s}/Trinity_BLAST_result.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length evalue stitle qcovhsp qstart qend sstart send' -remote".format(s=s), shell=True)


def de_novo_assemble_mapped_reads():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		
		# set the trimmed fastq files to use for de novo assembly
		trimmed_reads = []
			
		for f in os.listdir(s):
			if f.endswith(".trimmed.fastq") == True:
				trimmed = f
				trimmed_reads.append(trimmed)
		fastqs_to_assemble = ",".join(trimmed_reads)

		# cp the original and trimmed fastq files into the new directory, remove the white spaces and replace with _s
		call("mkdir {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True)
		
		for i in range(0, len(trimmed_reads)):
			fastq = trimmed_reads[i]
			call("cp {s}/{fastq} {s}/trinity_de_novo_assembly_mapped_reads_only/{fastq}.nospaces.fastq".format(s=s, fastq=fastq), shell=True)
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.fastq; do sed -i '' $'s/\ /\_/g' $f; done".format(s=s), shell=True)

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
			no_spaces_read = s + "/trinity_de_novo_assembly_mapped_reads_only/" + i + ".nospaces.fastq"
			trimmed_reads_no_spaces.append(no_spaces_read)
		fastqs_no_spaces = ",".join(trimmed_reads_no_spaces)
		
		call("rm {s}/bowtie_reference_files; rm {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True)
		call("mkdir {s}/bowtie_reference_files; cp {reference_sequence} {s}/bowtie_reference_files/{reference_name}".format(s=s,reference_sequence=reference_sequence, reference_name=reference_name), shell=True)
		call("bowtie2-build {s}/bowtie_reference_files/{reference_name} {s}/bowtie_reference_files/{reference_name}".format(s=s, reference_name=reference_name), shell=True)
		#call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
		call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {fastqs_no_spaces} -S {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam --local".format(s=s, fastqs_no_spaces=fastqs_no_spaces,reference_name=reference_name), shell=True)

		# put back the spaces so that trinity can use it
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_1\:N/\ 1\:N/g' $f; done".format(s=s), shell=True)
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_2\:N/\ 2\:N/g' $f; done".format(s=s), shell=True)

		# split sam file into forward and reverse and then extract fastqs
		call("splitsam.sh {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.unmapped.sam".format(s=s), shell=True)
		call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.fastq ow=t".format(s=s), shell=True)
		call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.fastq ow=t".format(s=s), shell=True)

		# perform de novo assembly
		call("/usr/local/bin/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --single {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.fastq,{s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.fastq --max_memory 3G --output {s}/trinity_de_novo_assembly_mapped_reads_only".format(s=s), shell=True)

		# pipe the output to blastn to return the identity of the contigs
		print "now BLAST searching contigs from Trinity assembly of mapped reads from %s" % s
		call("/usr/local/bin/ncbi-blast-2.6.0+/bin/blastn -db nt -query {s}/trinity_de_novo_assembly_mapped_reads_only/Trinity.fasta -out {s}/trinity_de_novo_assembly_mapped_reads_only/Trinity_BLAST_result.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length evalue stitle qcovhsp qstart qend sstart send' -remote".format(s=s), shell=True)



# calculate pi with popoolation
def pi_analyses():
	for s in sample_dict:

		if cfg.use_different_reference_for_each_sample == True:
			for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					gtf_name = file
					reference_name = file

			if ".fasta" in gtf_name:
				gtf_name = gtf_name.replace(".fasta", "")
			if ".fa" in gtf_name:
				gtf_name = gtf_name.replace(".fa", "")


			print gtf_name
			if "_H5_partial" in gtf_name:
				gtf_name = gtf_name.replace("_H5_partial","")
			if "full_genome" in gtf_name:
				gtf_name = gtf_name.replace("full_genome","")
			if "_mixed" in gtf_name:
				gtf_name = gtf_name.replace("_mixed","")
			print gtf_name


		elif cfg.use_different_reference_for_each_sample == False:
			reference_name = cfg.reference_sequence
			if ".fasta" in reference_name:
				gtf_name = reference_name.replace(".fasta", "")
			elif ".fa" in reference_name:
				gtf_name = reference_name.replace(".fa", "")

		# define sam file based on whether duplicate read removal and coverage normalization are turned on or off
		if cfg.remove_duplicate_reads == False and cfg.normalize_coverage == False:
			sam_file = s + ".sam"
		elif cfg.remove_duplicate_reads == True and cfg.normalize_coverage == False:
			sam_file = s + ".nodups.sam"
		elif cfg.remove_duplicate_reads == True and cfg.normalize_coverage == True:
			sam_file = s + ".nodups.normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
		elif cfg.remove_duplicate_reads == False and cfg.normalize_coverage == True:
			s + ".normalized." + str(cfg.coverage_normalization_depth) + "x.sam"
		base_name = sam_file.replace(".sam", "")
				
		call("mkdir {s}/popoolation_analyses; cp {s}/{sam_file} {s}/popoolation_analyses/{sam_file}".format(s=s, sam_file=sam_file, base_name=base_name), shell=True)
		call("samtools view -bS {s}/popoolation_analyses/{sam_file} > {s}/popoolation_analyses/{base_name}.bam; samtools sort {s}/popoolation_analyses/{base_name}.bam > {s}/popoolation_analyses/{base_name}.sorted.bam".format(s=s, sam_file=sam_file, base_name=base_name), shell=True)
		call("samtools mpileup -d 1000000 {s}/popoolation_analyses/{base_name}.sorted.bam > {s}/popoolation_analyses/{base_name}.pileup".format(s=s, base_name=base_name), shell=True)


		if cfg.calculate_genewise_pi == True:
			call("perl /usr/local/bin/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 500 --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --dissable-corrections --gtf /usr/local/bin/snpEff_latest_core/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{base_name}.pileup --output {s}/popoolation_analyses/{base_name}.pi.txt".format(s=s, base_name=base_name, reference_name=reference_name,gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count,min_quality=cfg.min_quality), shell=True)

		if cfg.calculate_genewise_piNpiS == True:
			call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-at-position.pl --measure pi --pool-size 500 --codon-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/nsl_p1.txt --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --dissable-corrections --gtf /usr/local/bin/snpEff_latest_core/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{base_name}.pileup --output {s}/popoolation_analyses/{base_name}.syn-nonsyn.txt".format(s=s, base_name=base_name, reference_name=reference_name, gtf_name=gtf_name, min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality), shell=True)

		if cfg.calculate_sliding_window_piNpiS == True:
			call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 500 --codon-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/nsl_p1.txt --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --window-size {pi_window_size} --step-size {pi_step_size} --dissable-corrections --gtf /usr/local/bin/snpEff_latest_core/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{base_name}.pileup --output {s}/popoolation_analyses/{base_name}.sliding.txt".format(s=s, base_name=base_name, min_coverage=cfg.min_coverage,  gtf_name=gtf_name, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality,pi_window_size=cfg.pi_window_size,pi_step_size=cfg.pi_step_size), shell=True)

		if cfg.perform_subsampling == True:
			call("mkdir {s}/popoolation_analyses/subsampled".format(s=s), shell=True)
			call("perl /usr/local/bin/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input {s}/popoolation_analyses/{base_name}.pileup --output {s}/popoolation_analyses/subsampled/{base_name}.{subsample_level}x.pileup --target-coverage {subsample_level} --max-coverage {max_coverage} --min-qual {min_quality} --method withoutreplace --fastq-type illumina".format(s=s,base_name=base_name, subsample_level=cfg.subsample_level,min_coverage=cfg.min_coverage,max_coverage=cfg.max_coverage, min_quality=cfg.min_quality), shell=True)

			if cfg.calculate_subsampled_pi == True:
				call("perl /usr/local/bin/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 500 --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --dissable-corrections --gtf /usr/local/bin/snpEff_latest_core/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/subsampled/{base_name}.{subsample_level}x.pileup --output {s}/popoolation_analyses/subsampled/{base_name}.{subsample_level}x.pi.txt".format(s=s, base_name=base_name, subsample_level=cfg.subsample_level,reference_name=reference_name,gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality), shell=True)

			if cfg.calculate_subsampled_piNpiS == True:
				call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-at-position.pl --measure pi --pool-size 500 --codon-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/nsl_p1.txt --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --dissable-corrections --gtf /usr/local/bin/snpEff_latest_core/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/subsampled/{base_name}.{subsample_level}x.pileup --output {s}/popoolation_analyses/subsampled/{base_name}.{subsample_level}x.syn-nonsyn.txt".format(s=s, base_name=base_name, subsample_level=cfg.subsample_level,reference_name=reference_name, gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality), shell=True)


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

combine_fastqfiles()

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

if cfg.run_popoolation == True:
	pi_analyses()
