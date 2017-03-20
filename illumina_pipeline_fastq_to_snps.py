#!/usr/bin/env python

import sys, subprocess, glob, os, shutil, re, importlib
from subprocess import call

config_filename = "config"   #sys.argv[1]
config = __import__(config_filename)						# imports config_filename as a module
cfg = config.configuration()								# make cfg an instance of class configuration accessed from the config module

file_list = []
sample_list = []
sample_dict = {}
reference_list = []


# make a list of all the fastq files in the directory and store their sample name as the fastq file name without the R1_001 or R2_001. 
for file in glob.glob("*.fastq"):						# glob will find instances of all files ending in .fastq as shell would
	file_list.append(file)
	
for file in file_list:									# find all reads that are pairs
	if "R1" in file:
		samplename = file.replace("_R1_001", "")		# rename R1 reads to take out the R1_001.fastq
		samplename = samplename.replace(".fastq", "")
	elif "R2" in file:
		samplename = file.replace("_R2_001", "")		# rename R2 reads to take out the R2_001.fastq
		samplename = samplename.replace(".fastq", "")	
	else:
		samplename=file									# if there are only nonpaired reads, just keep the same name
			
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
			if cfg.paired_trim == False:
				call("java -jar /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE {s}/{value} {s}/{value}.trimmed.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True)
			
			elif cfg.paired_trim == True: 
				call("java -jar /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE {s}/{value} {s}/{value} -baseout {s}/{s}.fastq SLIDINGWINDOW:{window}:{qscore} MINLEN:{MINLEN}".format(s=s, value=value, MINLEN=cfg.minlength, window=cfg.window_size,qscore=cfg.trim_qscore), shell=True)

	

# perform mapping
def map(sample_list):
	
	if cfg.use_different_reference_for_each_sample == False:
		call("bowtie2-build {reference_sequence} {reference_sequence_name}".format(reference_sequence=cfg.reference_sequence, reference_sequence_name=cfg.reference_sequence_name), shell=True)
		for s in sample_dict:
			call("bowtie2 -x {reference_sequence} -U {s}/{f1}.trimmed.fastq,{s}/{f2}.trimmed.fastq -S {s}/{s}.sam --local".format(s=s, reference_sequence=cfg.reference_sequence, f1=sample_dict[s][0], f2=sample_dict[s][1]), shell=True)
	
	
	elif cfg.use_different_reference_for_each_sample == True:		
		for s in sample_dict:
			for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					reference_name = file
			call("rm {s}/bowtie_reference_files".format(s=s), shell=True)
			call("bowtie2-build {s}/{reference_name} {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True)
			call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
			call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {s}/{f1}.trimmed.fastq,{s}/{f2}.trimmed.fastq -S {s}/{s}.sam --local".format(s=s, reference_name=reference_name, f1=sample_dict[s][0], f2=sample_dict[s][1]), shell=True)
		

def call_SNPs():
	for s in sample_dict:
		#convert sam to bam, then sort bam file
		call("samtools view -bS {s}/{s}.sam > {s}/{s}.bam".format(s=s), shell=True)
		call("samtools sort {s}/{s}.bam > {s}/{s}.sorted.bam".format(s=s), shell=True)
		call("rm {s}/{s}.lofreq.{freq}.vcf; rm {s}/{s}.filtered.{freq}.vcf".format(s=s,freq=cfg.snp_frequency), shell=True)
		
		# assign names to snpEff_ref_name variable, which will be used for amino acid annotation with snpEff
		if ".fasta" in cfg.reference_sequence_name: 
			snpEff_ref_name = cfg.reference_sequence_name.replace(".fasta", "")
		if ".fa" in cfg.reference_sequence_name: 
			snpEff_ref_name = cfg.reference_sequence_name.replace(".fa", "")

	
		# run SNP calling and amino acid change annotations if samples are all mapped to the same reference
		if cfg.use_different_reference_for_each_sample == False:
			# call variants with lofreq and filter them
			if cfg.use_lofreq == True:
				print "now calling SNPs on %s using lofreq" % s
				call("lofreq call -f {reference_sequence} -o {s}/{s}.lofreq.{snp_frequency}.vcf {s}/{s}.sorted.bam".format(s=s, snp_frequency=cfg.snp_frequency, reference_sequence=cfg.reference_sequence), shell=True)	
				call("lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {s}/{s}.lofreq.{snp_frequency}.vcf -o {s}/{s}.filtered.{snp_frequency}.vcf".format(s=s, min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold, snp_frequency=cfg.snp_frequency), shell=True)
				
			if cfg.annotate_aa_changes == True:
				call("java -jar /usr/local/bin/snpEff/snpEff.jar {snpEff_ref_name} {s}/{s}.lofreq.filtered.{snp_frequency}.vcf > {s}/{s}.lofreq.annotated.{snp_frequency}.vcf".format(snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True)

			if cfg.use_varscan == True:
				print "now calling SNPs on %s using varscan" % s
				call("samtools mpileup -d1000000 {s}/{s}.sorted.bam > {s}/{s}.pileup -f {reference_sequence}".format(s=s, reference_sequence=cfg.reference_sequence), shell=True)
				call("java -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp {s}/{s}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{s}.varscan.snps.{snp_frequency}.vcf".format(s=s,snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True)		
				
			if cfg.annotate_aa_changes == True:
				call("java -jar /usr/local/bin/snpEff/snpEff.jar {snpEff_ref_name} {s}/{s}.varscan.snps.{snp_frequency}.vcf > {s}/{s}.varscan.snps.annotated.{snp_frequency}.vcf".format(snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True)
	
		
		
		# run SNP calling and amino acid change annotations if samples are mapped to their own references
		elif cfg.use_different_reference_for_each_sample == True:		
			for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					reference_sequence = file
			
			# assign snpEff reference names
			if ".fasta" in reference_sequence: 
				snpEff_ref_name = reference_sequence.replace(".fasta", "")
			elif ".fa" in reference_sequence: 
				snpEff_ref_name = reference_sequence.replace(".fa", "")				
			
			if cfg.use_lofreq == True: 
				print "now calling SNPs on %s using lofreq" % s
				call("lofreq call -f {s}/{reference_sequence} -o {s}/{s}.lofreq.{snp_frequency}.vcf {s}/{s}.sorted.bam".format(s=s, snp_frequency=cfg.snp_frequency, reference_sequence=reference_sequence), shell=True)	
				call("lofreq filter --cov-min {min_cov} --snvqual-thresh {snp_qual_threshold} --af-min {snp_frequency} -i {s}/{s}.lofreq.{snp_frequency}.vcf -o {s}/{s}.lofreq.filtered.{snp_frequency}.vcf".format(s=s, min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold, snp_frequency=cfg.snp_frequency), shell=True)
				
			if cfg.annotate_aa_changes == True:
				call("java -jar /usr/local/bin/snpEff/snpEff.jar {snpEff_ref_name} {s}/{s}.lofreq.filtered.{snp_frequency}.vcf > {s}/{s}.lofreq.annotated.{snp_frequency}.vcf".format(snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True)
	
			if cfg.use_varscan == True: 
				print "now calling SNPs on %s using varscan" % s
				call("samtools mpileup -d 1000000 {s}/{s}.sorted.bam > {s}/{s}.pileup -f {s}/{reference_sequence}".format(s=s, reference_sequence=reference_sequence), shell=True)
				call("java -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp {s}/{s}.pileup --min-coverage {min_cov} --min-avg-qual {snp_qual_threshold} --min-var-freq {snp_frequency} --strand-filter 1 --output-vcf 1 > {s}/{s}.varscan.snps.{snp_frequency}.vcf".format(s=s,snp_frequency=cfg.snp_frequency,min_cov=cfg.min_cov, snp_qual_threshold=cfg.snp_qual_threshold), shell=True)
				
			if cfg.annotate_aa_changes == True:
				call("java -jar /usr/local/bin/snpEff/snpEff.jar {snpEff_ref_name} {s}/{s}.varscan.snps.{snp_frequency}.vcf > {s}/{s}.varscan.snps.annotated.{snp_frequency}.vcf".format(snpEff_ref_name=snpEff_ref_name,s=s,snp_frequency=cfg.snp_frequency), shell=True)

				
		#call("cd {s}; rm -rf snp_calls; mkdir snp_calls; for f in *.vcf; do mv $f snp_calls/$f; done".format(s=s), shell=True)
		
	if cfg.use_lofreq == True and cfg.annotate_aa_changes == True: 			
		call("rm combined.{snp_frequency}.snps.txt; grep -r --include='*.lofreq.annotated.{snp_frequency}.vcf' . * >> combined.lofreq.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'/#/d' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\;/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\:/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\t*\=/\t/g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\.fastq.*\.vcf//g' combined.{snp_frequency}.snps.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		
	
	if cfg.use_varscan == True and cfg.annotate_aa_changes == True: 
		call("rm combined.{snp_frequency}.snps.txt; grep -r --include='*.varscan.snps.annotated.{snp_frequency}.vcf' . * >> varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'/Position/d' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\:/\t/g' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)
		call("sed -i '' $'s/\.fastq.*\.txt//g' varscan.snps.{snp_frequency}.txt".format(snp_frequency=cfg.snp_frequency), shell=True)



# perform de novo assembly with Trinity
def de_novo_assembly():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		fastq1 = sample_dict[s][0]+".trimmed.fastq"
		fastq2 = sample_dict[s][1]+".trimmed.fastq"
		call("/usr/local/bin/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --single {s}/{fastq1},{s}/{fastq2} --max_memory 3G --output {s}/trinity_output".format(s=s,fastq1=fastq1, fastq2=fastq2), shell=True)

		# pipe the output to blastn to return the identity of the contigs
		call("/usr/local/bin/ncbi-blast-2.6.0+/bin/blastn -db nt -query {s}/trinity_output/Trinity.fasta -out {s}/Trinity_BLAST_result.txt -max_target_seqs 10 -outfmt '7 qseqid sseqid pident length evalue stitle qcovhsp qstart qend sstart send' -remote".format(s=s), shell=True)
	

def de_novo_assemble_mapped_reads():
	# perform de novo assembly using trinity on the trimmed fastq files
	for s in sample_dict:
		fastq1 = sample_dict[s][0]+".trimmed.fastq"
		fastq2 = sample_dict[s][1]+".trimmed.fastq"	
		
		# cp the original and trimmed fastq files into the new directory, remove the white spaces and replace with _s
		call("mkdir {s}/trinity_de_novo_assembly_mapped_reads_only; cp {s}/{fastq1} {s}/trinity_de_novo_assembly_mapped_reads_only/{fastq1}.nospaces.fastq; cp {s}/{fastq2} {s}/trinity_de_novo_assembly_mapped_reads_only/{fastq2}.nospaces.fastq".format(s=s, fastq1=fastq1, fastq2=fastq2), shell=True)
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.fastq; do sed -i '' $'s/\ /\_/g' $f; done".format(s=s), shell=True)
		
		# redo bowtie mapping and put the white spaces back into the sam file
		for file in os.listdir(s):
				if file.endswith(".fasta") or file.endswith(".fa"):
					reference_name = file
		call("rm {s}/bowtie_reference_files".format(s=s), shell=True)
		call("bowtie2-build {s}/{reference_name} {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True)
		call("mkdir {s}/bowtie_reference_files; cd {s}/; for f in *.bt2; do mv $f bowtie_reference_files/$f; done".format(s=s), shell=True)
		call("bowtie2 -x {s}/bowtie_reference_files/{reference_name} -U {s}/trinity_de_novo_assembly_mapped_reads_only/{fastq1}.nospaces.fastq,{s}/trinity_de_novo_assembly_mapped_reads_only/{fastq2}.nospaces.fastq -S {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam --local".format(s=s, fastq1=fastq1, fastq2=fastq2,reference_name=reference_name), shell=True)
		
		# put back the spaces so that trinity can use it
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_1\:N/\ 1\:N/g' $f; done".format(s=s), shell=True)
		call("for f in {s}/trinity_de_novo_assembly_mapped_reads_only/*.sam; do sed -i '' $'s/\_2\:N/\ 2\:N/g' $f; done".format(s=s), shell=True)
		
		# split sam file into forward and reverse and then extract fastqs
		call("splitsam.sh {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam {s}/trinity_de_novo_assembly_mapped_reads_only/{s}.unmapped.sam".format(s=s), shell=True)
		call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.forward.fastq".format(s=s), shell=True)		
		call("reformat.sh in={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.sam out={s}/trinity_de_novo_assembly_mapped_reads_only/{s}.reverse.fastq".format(s=s), shell=True)		

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
				if file.endswith(".fasta"):
					gtf_name = file.replace(".fasta", "")
					reference_name = file
				elif file.endswith(".fa"):
					gtf_name = file.replace(".fa", "")
					reference_name = file
					
		elif cfg.use_different_reference_for_each_sample == False: 
			reference_name = cfg.reference_sequence
			if ".fasta" in reference_name:
				gtf_name = reference_name.replace(".fasta", "")
			elif ".fa" in reference_name: 
				gtf_name = reference_name.replace(".fa", "")
				
		call("mkdir {s}/popoolation_analyses; cp {s}/{s}.sam {s}/popoolation_analyses/{s}.sam".format(s=s), shell=True)
		call("samtools view -bS {s}/popoolation_analyses/{s}.sam > {s}/popoolation_analyses/{s}.bam; samtools sort {s}/popoolation_analyses/{s}.bam > {s}/popoolation_analyses/{s}.sorted.bam".format(s=s), shell=True)
		call("samtools mpileup -d 1000000 {s}/popoolation_analyses/{s}.sorted.bam > {s}/popoolation_analyses/{s}.pileup -f {s}/{reference_name}".format(s=s, reference_name=reference_name), shell=True)


		if cfg.subsample == False: 
			if cfg.calculate_genewise_pi == True:
				call("perl /usr/local/bin/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 500 --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --gtf /usr/local/bin/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{s}.pileup --output {s}/popoolation_analyses/{s}.pi.txt".format(s=s, reference_name=reference_name,gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count,min_quality=cfg.min_quality), shell=True)
			
			if cfg.calculate_genewise_piNpiS == True: 
				call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-at-position.pl --measure pi --pool-size 500 --codon-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/{nonsyn_length_table} --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --dissable-corrections --gtf /usr/local/bin/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{s}.pileup --output {s}/popoolation_analyses/{s}.syn-nonsyn.txt".format(s=s, reference_name=reference_name, gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality, nonsyn_length_table=cfg.nonsyn_length_table), shell=True)
		
			if cfg.calculate_sliding_window_piNpiS == False:
				call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 500 --codon-table syn-nonsyn/codon-table.txt --nonsyn-length-table syn-nonsyn/{nonsyn_length_table} --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --window-size {pi_window_size} --step-size {pi_step-size} --dissable-corrections --gtf /usr/local/bin/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/{s}.pileup --output {s}/popoolation_analyses/{s}.sliding.txt".format(s=s, min_coverage=cfg.min_coverage,  gtf_name=gtf_name, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality,nonsyn_length_table=cfg.nonsyn_length_table,pi_window_size=cfg.pi_window_size,pi_step_size=cfg.pi_step_size), shell=True)

		if cfg.subsample == True: 
			call("mkdir {s}/popoolation_analyses/subsampled".format(s=s), shell=True)
			call("perl /usr/local/bin/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input {s}/popoolation_analyses/{s}.pileup --output {s}/popoolation_analyses/subsampled/{s}.{subsample_level}x.pileup --target-coverage {subsample_level} --max-coverage {max_coverage} --min-qual {min_quality} --method withoutreplace --fastq-type illumina".format(s=s,subsample_level=cfg.subsample_level,min_coverage=cfg.min_coverage,max_coverage=cfg.max_coverage, min_quality=cfg.min_quality), shell=True)
			
			if cfg.calculate_genewise_pi == True:
				call("perl /usr/local/bin/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 500 --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual {min_quality} --gtf /usr/local/bin/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/subsampled/{s}.{subsample_level}x.pileup --output {s}/popoolation_analyses/subsampled/{s}.{subsample_level}x.pi.txt".format(s=s, subsample_level=cfg.subsample_level,reference_name=reference_name,gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality), shell=True)

			if cfg.calculate_genewise_piNpiS == True:
				call("perl /usr/local/bin/popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-at-position.pl --measure pi --pool-size 500 --codon-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /usr/local/bin/popoolation_1.2.2/syn-nonsyn/{nonsyn_length_table} --min-count {min_count} --min-coverage {min_coverage} --max-coverage {max_coverage} --min-qual 20 --dissable-corrections --gtf /usr/local/bin/snpEff/data/{gtf_name}/genes.gtf --pileup {s}/popoolation_analyses/subsampled/{s}.{subsample_level}x.pileup --output {s}/popoolation_analyses/subsampled/{s}.{subsample_level}.syn-nonsyn.txt".format(s=s, subsample_level=cfg.subsample_level,reference_name=reference_name, gtf_name=gtf_name,min_coverage=cfg.min_coverage, max_coverage=cfg.max_coverage,min_count=cfg.min_count, min_quality=cfg.min_quality, nonsyn_length_table=cfg.nonsyn_length_table), shell=True)
	

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

if cfg.call_snps == True:
	call_SNPs()
	
if cfg.de_novo_assemble_mapped_reads == True:
	de_novo_assemble_mapped_reads()

create_parameter_file()

if cfg.de_novo_assembly == True:
	de_novo_assembly()
	
	
if cfg.calculate_genewise_pi == True or cfg.calculate_genewise_piNpiS == True or cfg.calculate_sliding_window_piNpiS == True:
	pi_analyses()