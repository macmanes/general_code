#!/usr/bin/make -rRsf

###########################################
###        -usage 'assembly.mk RUN=run CPU=8 READ1=/location/of/read1.fastq READ2=/location/of/read2.fastq'
###         -RUN= name of run
###
###        -limitations=  must use PE files.. no support for SE...
###
###         -Make sure your Trinity base directory 
###         	is set properly
###         -Make sure barcode file is located with
###           BCODES= tag
###          -Make sure you pull the config.analy file from GIT, or make you own.
############################################
TRINITY := /home/macmanes/trinityrnaseq-code/trunk
BCODES := /home/macmanes/Dropbox/barcodes.fa
CONFIG:= /home/macmanes/Dropbox/config.analy



##### No Editing should be necessary below this line  #####




TRIM=5
CPU=2
RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

.PHONY: check clean
all: check $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq both.reptile.err $(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa $(RUN).Trinity.fasta RSEM.genes.results
trim: check $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq
correct: check both.reptile.err
merge: check $(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa
assemble: check $(RUN).Trinity.fasta
rsem: check RSEM.genes.results

check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatic is Installed
	command -v $(TRINITY/Trinity.pl) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	command -v fastq-converter-v2.0.pl >/dev/null 2>&1 || { echo >&2 "I require fastq-converter-v2.0.pl (Reptile package) but it's not installed.  Aborting."; exit 1; }
	command -v reptile-omp >/dev/null 2>&1 || { echo >&2 "I require reptile-omp but it's not installed.  Aborting."; exit 1; }
	command -v reptile_merger >/dev/null 2>&1 || { echo >&2 "I require reptile_merger but it's not installed.  Aborting."; exit 1; }
	@echo Reptile is installed"\n"
	if [ -f $(READ1) ]; then echo 'left fastQ exists'; else echo 'Im having trouble finding your left fastQ file, check PATH \n'; exit 1; fi;
	if [ -f $(READ2) ]; then echo 'right fastQ exists \n'; else echo 'Im having trouble finding your right fastQ file, check PATH \n'; exit 1; fi;
$(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq: $(READ1) $(READ2)
	@echo About to start trimming
		java -Xmx30g -jar $(TRIMMOMATIC) PE -phred33 -threads $(CPU) \
		$(READ1) \
		$(READ2) \
		$(RUN).pp.1.fq \
		$(RUN).up.1.fq \
		$(RUN).pp.2.fq \
		$(RUN).up.2.fq \
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$(TRIM) TRAILING:$(TRIM) SLIDINGWINDOW:4:$(TRIM) MINLEN:25 ; 
		cat $(RUN).pp.1.fq $(RUN).up.1.fq > $(RUN)_left.$(TRIM).fastq ; 
		cat $(RUN).pp.2.fq $(RUN).up.2.fq > $(RUN)_right.$(TRIM).fastq ; 


both.reptile.err: $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq
	@echo About to start error correction
	fastq-converter-v2.0.pl ./ ./ 1 #files MUST have fastq extension
	sed -i 's_[0-9]$$_&/1_' $(RUN)_left.fa #add /1 to ID reads as left
	sed -i 's_[0-9]$$_&/2_' $(RUN)_right.fa #add /2 to ID reads as right
	sed -i 's_^>.*[0-9]$$_&/1_' $(RUN)_left.q
	sed -i 's_^>.*[0-9]$$_&/2_' $(RUN)_right.q
	cat $(RUN)_left.fa $(RUN)_right.fa > both.fa
	cat $(RUN)_left.q $(RUN)_right.q > both.q
	reptile-omp $(CONFIG) #Do error corection


$(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa: both.reptile.err
	reptile_merger both.fa $< both.reptile.corr.fa #make error corrected fasta file
	grep -aA1 '/1' both.reptile.corr.fa > $(RUN).left.rept.corr.fa
	grep -aA1 '/2' both.reptile.corr.fa > $(RUN).right.rept.corr.fa
	

$(RUN).Trinity.fasta:  $(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fa --JM 10G \
	--left $(RUN).left.rept.corr.fa --right $(RUN).right.rept.corr.fa --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)
	
RSEM.genes.results: $(RUN).Trinity.fasta
	$(TRINITY)/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $< --seqType fq --left $(READ1) \
	--right $(READ2) --thread_count $(CPU) --SS_lib_type RF -- --bowtie-chunkmbs 512

clean: 
	rm TRANS*
	rm *bam
	rm both.reptile.corr.fa

