#!/usr/bin/make -rRsf

###########################################
###        -usage 'slicing.mk RUN=run CPU=8 READ1=/location/of/read1.fastq READ2=/location/of/read2.fastq'
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
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
BCODES := /home/macmanes/Dropbox/barcodes.fa
CONFIG:= /home/macmanes/Dropbox/config.analy



##### No Editing should be necessary below this line  #####
CPU=2
RUN=run
PWD = `pwd`
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

all: check trim correct merge assemble rsem

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

trim: $(READ1) $(READ2)
	
	@echo About to start trimming
	for SLICE in 3 4 5 6 7 8 9 10; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE -phred33 -threads $(CPU) \
		bin_$(SLICE)_R1.fastq.gz \
		bin_$(SLICE)_R2.fastq.gz \
		T.$$SLICE.pp.1.fq \
		T.$$SLICE.up.1.fq \
		T.$$SLICE.pp.2.fq \
		T.$$SLICE.up.2.fq \
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:25 ; \
		cat T.$$SLICE.pp.1.fq T.$$SLICE.up.1.fq > left.$$SLICE.fastq ; \
		cat T.$$SLICE.pp.2.fq T.$$SLICE.up.2.fq > right.$$SLICE.fastq ; \
		rm T.$$SLICE.pp.2.fq T.$$SLICE.up.2.fq T.$$SLICE.pp.1.fq T.$$SLICE.up.1.fq ; \
	done

correct: left.$$SLICE.fastq right.$$SLICE.fastq
	@echo About to start error correction
	perl fastq-converter-v2.0.pl ./ ./ 1 #files MUST have fastq extension
	sed -i 's_[0-9]$_&/1_' $(SLICE)_left.fa #add /1 to ID reads as left
	sed -i 's_[0-9]$_&/2_' $(SLICE)_right.fa #add /2 to ID reads as right
	sed -i 's_^>.*[0-9]$_&/1_' $(SLICE)_left.q
	sed -i 's_^>.*[0-9]$_&/2_' $(SLICE)_right.q
	cat $(SLICE)_left.fa $(SLICE)_right.fa > both.fa
	cat $(SLICE)_left.q $(SLICE)_right.q > both.q
	reptile-omp $(CONFIG) #Do error corection

merge: both.reptile.err
	reptile_merger both.fa $< both.reptile.corr.fa #make error corrected fasta file
	grep -aA1 '/1' both.reptile.corr.fa > $(SLICE).left.rept.corr.fa
	grep -aA1 '/2' both.reptile.corr.fa > $(SLICE).right.rept.corr.fa
	

assemble:  $(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa
	$(TRINITY)/Trinity.pl --full_cleanup --SS_lib_type RF --min_kmer_cov 2 --seqType fa --JM 30G \
	--left $(RUN).left.rept.corr.fa --right $(RUN).right.rept.corr.fa --CPU $(CPU) --output $(RUN)
	
rsem: $(RUN).Trinity.fasta
	$(TRINITY)/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $< --seqType fq --left $(READ1) \
	--right $(READ2) --thread_count $(CPU) --SS_lib_type RF -- --bowtie-chunkmbs 512


