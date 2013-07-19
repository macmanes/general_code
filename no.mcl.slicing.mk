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

ARRAY = 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 Pool12

##### No Editing should be necessary below this line  #####
CPU=2
RUN=run
PWD = `pwd`
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')
#SHELL=/bin/bash -o pipefail

all: check b2 trim correct merge assemble rsem


check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatic is Installed
	command -v $(TRINITY/Trinity.pl) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	command -v fastq-converter-v2.0.pl >/dev/null 2>&1 || { echo >&2 "I require fastq-converter-v2.0.pl (Reptile package) but it's not installed.  Aborting."; exit 1; }
	command -v reptile-omp >/dev/null 2>&1 || { echo >&2 "I require reptile-omp but it's not installed.  Aborting."; exit 1; }
	command -v reptile_merger >/dev/null 2>&1 || { echo >&2 "I require reptile_merger but it's not installed.  Aborting."; exit 1; }
	@echo Reptile is installed
	command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "I require Bowtie2 but it's not installed.  Aborting."; exit 1; }
	@echo Bowtie2 is installed"\n"

b2: 
	for SLICE in $(ARRAY); do \
		bowtie2 -t -p $(CPU) -x Dmel/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome \
		-U bin.$$SLICE.R1.fastq.gz \
		--al-gz pp.$$SLICE.1.fastq.gz \
		> /dev/null ; \
		bowtie2 -t -p $(CPU) -x Dmel/Ensembl/BDGP5.25/Sequence/Bowtie2Index/genome \
		-U bin.$$SLICE.R2.fastq.gz \
		--al-gz pp.$$SLICE.2.fastq.gz \
		> /dev/null ; \
	done

trim:
	@echo About to start trimming
	for SLICE in $(ARRAY) ; do \
		java -Xmx20g -jar $(TRIMMOMATIC) PE -phred33 -threads $(CPU) \
		pp.$$SLICE.1.fastq.gz \
		pp.$$SLICE.2.fastq.gz \
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

correct:
	@echo $$SHELL
	for SLICE in $(ARRAY); do \
		fastq-converter-v2.0.pl ./ ./ 1 ;  \
		sed -i 's\[0-9]$$\&/1\' left.fa ; \
		sed -i 's\[0-9]$$\&/2\g' right.fa ;\
		sed -i 's\^>.*[0-9]$$\&/1\g' left.q ;\
		sed -i 's\^>.*[0-9]$$\&/2\g' right.q ;\
		cat left.fa right.fa > both.fa ;\
		cat left.q right.q > both.q ;\
		reptile-omp $(CONFIG) ;\
		rm left.fa right.fa left.q right.q ;\
	done

merge: both.reptile.err
	for SLICE in $(ARRAY); do \
		reptile_merger both.fa $< both.reptile.corr.fa ; \
		grep -aA1 '/1' both.reptile.corr.fa > left.$$SLICE.rept.corr.fa ; \
		grep -aA1 '/2' both.reptile.corr.fa > right.$$SLICE.rept.corr.fa ; \
		rm both.reptile.corr.fa ;\
		gzip left.$$SLICE.rept.corr.fa ;\
		gzip right.$$SLICE.rept.corr.fa ;\
	done	

assemble: 
	for SLICE in $(ARRAY); do \
			$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 \
			--seqType fq --JM 20G \
			--left left.$$SLICE.fastq \
			--right right.$$SLICE.fastq --CPU $(CPU) --output $$SLICE; done

rsem: $(SLICE).Trinity.fasta
	$(TRINITY)/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $< --seqType fq --left $(READ1) \
	--right $(READ2) --thread_count $(CPU) --SS_lib_type RF -- --bowtie-chunkmbs 512


