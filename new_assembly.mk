#!/usr/bin/make -rRsf

###########################################
###        -usage 'assembly.mk RUN=run CPU=8 MEM=15 READ1=/location/of/read1.fastq READ2=/location/of/read2.fastq'
###         -RUN= name of run
###
###        -limitations=  must use PE files.. no support for SE...
###
###         -Make sure your Trinity base directory 
###         	is set properly
###         -Make sure barcode file is located with
###           BCODES= tag
###          
############################################
TRINITY := /home/macmanes/trinityrnaseq-code/trunk
BCODES := /home/macmanes/Dropbox/barcodes.fa



##### No Editing should be necessary below this line  #####

MINLEN=25
PHRED=33
SEQ=fq
MINK=1
MEM=2
TRIM=5
CPU=2
BCPU=$(CPU)
RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

.PHONY: check clean
all: check $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq $(RUN).Trinity.fasta $(RUN).xprs
trim: check $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq
assemble: check $(RUN).Trinity.fasta
express: check $(RUN).xprs

check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatic is Installed
	command -v $(TRINITY/Trinity.pl) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	if [ -f $(READ1) ]; then echo 'left fastQ exists'; else echo 'Im having trouble finding your left fastQ file, check PATH \n'; exit 1; fi;
	if [ -f $(READ2) ]; then echo 'right fastQ exists \n'; else echo 'Im having trouble finding your right fastQ file, check PATH \n'; exit 1; fi;

$(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq: $(READ1) $(READ2)
	@echo About to start trimming
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE -phred$(PHRED) -threads $(CPU) \
		$(READ1) \
		$(READ2) \
		$(RUN).pp.1.fq \
		$(RUN).up.1.fq \
		$(RUN).pp.2.fq \
		$(RUN).up.2.fq \
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$(TRIM) TRAILING:$(TRIM) SLIDINGWINDOW:4:$(TRIM) MINLEN:$(MINLEN) ; 
		cat $(RUN).pp.1.fq $(RUN).up.1.fq > $(RUN)_left.$(TRIM).fastq ; 
		cat $(RUN).pp.2.fq $(RUN).up.2.fq > $(RUN)_right.$(TRIM).fastq ; 
	
$(RUN).Trinity.fasta: $(RUN)_left.$(TRIM).fastq $(RUN)_right.$(TRIM).fastq
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov $(MINK) --seqType $(SEQ) --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G --bflyCPU $(BCPU) \
	--left $(RUN)_left.$(TRIM).fastq --right $(RUN)_right.$(TRIM).fastq --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)
	
$(RUN).xprs: $(RUN).Trinity.fasta
		@echo Quantitiating Transcripts
		bwa index -p index $(RUN).Trinity.fasta
		bwa mem -t$(CPU) index $(READ1) $(READ2) \
		2>>first.log | samtools view -@6 -Sub - | tee >(samtools flagstat - > $(RUN).mapping.stats) |  express -o $(RUN).xprs -p $(CPU) $(RUN).Trinity.fasta 2>>first.log

clean: 
	rm *TRANS*
	rm *bam *.q *err
