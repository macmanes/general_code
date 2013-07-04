#! /bin/bash

ARRAY=5

##### No Editing should be necessary below this line  #####
CPU=8
RUN=run


for SLICE in $ARRAY; do
	if [ $SLICE -eq 5 ];
	then
		MIN=300
		/home/macmanes/trinityrnaseq_r2013-02-25/Trinity.pl --full_cleanup --min_kmer_cov 2 \
		--seqType fq --JM 10G --min_contig_length $MIN \
		--left ~/trinityrnaseq-code/trunk/sample_data/test_Trinity_Assembly/reads.left.fq \
		--right ~/trinityrnaseq-code/trunk/sample_data/test_Trinity_Assembly/reads.right.fq --CPU $CPU --output $SLICE
	else
		echo ""
	fi \
done

