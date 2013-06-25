#!/bin/bash

##Raw Reads
export FLUX_MEM="25G"
cd /media/macmanes/hd/simulated_reads
if [ ! -f right.fq ] ; 
then
	echo "Flux Sim @ `date`."
	~/flux-simulator-1.2.1/bin/flux-simulator -p test5.par -l -s -x --threads 8
	cat test5.fastq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/1$ > left.fq
	cat test5.fastq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/2$ > right.fq
else
	echo "I've already simulated the reads"
fi
chmod -w right.fq
chmod -w left.fq 
mkdir raw.reads
cd raw.reads/
if [ ! -f raw.Trinity.fasta ] ; 
then
	rm -fr raw
	echo "I'm starting the raw trin assembly"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left ../left.fq --right ../right.fq --full_cleanup --CPU 8 --output raw > raw.trin.out
else
	echo "I've already done the raw read assembly"
fi	
###Reptile
echo "##################################################"
echo "Starting Reptile @ `date`."
echo "##################################################"
mkdir ../reptile
cd ../reptile
mkdir data
cp ../right.fq data/right.fastq
cp ../left.fq data/left.fastq
if [ ! -f right.reptile.corr.fa.gz ] ; 
then
	perl ~/reptile-v1.1/reptile-v1.1/utils/fastq-converter-v2.0.pl data/ data/ 1 #Must have .fastQ extension
	~/reptile-omp/reptile-omp config.run
	~/reptile-v1.1/reptile-v1.1/utils/reptile_merger/reptile_merger data/right.fa data/right.reptile.err right.reptile.corr.fa
	~/reptile-omp/reptile-omp config.run1
	~/reptile-v1.1/reptile-v1.1/utils/reptile_merger/reptile_merger data/left.fa data/left.reptile.err left.reptile.corr.fa
else
	echo "I've already done the Rept corr"
fi
if [ ! -f reptile.Trinity.fasta ] ; 
then
	rm -fr reptile
	echo "Starting trinity"
	sed -i 's_^>[0-9].*_&/2_g' right.reptile.corr.fa
	sed -i 's_^>[0-9].*_&/1_g' left.reptile.corr.fa
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left left.reptile.corr.fa --right right.reptile.corr.fa --full_cleanup --CPU 8 --output reptile > reptile.trin.out
	gzip right.reptile.corr.fa left.reptile.corr.fa
#	rm data/*q data/*err
else
	echo "I've already done the Rept assembo"
fi
for i in `find * -type f -size -95c`; do rm $i; done

###SGA #Need to fo trinity assembly

echo "##################################################"
echo "Starting SGA @ `date`."
echo "##################################################"

mkdir ../sga
cd ../sga
if [ ! -f sga.corr.left.fq.gz ] ; 
then
	echo "I'm starting the SGA correction"	
	/home/macmanes/sga/src/SGA/sga preprocess --pe-mode 1 -m 25 -o sga.fq ../left.fq ../right.fq
	/home/macmanes/sga/src/SGA/sga index -a ropebwt -t 8 --no-reverse sga.fq
	/home/macmanes/sga/src/SGA/sga correct -k 25 --learn -t 8 --metrics=corr.metrics -o sga.corr.fq sga.fq > raw.sga.out
	cat sga.corr.fq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/1$ > sga.corr.left.fq
	cat sga.corr.fq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/2$ > sga.corr.right.fq
	rm sga.fq sga.corr.fq sga.sai sga.bwt
else
	echo "I've already done the SGA correction"
fi
gzip sga.corr.right.fq sga.corr.left.fq &
if [ ! -f sga.Trinity.fasta ] ; 
then
	rm -fr sga
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left sga.corr.left.fq.gz --right sga.corr.right.fq.gz --full_cleanup --CPU 8 --output sga > sga.trin.out
else
	echo "I've already done the SGA assemblt"	
fi

####Seecer
echo "##################################################"
echo "Starting seecer @ `date`."
echo "##################################################"
mkdir ../seecer
cd ../seecer


if [ ! -f right.fq_corrected.fa.gz ] ; 
then
	sh ~/seecer/SEECER-0.1.2/SEECER/bin/run_seecer.sh -t /media/macmanes/hd/flux/seecer/ ../left.fq ../right.fq > seecer.out 2>&1
	mv ../right.fq_corrected.fa .
	mv ../left.fq_corrected.fa .
else
	echo "I've already done the seecer correction"
fi
if [ ! -f seecer.Trinity.fasta ] ; 
then
	rm -fr seecer
	echo "I'n starting the seecer assembly"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left left.fq_corrected.fa --right right.fq_corrected.fa --full_cleanup --CPU 8 --output seecer > seecer.trin.out
else
	echo "I've already done the seecer assembly"
fi
for i in `find * -type f -size -95c`; do rm $i; done
gzip left.fq_corrected.fa right.fq_corrected.fa &
#rm counts_17_3 corrected.fasta ../*N


####AllPaths #Need to do trinity assembly
echo "##################################################"
echo "Starting AllPaths @ `date`."
echo "##################################################"
mkdir ../AllPaths
cd ../AllPaths

if [ ! -f corr.left.fasta.gz ] ; 
then
	python ~/trinityrnaseq_r2013-02-25/util/ec_norm/preproc.py -m 30 -t 8 -l ../left.fq -r ../right.fq -o corr -H True --error_corr True
	sed -i 's_read_unpaired_g' corr.unpaired.fa
	cat corr.left.fa corr.unpaired.fa > corr.left.fasta
	mv corr.right.fa corr.right.fasta
	rm corr.left.fa corr.unpaired.fa
else
	echo "I've already done the AllPaths Corr"
fi
if [ ! -f allp.Trinity.fasta ] ; 
then
	rm -fr allp
	echo "I'm starting the allp assembly"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left corr.left.fasta --right corr.right.fasta --full_cleanup --CPU 8 --output allp > allp.trin.out
else
	echo "I've already done the AllPaths assembly"
fi
for i in `find * -type f -size -95c`; do rm $i; done
gzip corr.right.fasta corr.left.fasta
####Echo
echo "##################################################"
echo "Starting Echo @ `date`."
echo "##################################################"

#cd ../echo

#if [ ! -f right.echo.fastq.gz ] ; 
#then
#	python /home/macmanes/echo_v1_12/ErrorCorrection.py -b 2000000 --nh 2048 --ncpu 8  -o right.echo.fastq right.fastq #input has to be .fastQ, not fq.. grrrrr!
#	gzip right.fastq
#else
#	echo "echo right done"
#fi
#if [ ! -f left.echo.fastq.gz ] ; 
#then
#	python /home/macmanes/echo_v1_12/ErrorCorrection.py -b 2000000 --nh 2048 --ncpu 8  -o left.echo.fastq left.fastq
#	gzip left.fastq
#else
#	echo "left done"
#fi
#if [ ! -f corr.echo.Trinity.fasta ] ; 
#then
#	rm -fr echo.corr
#	echo "Starting trinity"
#	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left left.echo.fastq --right right.echo.fastq --full_cleanup --CPU 8 --output echo > echo.trin.out
#	gzip left.echo.fastq right.echo.fastq
#else
#	echo "echo assembly done"
#fi


##Raw Reads MOUSE
echo "################################################################################################################################################################"
echo "Starting MOUSE Raw @ `date`."
echo "################################################################################################################################################################"

#cd ../raw.reads/
#if [ ! -f raw.mouse.reads.Trinity.fasta ] ; 
#then
#	rm -fr raw.mouse.reads
#	echo "Starting trinity"
#	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left /media/macmanes/hd1/mouse/raw_reads/out_1.fastq --right /media/macmanes/hd1/mouse/raw_reads/out_2.fastq --full_cleanup --CPU 8 --output raw.mouse.reads > mouse.trin.raw.out
#else
#	echo "Ive done the raw mouse read assembly"
#fi
#for i in `find * -type f -size -95c`; do rm $i; done
echo "##################################################"
echo "Starting mouse SGA @ `date`."
echo "##################################################"


cd ../sga
if [ ! -f sga.mouse.corr.right.fq.gz ] ; 
then
	/home/macmanes/sga/src/SGA/sga preprocess --pe-mode 1 -o sga.fq /media/macmanes/hd1/mouse/raw_reads/out_1.fastq /media/macmanes/hd1/mouse/raw_reads/out_2.fastq
	/home/macmanes/sga/src/SGA/sga index -a ropebwt -t 8 --no-reverse sga.fq
	/home/macmanes/sga/src/SGA/sga correct -k 25 --learn -t 8 --metrics=corr.mouse.metrics -o sga.mouse.corr.fq sga.fq > mouse.sga.out
	cat sga.mouse.corr.fq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/1$ > sga.mouse.corr.left.fq
	cat sga.mouse.corr.fq | grep -A3 --no-group-separator -E :[0-9]\{1,7\}:[0-9]\{1,7\}:[0-9]\{1,7\}/2$ > sga.mouse.corr.right.fq
	rm sga.fq sga.mouse.corr.fq sga.sai sga.bwt
else
	echo "Ive done mouse SGA"
fi
if [ ! -f sga.mouse.Trinity.fasta ] ; 
then
	rm -fr sga.mouse
	echo "Starting trinity"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left sga.mouse.corr.left.fq --right sga.mouse.corr.right.fq --full_cleanup --CPU 8 --output sga.mouse > mouse.sga.trin.out 2>&1
	gzip sga.mouse.corr.left.fq sga.mouse.corr.right.fq &
else
	echo "ive done the mouse SGA assembo"
fi
for i in `find * -type f -size -95c`; do rm $i; done
####Seecer
echo "##################################################"
echo "Starting MOUSE seecer @ `date`."
echo "##################################################"

cd ../seecer
if [ ! -f out_1.fastq_corrected.fa ] ; 
then
	sh ~/seecer/SEECER-0.1.2/SEECER/bin/run_seecer.sh -t /media/macmanes/hd/flux/seecer/ /media/macmanes/hd1/mouse/raw_reads/out_1.fastq /media/macmanes/hd1/mouse/raw_reads/out_2.fastq > mouse.seecer.out 2>&1
	mv /media/macmanes/hd1/mouse/raw_reads/out_2.fastq_corrected.fa .
	mv /media/macmanes/hd1/mouse/raw_reads/out_1.fastq_corrected.fa .
#	rm counts_17_3 corrected.fasta /media/macmanes/hd1/mouse/raw_reads/*N
else
	echo "mouse seecer done"
fi
for i in `find * -type f -size -95c`; do rm $i; done
if [ ! -f mouse.seecer.Trinity.fasta ] ; 
then
	rm -fr mouse.seecer
	echo "Starting trinity"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left out_1.fastq_corrected.fa --right out_2.fastq_corrected.fa --full_cleanup --CPU 8 --output mouse.seecer > mouse.seecer.trin.out
	gzip out_2.fastq_corrected.fa out_1.fastq_corrected.fa &
else
	echo "mouse seecer asse. done"
fi
for i in `find * -type f -size -95c`; do rm $i; done
####AllPaths #Need to do trinity assembly
echo "##################################################"
echo "Starting MOUSE AllPaths @ `date`."
echo "##################################################"

cd ../AllPaths
if [ ! -f mouse.allp.right.fasta ] ; 
then
	python ~/trinityrnaseq_r2013-02-25/util/ec_norm/preproc.py -m 30 -t 8 -l /media/macmanes/hd1/mouse/raw_reads/out_1.fastq -r /media/macmanes/hd1/mouse/raw_reads/out_2.fastq -o mouse.corr -H True --error_corr True
	sed -i 's_read_unpaired_g' corr.unpaired.fa
	cat corr.left.fa corr.unpaired.fa > mouse.allp.left.fasta
	mv corr.right.fa mouse.allp.right.fasta
#	rm corr.fastq.ids corr.left.fa corr.unpaired.fa
else
	echo " Allp mouse done"
fi
if [ ! -f mouse.allp.Trinity.fasta ] ; 
then
	rm -fr mouse.allp
	echo "Starting Mouse trinity"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left mouse.allp.left.fasta --right mouse.allp.right.fasta --full_cleanup --CPU 8 --output mouse.allp > mouse.allp.trin.out
else
	echo "mouse allp assem done"
fi
for i in `find * -type f -size -95c`; do rm $i; done
###Reptile
echo "##################################################"
echo "Starting MOUSE Reptile @ `date`."
echo "##################################################"

cd ../reptile

if [ ! -f right.mouse.reptile.corr.fa ] ; 
then
	echo "Mouse rept"
	cp /media/macmanes/hd1/mouse/raw_reads/out_1.fastq data/
	cp /media/macmanes/hd1/mouse/raw_reads/out_2.fastq data/
	perl ~/reptile-v1.1/reptile-v1.1/utils/fastq-converter-v2.0.pl data/ data/ 1 #Must have .fastQ extension
	~/reptile-omp/reptile-omp config.run.mouse
	~/reptile-v1.1/reptile-v1.1/utils/reptile_merger/reptile_merger data/out_2.fa data/out_2.err right.mouse.reptile.corr.fa
	~/reptile-omp/reptile-omp config.run1.mouse
	~/reptile-v1.1/reptile-v1.1/utils/reptile_merger/reptile_merger data/out_1.fa data/out_1.err left.mouse.reptile.corr.fa
	sed -i 's_^>[0-9]:.*_&/2_g' right.mouse.reptile.corr.fa
	sed -i 's_^>[0-9]:.*_&/1_g' left.mouse.reptile.corr.fa
else
	echo "I've already done the Rept corr"
fi
if [ ! -f mouse.reptile.Trinity.fasta ] ; 
then
	rm -fr reptile
	echo "Starting trinity"
	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fa --JM 30G --left left.mouse.reptile.corr.fa --right right.mouse.reptile.corr.fa --full_cleanup --CPU 8 --output mouse.reptile > mouse.reptile.trin.out
	gzip data/right.mouse.reptile.corr.fa data/left.mouse.reptile.corr.fa data/*fa
#	rm data/*q data/*err
else
	echo "I've already done the mouse Rept assembo"
fi
for i in `find * -type f -size -95c`; do rm $i; done


####Echo
echo "##################################################"
echo "Starting MOUSE Echo @ `date`."
echo "##################################################"

#cd ../echo
#if [ ! -f mouse.right.echo.fastq ] ; 
#then
#	python /home/macmanes/echo_v1_12/ErrorCorrection.py -b 2000000 --nh 2048 --ncpu 8  -o mouse.right.echo.fastq /media/macmanes/hd1/mouse/raw_reads/out_2.fastq
#	python /home/macmanes/echo_v1_12/ErrorCorrection.py -b 2000000 --nh 2048 --ncpu 8  -o mouse.left.echo.fastq /media/macmanes/hd1/mouse/raw_reads/out_1.fastq
#else
#	echo "echo mouse done"
#fi
#if [ ! -f mouse.echo.Trinity.fasta ] ; 
#then
#	rm -fr mouse.echo
#	echo "Starting trinity"
#	~/trinityrnaseq_r2013-02-25/Trinity.pl --seqType fq --JM 30G --left mouse.left.echo.fastq --right mouse.right.echo.fastq --full_cleanup --CPU 8 --output mouse.echo > mouse.echo.trin.out
#	gzip left.echo.fastq right.echo.fastq
#else
#	echo "mouse echo assembo done"
#fi
for i in `find * -type f -size -95c`; do rm $i; done

echo "##################################################"
echo "Starting SolexaQA @ `date`."
echo "##################################################"

cd ../solexaQA
perl /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1/SolexaQA.pl ../right.fq &
perl /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1/SolexaQA.pl ../left.fq
perl /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1/SolexaQA.pl /media/macmanes/hd1/mouse/raw_reads/out_1.fastq &
perl /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1/SolexaQA.pl /media/macmanes/hd1/mouse/raw_reads/out_2.fastq
for i in `find * -type f -size -95c`; do rm $i; done




