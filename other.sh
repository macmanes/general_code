#!/bin/bash



 
java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q0/celegans.Left.Q0.pp.fq.gz \
q0/celegans.Left.Q0.up.fq.gz \
q0/celegans.Right.Q0.pp.fq.gz \
q0/celegans.Right.Q0.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:0 \
TRAILING:0 \
SLIDINGWINDOW:8:0 \
MINLEN:50; 


java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q2/celegans.Left.Q2.pp.fq.gz \
q2/celegans.Left.Q2.up.fq.gz \
q2/celegans.Right.Q2.pp.fq.gz \
q2/celegans.Right.Q2.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:2 \
TRAILING:2 \
SLIDINGWINDOW:8:2 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q5/celegans.Left.Q5.pp.fq.gz \
q5/celegans.Left.Q5.up.fq.gz \
q5/celegans.Right.Q5.pp.fq.gz \
q5/celegans.Right.Q5.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:8:5 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q10/celegans.Left.Q10.pp.fq.gz \
q10/celegans.Left.Q10.up.fq.gz \
q10/celegans.Right.Q10.pp.fq.gz \
q10/celegans.Right.Q10.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:10 \
TRAILING:10 \
SLIDINGWINDOW:8:10 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q15/celegans.Left.Q15.pp.fq.gz \
q15/celegans.Left.Q15.up.fq.gz \
q15/celegans.Right.Q15.pp.fq.gz \
q15/celegans.Right.Q15.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:15 \
TRAILING:15 \
SLIDINGWINDOW:8:15 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR554453_1.fastq \
SRR554453_2.fastq
q20/celegans.Left.Q20.pp.fq.gz \
q20/celegans.Left.Q20.up.fq.gz \
q20/celegans.Right.Q20.pp.fq.gz \
q20/celegans.Right.Q20.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:8:20 \
MINLEN:50; 


cd q0
cat celegans.Left.Q0.up.fq.gz celegans.Right.Q0.up.fq.gz > celegans.Q0.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q0.pp.fq.gz -r celegans.Right.Q0.pp.fq.gz -u celegans.Q0.up.fq.gz -o q0 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q0 --bfly_opts --REDUCE


cd ../q10
cat celegans.Left.Q10.up.fq.gz celegans.Right.Q10.up.fq.gz > celegans.Q10.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q10.pp.fq.gz -r celegans.Right.Q10.pp.fq.gz -u celegans.Q10.up.fq.gz -o q10 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q10 --bfly_opts --REDUCE

cd ../q15
#cat celegans.Left.Q15.up.fq.gz celegans.Right.Q15.up.fq.gz > celegans.Q15.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q15.pp.fq.gz -r celegans.Right.Q15.pp.fq.gz -u celegans.Q15.up.fq.gz -o q15 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q15 --bfly_opts --REDUCE

cd ../q20
cat celegans.Left.Q20.up.fq.gz celegans.Right.Q20.up.fq.gz > celegans.Q20.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q20.pp.fq.gz -r celegans.Right.Q20.pp.fq.gz -u celegans.Q20.up.fq.gz -o q20 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q20 --bfly_opts --REDUCE

cd ../q2
cat celegans.Left.Q2.up.fq.gz celegans.Right.Q2.up.fq.gz > celegans.Q2.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q2.pp.fq.gz -r celegans.Right.Q2.pp.fq.gz -u celegans.Q2.up.fq.gz -o q2 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q2 --bfly_opts --REDUCE

cd ../q5
cat celegans.Left.Q5.up.fq.gz celegans.Right.Q5.up.fq.gz > celegans.Q5.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l celegans.Left.Q5.pp.fq.gz -r celegans.Right.Q5.pp.fq.gz -u celegans.Q5.up.fq.gz -o q5 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q5 --bfly_opts --REDUCE



cd /media/macmanes/hd/mouse/homo
 
java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q0/homo.Left.Q0.pp.fq.gz \
q0/homo.Left.Q0.up.fq.gz \
q0/homo.Right.Q0.pp.fq.gz \
q0/homo.Right.Q0.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:0 \
TRAILING:0 \
SLIDINGWINDOW:8:0 \
MINLEN:50; 


java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q2/homo.Left.Q2.pp.fq.gz \
q2/homo.Left.Q2.up.fq.gz \
q2/homo.Right.Q2.pp.fq.gz \
q2/homo.Right.Q2.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:2 \
TRAILING:2 \
SLIDINGWINDOW:8:2 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q5/homo.Left.Q5.pp.fq.gz \
q5/homo.Left.Q5.up.fq.gz \
q5/homo.Right.Q5.pp.fq.gz \
q5/homo.Right.Q5.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:5 \
TRAILING:5 \
SLIDINGWINDOW:8:5 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q10/homo.Left.Q10.pp.fq.gz \
q10/homo.Left.Q10.up.fq.gz \
q10/homo.Right.Q10.pp.fq.gz \
q10/homo.Right.Q10.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:10 \
TRAILING:10 \
SLIDINGWINDOW:8:10 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q15/homo.Left.Q15.pp.fq.gz \
q15/homo.Left.Q15.up.fq.gz \
q15/homo.Right.Q15.pp.fq.gz \
q15/homo.Right.Q15.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:15 \
TRAILING:15 \
SLIDINGWINDOW:8:15 \
MINLEN:50; 

java -Xmx10g -classpath /home/macmanes/software/trimmomatic-0.22.jar \
org.usadellab.trimmomatic.TrimmomaticPE \
-phred33 -threads 8 \
SRR597900_1.fastq \
SRR597900_2.fastq
q20/homo.Left.Q20.pp.fq.gz \
q20/homo.Left.Q20.up.fq.gz \
q20/homo.Right.Q20.pp.fq.gz \
q20/homo.Right.Q20.up.fq.gz \
ILLUMINACLIP:barcodes.fa:2:40:15 \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:8:20 \
MINLEN:50; 


cd q0
cat homo.Left.Q0.up.fq.gz homo.Right.Q0.up.fq.gz > homo.Q0.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q0.pp.fq.gz -r homo.Right.Q0.pp.fq.gz -u homo.Q0.up.fq.gz -o q0 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q0 --bfly_opts --REDUCE


cd ../q10
cat homo.Left.Q10.up.fq.gz homo.Right.Q10.up.fq.gz > homo.Q10.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q10.pp.fq.gz -r homo.Right.Q10.pp.fq.gz -u homo.Q10.up.fq.gz -o q10 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q10 --bfly_opts --REDUCE

cd ../q15
#cat homo.Left.Q15.up.fq.gz homo.Right.Q15.up.fq.gz > homo.Q15.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q15.pp.fq.gz -r homo.Right.Q15.pp.fq.gz -u homo.Q15.up.fq.gz -o q15 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q15 --bfly_opts --REDUCE

cd ../q20
cat homo.Left.Q20.up.fq.gz homo.Right.Q20.up.fq.gz > homo.Q20.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q20.pp.fq.gz -r homo.Right.Q20.pp.fq.gz -u homo.Q20.up.fq.gz -o q20 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q20 --bfly_opts --REDUCE

cd ../q2
cat homo.Left.Q2.up.fq.gz homo.Right.Q2.up.fq.gz > homo.Q2.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q2.pp.fq.gz -r homo.Right.Q2.pp.fq.gz -u homo.Q2.up.fq.gz -o q2 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q2 --bfly_opts --REDUCE

cd ../q5
cat homo.Left.Q5.up.fq.gz homo.Right.Q5.up.fq.gz > homo.Q5.up.fq.gz
python ~/git/trinityrnaseq/util/ec_norm/preproc.py -m30 -p33 -t8 -l homo.Left.Q5.pp.fq.gz -r homo.Right.Q5.pp.fq.gz -u homo.Q5.up.fq.gz -o q5 -H True --full True
~/trinityseq/trunk/Trinity.pl --seqType fa --SS_lib_type RF \
--JM 30G --left norm.corr.left.fa --right norm.corr.right.fa \
--CPU 8 --full_cleanup --output q5 --bfly_opts --REDUCE

























