#!/usr/bin/python

from Bio import SeqIO
import subprocess
import sys
import os

short_qual = [] 
for record in SeqIO.parse(open("sim.1a.fastq", "rU"), "fastq") :
    print record.format("fastq") 
    #print record.letter_annotations["phred_quality"]
    if len(record.letter_annotations["phred_quality"]) < 76 :
        short_qual.append(record)
 
print "Found %i short quals" % len(short_qual)
 
output_handle = open("fixed_seqs.1a.fastQ", "w")
SeqIO.write(short_quals, output_handle, "fasta")
output_handle.close()
