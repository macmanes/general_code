#!/usr/bin/python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import subprocess
import sys
import os


#this takes a list of numbers as a CSV file, 1 numbe per line. these numbers correspond >NUM in a multifasta file.

fasta_file = "/media/hd2/sparrow/abyssk30-long.contigs.fa" # Input fasta file
number_file = "/media/hd2/sparrow/numbers.txt.csv" # Input interesting numbers file, one per line
result_file = "/media/hd2/sparrow/sparrow.only.bird.fa" # Output fasta file
output_xml = "/media/hd2/sparrow/bird.blast.txt"

numbers = set()
with open(number_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            try:
                num = int(line)
                numbers.add(num)
            except:
                print line, "Line cannot be converted to integer"

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
end = False
with open(result_file, "w") as f:
    while end != True:
        try:
            seq = fasta_sequences.next()
        except:
            end = True
        try:
            name = int(seq.id)
        except:
            print seq.id, "Fasta name cannot be converted to integer"
        if name in numbers:
            SeqIO.write([seq], f, "fasta")



DATABASE='/media/hd/blastdb/finch.fa'
E='1e-10'

        #create the commandline string


cl = NcbiblastnCommandline(cmd='/home/matthew/ncbi-blast/bin/blastn', query=result_file, db=DATABASE, evalue=E, outfmt=6, out=output_xml, num_alignments=1, num_descriptions=1, num_threads=8)

        #actually run BLAST

stdout, stderr = cl()

for record in NCBIXML.parse(open(output_xml)) :
    print "QUERY: %s" % record.query
    for align in record.alignments :
        print " MATCH: %s" % align.title[:160]
        #for hsp in align.hsps :
            #print " HSP, e=%f, from position %i to %i" \
                #% (hsp.expect, hsp.query_start, hsp.query_end)
            #if hsp.align_length < 60 :
                 #print "  Query: %s" % hsp.query
                 #print "  Match: %s" % hsp.match
                 #print "  Sbjct: %s" % hsp.sbjct
            #else :
                 #print "  Query: %s..." % hsp.query[:57]
                 #print "  Match: %s..." % hsp.match[:57]
                 #print "  Sbjct: %s..." % hsp.sbjct[:57]



