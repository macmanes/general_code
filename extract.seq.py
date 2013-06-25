#!/usr/bin/python

from Bio import SeqIO
import subprocess
import sys
import os


#this takes a list of numbers as a CSV file, 1 numbe per line. these numbers correspond >NUM in a multifasta file.

fasta_file = sys.argv[1] # Input fasta file
number_file = "list" # Input interesting numbers file, one per line
result_file = sys.argv[2] # Output fasta file



numbers = set()
with open(number_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            try:
                num = str(line)
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
            name = str(seq.id)
        except:
            print seq.id, "Fasta name cannot be converted to integer"
        if name in numbers:
            SeqIO.write([seq], f, "fasta")


