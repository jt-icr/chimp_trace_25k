#!/usr/bin/env python2.7
"""parse_multi_fasta.py
	Creator :     Jeff Tomkins <jtomkins@icr.org>
	Create date :  April 27, 2015 
	Description :  Selects 25000 random seqs and parses into new fasta file.
	Modifications : JP Tomkins, improved filename var usage, 11/16/2015
"""

from __future__ import print_function
import sys
import random

from Bio import SeqIO


args = sys.argv[1:]

if len(args) != 2 :
    print("Error!  2 args: fasta file and number of seqs for new fasta file")
    sys.exit()
    
in_filename = args[0]
num_records = int(args[1])
new_filename = in_filename.rsplit('.',1)[0] + "_" + str(num_records) + "_seqs.fa"

#Open files for read/write
fi = open(in_filename, 'rU')
fo = open(new_filename, 'w')

#Parse individual fasta records into a list
records = list(SeqIO.parse(fi, "fasta"))

#Only use sequences of 100 bases or more
records = [rec for rec in records if len(rec) >= 100]

#Make sure there are enough seqs
if len(records) <= 25000:
    print("Error! Need a fasta file of at least 25,000 seqs!")
    sys.exit()

#Randomize
random.shuffle(records)

#Get specified number of fasta records and print to new file
for i in records[0:num_records] :
	SeqIO.write(i, fo, "fasta")

#Print the new file name
print("New file: ", new_filename)

#Close files and exit
fi.close()
fo.close()
sys.exit()
