#!/usr/bin/env python3.5

from os import listdir


__author__ = "Jeffrey P Tomkins, PhD"
__copyright__ = "Copyright 2016, Institute for Creation Research"
__email__ = "jtomkins@icr.org"


def get_non_hitters(in_csvfile, in_fastafile):
	'''Takes fasta file and corresponding blastn output file
	in csv format as input. Puts all non hitting sequences in 
	a new fasta file.'''

    fi_fasta = open(in_fastafile, 'r')
    fi_csv = open(in_csvfile, 'r')
    fo = open(in_fastafile.split('.')[0] + "_non_hitters.fa", 'w')

    # Get id nums and put in a set
    id_list = []
    for line in fi_csv:
        id_line = line.split(',')
        id_list.append(id_line[0])

    # Step through input fasta file and check ID against set
    # Write out fasta seqs to file not in set
    id_fset = frozenset(id_list)
    for line in fi_fasta:
        if line.startswith('>'):
            id_num = str(line.lstrip('>').split(' ')[0])
            if id_num not in id_fset:
                fo.write(line)
                reject_flag = False
            else:
                reject_flag = True
        else:
            if not reject_flag:
                fo.write(line)

    fi_fasta.close()
    fi_csv.close()
    fo.close()


def file_batch():
    '''Run this module as a batch script on all blastn csv and 
    fasta files in the working directory.'''

    csv_files = sorted([file for file in listdir('.')
                        if file.endswith(r'.csv')])
    if not csv_files:
        print("No files with *.csv found!")

    fasta_files = sorted([file for file in listdir('.')
                          if file.endswith(r'.fa')])
    if not fasta_files:
        print("No files with *.fa found!")

    print("Creating non-hitter fasta files for...")

    for csv_file, fa_file in zip(csv_files, fasta_files):
        get_non_hitters(csv_file, fa_file)
        print(csv_file, "and", fa_file)


if __name__ == '__main__':
    file_batch()
