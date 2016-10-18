#!/usr/bin/env python3.5

'''This module contains functions helpful in creating an SQL
database from fasta files and blastn data. This module was
edited for python standards compliance using pep8.'''

from shutil import copyfileobj
from shutil import move
import os
import sqlite3

import pandas as pd

__author__ = "Jeffrey P Tomkins, PhD"
__copyright__ = "Copyright 2016, Institute for Creation Research"
__email__ = "jtomkins@icr.org"


def add_seqfile_column(filename):
    '''Adds the name of the seqfile as 10th column in a csv file.'''

    seqfile_num = filename[4:7]
    with open(filename, "r") as fi:
        with open(filename.rstrip(".csv") + "_.csv", "w") as fo:
            for line in fi:
                fo.write(line.rstrip("\n") + "," + seqfile_num + "\n")


def concat_files(concat_outfilename, file_extension):
    '''Concatenates files in directory - outfile name and file extension
    of target files to concat as args.'''

    files = sorted([file for file in os.listdir('.')
                    if file.endswith(file_extension)])

    if not files:
        print("No files with " + file_extension + " found!")

    try:
        fo = open(concat_outfilename, 'w')
    except IOError:
        print('Cannot open', concat_outfilename)

    for file in files:
        try:
            fi = open(file, 'r')
            copyfileobj(fi, fo)
            fi.close()
        except IOError:
            print('Cannot open', file)
    fo.close()


def fasta_to_csv(fasta_filename, out_csvfilename):
    '''Puts a fasta format file into a csv file with two cols:
    the first being the seq id and then the DNA seq.'''

    def _get_fastadata(fastafile):
        '''Makes a list of tuples, each with a list containing fasta
        header data and a string containing DNA seq ([info], 'seq').'''
        try:
            with open(fastafile, 'r') as fi:
                return [(part[0],
                        part[2].replace('\n', ''))
                        for part in
                        [entry.partition('\n')
                        for entry in fi.read().split('>')[1:]]]
        except IOError:
            print('Cannot open', fasta_file)

    try:
        with open(out_csvfilename, 'w') as fo:
            for entry in _get_fastadata(fasta_filename):
                fo.write(entry[0].split()[0] + "," + entry[1] + "\n")
    except IOError:
        print('Cannot open', out_csvfilename)


def csv_to_db(db_name, csvfile, db_table, names):
    '''Populates an sqlite table from a csv file using the pandas library.'''
    conn = sqlite3.connect(db_name)

    df = pd.read_csv(csvfile, header=None)
    df.columns = names
    df.to_sql(db_table, conn, flavor='sqlite',
              schema=None, if_exists='append', index=None,
              index_label=None, chunksize=None, dtype=None)

    conn.close()


# CSV file prep and development of the chimp_trace_25k database
if __name__ == '__main__':

    # Target directory constants and db name
    DIR1 = '/Users/jtomkins/data/chimp_trace/25k_test'
    DIR2 = '/Users/jtomkins/data/chimp_trace/25k_test/non_hitters/on_chimp'
    DB_NAME = 'chimp_trace_25k.sqlite'

    # Directory safety check - should be in DIR1
    if os.getcwd() != DIR1:
        print('Not in proper directory!')
        exit()

    # Create the chimp_trace_25k db and it's tables.'''
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()

    c.execute('create table chimp_seq_data '
              '(gnl_num varchar(17) not null '
              ',dna_seq text not null '
              ',primary key (gnl_num))')

    c.execute('create table chimp_seq_year '
              '(seqfile_id varchar(3) not null '
              ',min_date varchar(4) '
              ',max_date varchar(4) '
              ',primary key (seqfile_id))')

    c.execute('create table chimp_blast_on_homo '
              '(qseqid varchar(17) not null '
              ',qstart integer not null '
              ',qend integer not null '
              ',mismatch integer not null '
              ',gapopen integer not null '
              ',pident decimal(4,2) not null '
              ',nident integer not null '
              ',length integer not null '
              ',qlen integer not null '
              ',seqfile varchar(3) not null '
              ',primary key (qseqid))')

    c.execute('create table nonhitter_blast_on_chimp '
              '(qseqid varchar(17) not null '
              ',qstart integer not null '
              ',qend integer not null '
              ',mismatch integer not null '
              ',gapopen integer not null '
              ',pident decimal(4,2) not null '
              ',nident integer not null '
              ',length integer not null '
              ',qlen integer not null '
              ',seqfile varchar(3) not null '
              ',primary key (qseqid))')

    c.execute('create index gnl_num_idx on chimp_seq_data (gnl_num)')
    c.execute('create index qseqid_idx_on_homo on chimp_blast_on_homo (qseqid)')
    c.execute('create index qseqid_idx_on_pan on nonhitter_blast_on_chimp (qseqid)')

    conn.commit()
    conn.close()

    # Add the seqfile name as 10th column in csv files
    csvfiles1 = [file for file in os.listdir(DIR1)
                 if file.endswith('.csv')]
    csvfiles2 = [file for file in os.listdir(DIR2)
                 if file.endswith('.csv')]

    if not csvfiles1:
        print("No files with *.csv found! in", DIR1)
        
    if not csvfiles2:
        print("No files with *.csv found! in", DIR2)

    os.chdir(DIR1)
    for file in csvfiles1:
        add_seqfile_column(file)

    os.chdir(DIR2)
    for file in csvfiles2:
        add_seqfile_column(file)

    # Process files in target directories based on file extension
    os.chdir(DIR1)
    concat_files('concat_chimp_blast_on_homo.csv', '_.csv')
    concat_files('concat_fasta_data.fa', '.fa')
    fasta_to_csv('concat_fasta_data.fa', 'concat_fasta_data.csv')
    os.chdir(DIR2)
    concat_files('concat_nonhitter_blastdat.csv', '_.csv')

    # Populate the chimp_blast_on_homo table
    os.chdir(DIR1)
    names1 = ['qseqid', 'qstart', 'qend', 'mismatch', 'gapopen',
              'pident', 'nident', 'length', 'qlen', 'seqfile']
    csv_to_db(DB_NAME, 'concat_chimp_blast_on_homo.csv', 'chimp_blast_on_homo',
              names=names1)

    # Populate the chimp_seq_data table
    names2 = ['gnl_num', 'dna_seq']
    csv_to_db(DB_NAME, 'concat_fasta_data.csv', 'chimp_seq_data', names=names2)

    # Populate the chimp_seq_year table from csv file created by
    # get_xml_seqyear.py script.
    names3 = ['seqfile_id', 'min_date', 'max_date']
    csv_to_db(DB_NAME, 'seq_year_csv', 'chimp_seq_year', names=names3)

    # Populate the nonhitter_blast_on_chimp table from csv file created by
    move(DIR2 + '/concat_nonhitter_blastdat.csv', DIR1 + '/concat_nonhitter_blastdat.csv')
    names4 = ['qseqid', 'qstart', 'qend', 'mismatch', 'gapopen',
              'pident', 'nident', 'length', 'qlen', 'seqfile']
    csv_to_db(DB_NAME, 'concat_nonhitter_blastdat.csv',
            'nonhitter_blast_on_chimp', names=names4)

    # Clean up the tmp files
    for file in os.listdir():
        if file.endswith('_.csv') or file.startswith('concat_'):
            os.unlink(file)
            
    exit()
