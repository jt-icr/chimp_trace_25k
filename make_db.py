#!/Applications/anaconda3/bin/python3.5

'''Makes an sqlite3 database with three tables from a large
number of blast output files in csv format and their
corresponding fasta and xml (seq info) files. This module was
edited for python standards compliance using pep8.'''

from shutil import copyfileobj
import os

import sqlite3
import pandas as pd


__author__ = "Jeffrey P Tomkins, PhD"
__contributor__ = ""
__copyright__ = "Copyright 2016, Institute for Creation Research"
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey P Tomkins"
__email__ = "jtomkins@icr.org"
__status__ = "Production"


def add_seqfile_column(filename):
    '''Adds the name of the seqfile as 10th column in a csv file.'''

    seqfile_num = filename[4:7]
    fo = open(filename.rstrip(".csv") + "_.csv", "w")
    with open(filename, "r") as fi:
        for line in fi:
            fo.write(line.rstrip("\n") + "," + seqfile_num + "\n")
    fo.close()


def concat_files(concat_outfilename, file_extension):
    '''Concatenates files in directory - outfile and file extension
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
    '''Puts a fasta format file into a csv format file with two cols:
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


def create_db():
    '''Creates the db and it's tables.'''

    conn = sqlite3.connect('chimp_trace_25k.sqlite')
    c = conn.cursor()

    # Create tables
    c.execute('create table blast_data '
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

    c.execute('create table seq_data '
              '(gnl_num varchar(17) not null '
              ',dna_seq text not null '
              ',foreign key (gnl_num) references blast_data(qseqid))')

    c.execute('create table seq_year '
              '(seqfile_id varchar(3) not null '
              ',min_date varchar(4) '
              ',max_date varchar(4) '
              ',foreign key (seqfile_id) references blast_data(seqfile))')

    conn.commit()
    conn.close()


def csv_to_db(csvfile, db_table, names):

    conn = sqlite3.connect('chimp_trace_25k.sqlite')

    try:
        with open(csvfile, 'r') as fi:
            # Open *csv as a pandas data frame
            df = pd.read_csv(fi, header=None)
            df.columns = names
    except IOError:
        print('Cannot open', csvfile)

    df.to_sql(db_table, conn, flavor='sqlite',
              schema=None, if_exists='append', index=None,
              index_label=None, chunksize=None, dtype=None)

    conn.close()


if __name__ == '__main__':

    csvfiles = [file for file in os.listdir('.')
                if file.endswith('.csv')]

    if not csvfiles:
        print("No files with " + "*.csv" + " found!")

    for file in csvfiles:
        add_seqfile_column(file)

    # Concatenate all files in directory based on file extension
    concat_files('concat_blast_data.csv', '_.csv')
    concat_files('concat_fasta_data.fa', '.fa')

    # Put fasta data into csv format
    fasta_to_csv('concat_fasta_data.fa', 'concat_fasta_data.csv')

    # Make the db
    create_db()

    # Populate the blast_data table
    names1 = ['qseqid', 'qstart', 'qend', 'mismatch', 'gapopen',
              'pident', 'nident', 'length', 'qlen', 'seqfile']
    csv_to_db('concat_blast_data.csv', 'blast_data', names=names1)

    # Populate the seq_data table
    names2 = ['gnl_num', 'dna_seq']
    csv_to_db('concat_fasta_data.csv', 'seq_data', names=names2)

    # Populate the seq_year table from csv file created by
    # get_xml_seqyear.py script.
    names3 = ['seqfile_id', 'min_date', 'max_date']
    csv_to_db('seq_year_csv', 'seq_year', names=names3)

    # Clean up the tmp files
    for file in os.listdir():
        if file.endswith('_.csv') or file.startswith('concat_'):
            os.unlink(file)
