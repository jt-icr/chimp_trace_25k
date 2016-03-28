#!/usr/bin/env python3.5

'''Summarizes the results of a blastn output and can be run
as a batch script to analyze each of an entire set of blastn
output files in csv format. Also generates a summary of the
entire set of blastn output files.'''

from os import listdir
from time import strftime

import pandas as pd
from matplotlib import pyplot as plt


__author__ = "Jeffrey P Tomkins, PhD"
__contributor__ = ""
__copyright__ = "Copyright 2016, Institute for Creation Research"
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey P Tomkins"
__email__ = "jtomkins@icr.org"
__status__ = "Production"


date = strftime('%m_%d_%Y')

# Set final csv summary outfile name with date stamp
csv_outfilename = 'blastn_summary_' + date + '.dat'


def get_dataframe(csvfile):
    '''Puts a csv file into a pandas dataframe'''
    try:
        with open(csvfile, 'r') as fi:
            # Open *csv as a pandas data frame
            return pd.read_csv(fi)
    except IOError:
        print('Cannot open', csvfile)


def set_csv_header():
    '''Set the csv data file header line for blastn results.'''
    try:
        with open(csv_outfilename, 'w') as fo:
            header = ['file_id', 'aln_ident', 'qseq_ident', 'aln_len',
                      'qseqret', 'qseqall', 'num_qseqs', 'num_hits',
                      'hitfreq', 'overall_ident']
            fo.write(",".join(header) + "\n")
    except IOError:
            print('Cannot open', csv_outfilename)


def get_blast_data(csv_filename, fasta_filename):
    '''Takes a blastn output file in csv format along with it's
    corresponding fasta file and generates a set of basic stats
    for the blast results and the sequences used to generate them.'''

    def _get_fastadata(fastafile):
        '''Makes a list of tuples, each with a list containing fasta
        header data and a string containing DNA seq ([info], 'seq')'''
        try:
            with open(fastafile, 'r') as fi:
                return [(part[0],
                        part[2].replace('\n', ''))
                        for part in
                        [entry.partition('\n')
                        for entry in fi.read().split('>')[1:]]]
        except IOError:
            print('Cannot open', fasta_file)

    # Get and setup the pandas dataframe
    df = get_dataframe(csv_filename)
    seqdata = _get_fastadata(fasta_filename)
    file_id = "pan" + csv_filename[4:7]

    num_hits = df.count()[0]
    num_qseqs = len(seqdata)

    # Assign column headings - replaces default [0,1,2...etc]
    df.columns = \
        'qseqid qstart qend mismatch gapopen pident nident length qlen'.split()

    # Speeds calc and removes this col as an obj dtype
    del df['qseqid']

    # Get overall qseq aln identity and put in new column called 'ident'
    df['ident'] = df.apply(lambda row: (row['nident']/row['length']
                           if row['length'] > row['qlen']
                           else row['nident']/row['qlen']), axis=1)

    # Create variables for outputs
    total_aligned = (df['ident'].mean() * df['qlen'].mean()) * num_hits
    total_seq_len = sum(len(s[1]) for s in (seqdata))
    ave_seqlen = ((total_seq_len/num_qseqs) * 100, 2)
    overall_ident = (total_aligned/total_seq_len) * 100
    ave_aln_ident = df['pident'].mean()
    ave_qseq_ident = df['ident'].mean() * 100
    ave_aln_len = df['length'].mean()
    ave_qseqret = df['qlen'].mean()
    ave_qseqall = total_seq_len/num_qseqs
    ave_hitfreq = num_hits/num_qseqs * 100

    # Send a human readable report to standard out
    print("Ave aln ident   : ", round(ave_aln_ident, 2))
    print("Ave qseq ident  : ", round(ave_qseq_ident, 2))
    print("Ave aln len     : ", round(ave_aln_len, 2))
    print("Ave qseqret len : ", round(ave_qseqret, 2))
    print("Ave qseqall len : ", round(ave_qseqall, 2))
    print("Num queryseqs   : ", num_qseqs)
    print("Num queryhits   : ", num_hits)
    print("Ave hit freq    : ", round(ave_hitfreq, 2))
    print("Overall ident   : ", round(overall_ident, 2))

    # Append csv data rows to header line
    try:
        with open(csv_outfilename, 'a') as fo:
            outline = [file_id, ave_aln_ident, ave_qseq_ident,
                       ave_aln_len, ave_qseqret, ave_qseqall,
                       num_qseqs, num_hits, ave_hitfreq,
                       overall_ident]
            fo.write(",".join(str(item) for item in outline) + "\n")
    except IOError:
        print('Cannot open', csv_outfilename)


def get_finalstats(csv_filename):
    '''Gets the final summary stats of csv file produced from the
    get_blast_data function. Header line of csv file: file_id,
    aln_ident,qseq_ident,aln_len,qseqret,qseqall,num_qseqs,num_hits,
    hitfreq,overall_ident'''

    df = get_dataframe(csv_filename)

    # Create two new dataframes based on overall ident
    hi = df[df['overall_ident'] > 89.99]
    lo = df[df['overall_ident'] < 90.00]

    # Get the statistics for all of the datasets
    print("<Summary stats for all datasets>")
    print("Ave aln ident   : ", round(df['aln_ident'].mean(), 2))
    print("Ave qseq ident  : ", round(df['qseq_ident'].mean(), 2))
    print("Ave aln len     : ", round(df['aln_len'].mean(), 2))
    print("Ave qseqret len : ", round(df['qseqret'].mean(), 2))
    print("Ave qseqall len : ", round(df['qseqall'].mean(), 2))
    print("Num queryseqs   : ", round(df['num_qseqs'].mean(), 0))
    print("Num queryhits   : ", round(df['num_hits'].mean(), 0))
    print("Ave hit freq    : ", round(df['hitfreq'].mean(), 2))
    print("Overall ident   : ", round(df['overall_ident'].mean(), 2))
    print()

    # Get the number of data sets grouped by overall ident
    print("<Num entries by overall ident>")
    print("Over  89.99% :", hi.count()[0])
    print("Below 90.00% :", lo.count()[0])
    print()
    # Print the data to file
    hi_file = "high_overallident_" + date + ".dat"
    hi.to_csv(hi_file)
    lo_file = "low_overallident_" + date + ".dat"
    lo.to_csv(lo_file)

    # Get stats for high overall identity data
    print("<Summary stats for high identity data>")
    print("Ave aln ident   : ", round(hi['aln_ident'].mean(), 2))
    print("Ave qseq ident  : ", round(hi['qseq_ident'].mean(), 2))
    print("Ave aln len     : ", round(hi['aln_len'].mean(), 2))
    print("Ave qseqret len : ", round(hi['qseqret'].mean(), 2))
    print("Ave qseqall len : ", round(hi['qseqall'].mean(), 2))
    print("Num queryseqs   : ", round(hi['num_qseqs'].mean(), 0))
    print("Num queryhits   : ", round(hi['num_hits'].mean(), 0))
    print("Ave hit freq    : ", round(hi['hitfreq'].mean(), 2))
    print("Overall ident   : ", round(hi['overall_ident'].mean(), 2))
    print()

    # Get stats for low overall identity data
    print("<Summary stats for low identity data>")
    print("Ave aln ident   : ", round(lo['aln_ident'].mean(), 2))
    print("Ave qseq ident  : ", round(lo['qseq_ident'].mean(), 2))
    print("Ave aln len     : ", round(lo['aln_len'].mean(), 2))
    print("Ave qseqret len : ", round(lo['qseqret'].mean(), 2))
    print("Ave qseqall len : ", round(lo['qseqall'].mean(), 2))
    print("Num queryseqs   : ", round(lo['num_qseqs'].mean(), 0))
    print("Num queryhits   : ", round(lo['num_hits'].mean(), 0))
    print("Ave hit freq    : ", round(lo['hitfreq'].mean(), 2))
    print("Overall ident   : ", round(lo['overall_ident'].mean(), 2))
    print()

    # Line graph the overall ident of the fasta files
    fasta_files = df['file_id']
    fasta_nums = [int(item[3:7]) for item in fasta_files]
    ident_dat = df['overall_ident']
    plt.plot(fasta_nums, ident_dat, color='green',
             marker='o', linestyle='solid')
    plt.title("Overall Identity of Each Fasta Trace Read Data Set")
    plt.ylabel("Overall Data Set Percent Identity")
    plt.xlabel("Fasta Trace Read Data Sets (001 to 101)")
    plt.xlim(0, 102)
    plt.show()


def blastn_batch_proc():
    '''Run this module as a batch script on all blastn csv and fasta
    files in the working directory.'''

    set_csv_header()

    csv_files = sorted([file for file in listdir('.')
                        if file.endswith('csv')])
    if not csv_files:
        print("No files with *.csv found!")

    fasta_files = sorted([file for file in listdir('.')
                          if file.endswith('fa')])
    if not fasta_files:
        print("No files with *.fa found!")

    for csv_filename, fa_filename in zip(csv_files, fasta_files):
        print("<{}>".format(csv_filename))
        print("<{}>".format(fa_filename))
        get_blast_data(csv_filename, fa_filename)
        print()

    get_finalstats(csv_outfilename)


if __name__ == '__main__':
    blastn_batch_proc()
