#!/usr/bin/env python3.5

'''Gets basic statistics on DNA sequences in all multi-fasta
files in directory with extension *.seq.'''

from os import listdir


__author__ = "Jeffrey P Tomkins, PhD"
__contributor__ = ""
__copyright__ = "Copyright 2016, Institute for Creation Research"
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Jeffrey P Tomkins"
__email__ = "jtomkins@icr.org"
__status__ = "Production"


def get_fasta(fastafile):
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


def get_seqfile_info(fasta_filename):
    '''Returns info from multi-fasta file in a 5-element
    tuple: (num_seqs, ave_seqlen, min_seqlen, max_seqlen,
    total_seqlen) utilizing the get_fasta function.'''

    fasta_tuple = get_fasta(fasta_filename)
    num_seqs = len(fasta_tuple)
    seqlen_list = [len(item[1]) for item in fasta_tuple]
    total_seqlen = sum(seqlen_list)
    ave_seqlen = total_seqlen // num_seqs
    min_seqlen = min(seqlen_list)
    max_seqlen = max(seqlen_list)

    return (num_seqs, ave_seqlen, min_seqlen, max_seqlen, total_seqlen)


def median(target_list):
    '''Returns the median of a list of numbers.'''
    n = len(target_list)
    sorted_list = (sorted(target_list))
    midpoint = n // 2

    if n % 2 == 1:
        # If odd return middle value
        return sorted_list[midpoint]
    else:
        # If even return average of middle values
        lo = midpoint - 1
        hi = midpoint
        return (sorted_list[lo] + sorted_list[hi]) / 2


def comma(num):
    '''Add comma to every 3rd digit. Takes integer or float and
    returns string.'''
    if type(num) == int:
        return '{:,}'.format(num)
    elif type(num) == float:
        return '{:,.2f}'.format(num)  # Rounds to 2 decimal places
    else:
        print("Need int or float as input to function comma()!")


def genome_coverage(seqlen, genome_size):
    '''Returns genome coverage using total sequence length (seqlen)
    and estimated genome size.'''
    return str(round(seqlen/genome_size, 2))


def get_all_stats(file_extension='seq'):
    '''Writes out a basic statistical summary of DNA sequences and
    prints the file to standard out.'''

    # Put file names with extension *seq in list.
    files = [file for file in listdir('.') if
             file.endswith(file_extension)]

    if not files:
        print("No files with extension", file_extension, "found!")

    # Get the basic seq stats and print to file
    with open('seq_num_report.txt', 'w') as fo:

        seqnum_list = []
        minseqlen_list = []
        maxseqlen_list = []
        seqlen_list = []

        # Write out intial header line
        fo.write("file_name     num     ave  min  max\n")

        # Process each seqfile and write out stats
        for file in files:
            # Get data for each file from 3-element tuple
            seq_data = get_seqfile_info(file)

            # Add commas and write out str data for each file
            num_str = comma(seq_data[0])
            aveseqlen_str = comma(seq_data[1])
            minlen_str = str(seq_data[2])
            maxlen_str = comma(seq_data[3])
            fo.write(file + "   " + num_str + "  " + aveseqlen_str +
                     "  " + minlen_str + "  " + maxlen_str + "\n")

            # Generate lists for summary of all files
            seqnum_list.append(seq_data[0])
            minseqlen_list.append(seq_data[2])
            maxseqlen_list.append(seq_data[3])
            seqlen_list.append(seq_data[4])

        overall_seqlen_mean = sum(seqlen_list) // sum(seqnum_list)

        # Write out summary of all seq files
        fo.write("======================================\n")
        fo.write("Summary for all files\n")
        fo.write("======================================\n")
        fo.write("Min file size: " + comma(min(seqnum_list)) + " seqs\n")
        fo.write("Max file size: " + comma(max(seqnum_list)) + " seqs\n")
        fo.write("Mean file size: " +
                 comma(sum(seqnum_list)//len(seqnum_list)) + " seqs\n")
        fo.write("Median file size: " + comma(median(seqnum_list)) +
                 " seqs\n")
        fo.write("Total seqs: " + comma(sum(seqnum_list)) + "\n")
        fo.write("Min seqlen: " + str(min(minseqlen_list)) + "\n")
        fo.write("Max seqlen: " + comma(max(maxseqlen_list)) + "\n")
        fo.write("Mean seqlen: " + comma(overall_seqlen_mean) + "\n")
        fo.write("Genome coverage: " +
                 genome_coverage(sum(seqlen_list), 3000000000) + "\n")

    # Print the output file to screen
    with open('seq_num_report.txt', 'r') as fi:
        for line in fi:
            print(line)


if __name__ == '__main__':
    get_all_stats()
