#!/usr/bin/env python3.5
'''Extracts the minmum and maximum run years from the xml
sequencing info file and puts into a csv file with seqfile
name. If xml file contains no year tags, 'NULL' is used.
e.g.
001,2000,2002
005,2002,2002
012,NULL,NULL
022,2002,2002'''

from os import listdir
import re


__author__ = "Jeffrey P Tomkins, PhD"
__contributor__ = ""
__copyright__ = "Copyright 2016, Institute for Creation Research"
__email__ = "jtomkins@icr.org"


# Put all the DNA seq xml files in a list
xml_files = sorted([file for file in listdir('.')
                    if file.startswith('xml')])

regex_yr = re.compile(r'(\d{4})')

seq_file = []
run_years = []

# Populate the seq_file and run_years lists from each xml file
for filename in xml_files:
    seq_file.append(filename[-3:])
    years = []
    with open(filename, 'r') as fi:
    	for line in fi:
            if '<RUN_DATE>' in line:
                year = re.search(regex_yr, line).group(0)       
                if year not in years:
                    years.append(year)
    if not years:
    	run_years.append(['NULL', 'NULL'])
    else:
        run_years.append(years)

# Write out the csv file as depicted above
with open("seq_year.csv", "w") as fo:
	for seqfile, year in zip(seq_file, run_years):
	    fo.write(seqfile + "," + min(year) + "," + max(year) + "\n")
