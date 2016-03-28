#!/usr/bin/perl -w
#Produces a set of shell scripts for blasting 
#JP Tomkins, Nov 17, 2015

my @query_files = ();

#Generate an array of query filenames incrementing the numerical portion of the string.
for (my $i = 1; $i <= 101; $i += 1) {
    if ($i < 10) {
        push @query_files, "pan_00".$i."_25000_seqs\.fa";
    }
    elsif ($i < 100) {
        push @query_files, "pan_0".$i."_25000_seqs\.fa";
    }
    elsif($i < 102) {
        push @query_files, "pan_".$i."_25000_seqs\.fa";
    }      
}

foreach my $query_file (@query_files) {
    my $outfile = substr($query_file, 0, 7) .".sh";
	unless (open(OUT, ">$outfile")) {
        print "Error: Cannot open file \"$outfile\" to write to!";
    }
    my $csv = substr($query_file, 0, 7) ."_25k_on_homo.csv";

print OUT"#!/bin/sh
blastn \\
-query $query_file \\
-task blastn \\
-db /home/jtomkins/homo_genome/homo_genome_GRCh37.71.fa \\
-out $csv \\
-evalue 0.1 \\
-word_size 11 \\
-outfmt \"10 qseqid qstart qend mismatch gapopen pident nident length qlen\" \\
-max_target_seqs 1 \\
-max_hsps 1 \\
-dust no \\
-soft_masking false \\
-perc_identity 50 \\
-gapopen 3 \\
-gapextend 3 \\
-num_threads 10";

chmod 0755, $outfile;
}
