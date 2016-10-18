--SQL queries on the chimp_trace_25k db
--=====================================

--Check tables
select * from chimp_blast_on_homo limit 5;

select * from chimp_seq_data limit 5;

select * from chimp_seq_year limit 5;

select * from nonhitter_blast_on_chimp limit 5;

--Table row counts
select count(*)
from chimp_blast_on_homo;

select count(*)
from chimp_seq_data;

select count(*)
from chimp_seq_year;

select count(*)
from nonhitter_blast_on_chimp;

--====================================
--Analysis of all seq data sets
--====================================

--Get averages for all seqfiles
select seqfile
    ,round(avg(pident),1) as avg_pident
    ,min(pident) as min_pident
    ,round(avg(gapopen),1) as avg_gapopen
    ,round(avg(mismatch),1) as avg_mismatch
    ,round(avg(length),1) as avg_aln_len
    ,round(avg(qlen),1) as avg_qlen
    ,round(avg(nident),1) as avg_nident
    ,(sum(nident) / ((25000*101) * avg(qlen)))*100 as overall_ident
from chimp_blast_on_homo;

--Get averages for each seqfile
select seqfile
    ,round(avg(pident),1) as avg_pident
    ,min(pident) as min_pident
    ,round(avg(gapopen),1) as avg_gapopen
    ,round(avg(mismatch),1) as avg_mismatch
    ,round(avg(qlen),1) as avg_qlen
    ,round(avg(nident),1) as avg_nident
    ,(sum(nident) / (25000 * avg(qlen)))*100 as overall_ident
from chimp_blast_on_homo
group by seqfile
order by seqfile + 0 asc;

--Create a view with number of sequences
--where aln greater than 99.7
create view num_seqs_99_8 as
select chimp_blast_on_homo.seqfile, count(chimp_blast_on_homo.seqfile) as num_seqs_99_8
from chimp_blast_on_homo
where chimp_blast_on_homo.pident > 99.7
group by chimp_blast_on_homo.seqfile
order by chimp_blast_on_homo.seqfile + 0 asc;

drop view num_seqs_99_8;
select * from num_seqs_99_8;

--Percent hits and overall ident per file
create view perc_hits_overall_ident as
select chimp_blast_on_homo.seqfile as seq_file
    ,count(*)*100/25000 as perc_hits
    ,(sum(nident) / (25000 * avg(qlen)))*100  as overall_ident
from chimp_blast_on_homo
group by chimp_blast_on_homo.seqfile
order by chimp_blast_on_homo.seqfile + 0 asc;

drop view perc_hits_overall_ident;
select * from perc_hits_overall_ident;

--Num seqfiles < 90.5% overall ident
select count(*) from perc_hits_overall_ident
where round(overall_ident,2) < 90.5;

--Number of perc_hits with num seqs greater than 99.7% ident
--Join two views above
select seqfile, perc_hits, num_seqs_99_8
from num_seqs_99_8 join perc_hits_overall_ident
  on seqfile = seq_file
order by seqfile + 0 asc;

--Get avg aln data in specific ranges
SELECT count(*)
  ,round(avg(pident),1) as avg_pident
  ,round(avg(gapopen),1) as avg_gapopen
  ,round(avg(mismatch),1) as avg_mismatch
  ,round(avg(length),1) as avg_aln_len
from chimp_blast_on_homo
where pident between 90.00 and 100.00;

SELECT count(*)
  ,round(avg(pident),1) as avg_pident
  ,round(avg(gapopen),1) as avg_gapopen
  ,round(avg(mismatch),1) as avg_mismatch
  ,round(avg(length),1) as avg_aln_len
from chimp_blast_on_homo
where pident between 80.00 and 90.00;

SELECT count(*)
  ,round(avg(pident),1) as avg_pident
  ,round(avg(gapopen),1) as avg_gapopen
  ,round(avg(mismatch),1) as avg_mismatch
  ,round(avg(length),1) as avg_aln_len
from chimp_blast_on_homo
where pident between 70.00 and 80.00;

SELECT count(*)
  ,round(avg(pident),1) as avg_pident
  ,round(avg(gapopen),1) as avg_gapopen
  ,round(avg(mismatch),1) as avg_mismatch
  ,round(avg(length),1) as avg_aln_len
from chimp_blast_on_homo
where pident between 60.00 and 70.00;

--====================================
--Analysis with seq year data
--====================================

--Number of seq sets with year info
select count(*)
from chimp_seq_year
where min_date not null and
      max_date not null;

--Number of hits returned per file containing year data
create view perc_hits_overall_ident_yr as
  select seqfile
  ,avg(pident) as pident
  ,avg(gapopen) as gapopen
  ,avg(mismatch) as mismatch
  ,avg(length) as aln_len
  ,avg(qlen) as qlen
  ,avg(nident)/avg(qlen) as qident
  ,(count(*)*100/25000) as perc_hits
  ,(sum(nident) / (25000 * avg(qlen)))*100  as overall_ident
  ,cast(min_date as int) as begin_yr
  ,cast(max_date as int) as end_yr
from chimp_blast_on_homo join chimp_seq_year
  on seqfile = seqfile_id
where min_date not null and
      max_date not null
group by seqfile
order by seqfile + 0 asc;

--Check perc_hits_yr view
select * from perc_hits_overall_ident_yr;

drop view perc_hits_overall_ident_yr;

--Get averages prior to 2005
select 
      round(avg(pident), 2) as avg_pident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(aln_len),1) as avg_aln_len
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_yr
where end_yr < 2005;

--Get averages after 2004
select 
      round(avg(pident), 2) as avg_pident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(aln_len),1) as avg_aln_len
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_yr
where end_yr > 2004;

--Get files with perc hits > 95
select *
from perc_hits_overall_ident_yr
where perc_hits > 95;

--Get average perc_hits for those > 95
select
      round(avg(pident), 2) as avg_pident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
      ,round(avg(begin_yr),1) as avg_begin_yr
      ,round(avg(end_yr),1) as avg_end_yr
      ,min(begin_yr), max(begin_yr)
      ,min(end_yr), max(end_yr)
from perc_hits_overall_ident_yr
where perc_hits > 95;

--Get files with perc hits < 95
select *
from perc_hits_overall_ident_yr
where perc_hits < 95;

--Get average perc_hits for those < 96
select
      round(avg(pident), 2) as avg_pident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
      ,round(avg(begin_yr),1) as avg_begin_yr
      ,round(avg(end_yr),1) as avg_end_yr
      ,min(begin_yr), max(begin_yr)
      ,min(end_yr), max(end_yr)
from perc_hits_overall_ident_yr
where perc_hits < 96;

--Trends by end_yr
------------------

--By year
select end_yr
      ,count(*)
      ,round(avg(pident), 2) as avg_pident
      ,round(avg(qident) * 100, 2) as avg_qident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_yr
group by end_yr
order by end_yr;

-- Year < 2005
select 
      count(*)
      ,round(avg(pident), 2) as avg_pident
      ,round(avg(qident) * 100, 2) as avg_qident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_yr
where end_yr < 2005
order by end_yr;

-- Year > 2004
select 
      count(*)
      ,round(avg(pident), 2) as avg_pident
      ,round(avg(qident) * 100, 2) as avg_qident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_yr
where end_yr > 2004
order by end_yr;

--====================================
--Analysis for seq year = NULL
--====================================

--Number of seq sets with no year info
select count(*)
from chimp_seq_year
where min_date is null and
      max_date is null;

--Get filenames
select seqfile_id
from chimp_seq_year
where min_date is null and
      max_date is null
order by seqfile_id; 

--Number of hits returned per file containing NULL year data
create view perc_hits_overall_ident_null_yr as
  select seqfile
  ,avg(pident) as pident
  ,avg(gapopen) as gapopen
  ,avg(mismatch) as mismatch
  ,avg(length) as aln_len
  ,avg(qlen) as qlen
  ,avg(nident)/avg(qlen) as qident
  ,(count(*)*100/25000) as perc_hits
  ,(sum(nident) / (25000 * avg(qlen)))*100  as overall_ident
  ,cast(min_date as int) as begin_yr
  ,cast(max_date as int) as end_yr
from chimp_blast_on_homo join chimp_seq_year
  on seqfile = seqfile_id
where min_date is null and
      max_date is null
group by seqfile
order by seqfile + 0 asc;

--Check perc_hits_yr view
select * from perc_hits_overall_ident_null_yr;

drop view perc_hits_overall_ident_null_yr;

--Get averages
select 
      round(avg(pident), 2) as avg_pident
      ,round(avg(qident), 2) as avg_qident
      ,round(avg(gapopen), 2) as avg_gapopen
      ,round(avg(mismatch), 2) as avg_mismatch
      ,round(avg(aln_len),1) as avg_aln_len
      ,round(avg(perc_hits), 2) as avg_perc_hits
      ,round(avg(overall_ident), 2) as overall_ident
from perc_hits_overall_ident_null_yr;


----------------------------------------------------------------------
--Taking a look at the non-hitters (on human) and their BLAST results
--on chimpanzee
----------------------------------------------------------------------

--Number of non-hitters on human
select
    (select count(*) from chimp_seq_data) -
    (select count(*) from chimp_blast_on_homo)
    as num_nonhitters;

--Number of nonhitters hit on chimp
select count(*) from nonhitter_blast_on_chimp;

--Averages for non-hitters on chimp
select count(*)
      ,round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
from nonhitter_blast_on_chimp;

--Averages for non-hitters on chimp by seqfile
select count(*)
      ,round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
from nonhitter_blast_on_chimp
group by seqfile
order by seqfile + 0 asc;

--Number of hits returned per file containing year data
create view nonhitters_yr as
  select seqfile
  ,avg(pident) as pident
  ,avg(gapopen) as gapopen
  ,avg(mismatch) as mismatch
  ,cast(min_date as int) as begin_yr
  ,cast(max_date as int) as end_yr
from nonhitter_blast_on_chimp join chimp_seq_year
  on seqfile = seqfile_id
where min_date not null and
      max_date not null
group by seqfile
order by seqfile + 0 asc;

drop view nonhitters_yr;

select * from nonhitters_yr;

--Get aln averages for nonhitters on chimp with yr data
select 
      round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
from nonhitters_yr;

--Nonhitters on chimp avg by end_yr
select 
      round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
      ,end_yr
from nonhitters_yr
group by end_yr
order by end_yr;

--Nonhitters on chimp avg where end_yr < 2005
select count(*)
      ,round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
from nonhitters_yr
where end_yr < 2005;

--Nonhitters on chimp avg where end_yr > 2004
select count(*)
      ,round(avg(pident), 2) as pident
      ,round(avg(gapopen), 2) as gapopen
      ,round(avg(mismatch), 2) as mismatch
from nonhitters_yr
where end_yr > 2004;
