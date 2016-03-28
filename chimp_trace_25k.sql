--Various sql queries on the chimp_trace_25k db
--Written by JP Tomkins, March, 2016
--Email: jtomkins@icr.org
--Institute for Creation Research

--Head check tables
select * from blast_data limit 5;

select * from seq_data limit 5;

select * from seq_year limit 5;

--Table row count checks
select count(*)
from blast_data;

select count(*)
from seq_data;

select count(*)
from seq_year;

--Number of seq sets with year info
select count(*)
from seq_year
where min_date not null;

--Get the percent hits for each seq file
select round(count(*)*100/25000, 2) as perc_hits, seqfile
from blast_data
group by seqfile;

--Get percent rows less than 90% ident
SELECT round(count(*)*100/2525000, 2)
from blast_data
where pident < 90.00;

--Get percent rows less than 80% ident
SELECT count(*) as num_less_than_80
from blast_data
where pident < 80.00;

--Get avg data in specific ranges
SELECT count(*), round(avg(pident),1), round(avg(gapopen),1),
  round(avg(mismatch),1)
from blast_data
where pident between 90.00 and 100.00;

SELECT count(*), round(avg(pident),1), round(avg(gapopen),1)
from blast_data
where pident between 80.00 and 90.00;

SELECT count(*), round(avg(pident),1), round(avg(gapopen),1)
from blast_data
where pident between 70.00 and 80.00;

SELECT count(*), round(avg(pident),1), round(avg(gapopen),1)
from blast_data
where pident between 60.00 and 70.00;

--Get averages for each seqfile
select seqfile, round(avg(pident),1), round(avg(gapopen),1), 
  round(avg(mismatch),1), round(avg(qlen),1), round(avg(nident),1)
from blast_data
group by seqfile
order by avg(pident) desc;

--Create a seqfile_data view and use for getting grouped seqfile info
create view seqfile_data as
select seqfile 
      ,round(avg(pident),1) as pident 
      ,round(avg(gapopen),1) as gapopen 
      ,round(avg(mismatch),1) as mismatch 
      ,round(avg(length), 1) as alnlen
      ,round(avg(qlen),1) as qlen 
      ,round(avg(nident),1) as nident
from blast_data
group by seqfile
order by avg(pident) desc;

--Check the seq_data table
select *
from seqfile_data;

--Create seqdat_yr view joining the seqfile_data view and seq_year table
create view seqdat_yr as
select seqfile
      ,round((nident/qlen)*100, 1) as qlenident
      ,pident
      ,gapopen
      ,mismatch
      ,cast(min_date as int) as begin_yr 
      ,cast(max_date as int) as end_yr
from seqfile_data join seq_year
  on seqfile = seqfile_id
where min_date not null
order by pident desc;

--Check the seq_yr table
select * 
from seqdat_yr;

--Trends by end_yr
select round(avg(qlenident),1), round(avg(pident),1), 
  round(avg(gapopen),1), round(avg(mismatch),1), end_yr,
  count(seqfile)
from seqdat_yr
group by end_yr
order by end_yr;

--Trends by year, sorted by begin_yr
select round(avg(qlenident),1), round(avg(pident),1), 
  round(avg(gapopen),1), round(avg(mismatch),1), begin_yr,
  end_yr, count(seqfile)
from seqdat_yr
group by begin_yr
order by begin_yr;

--Trends by begin_yr between 2000 and 2003
select round(avg(qlenident),1), round(avg(pident),1), 
  round(avg(gapopen),1), round(avg(mismatch),1), count(seqfile)
from seqdat_yr
where begin_yr between 2000 and 2003;

--Trends by begin_yr between 2004 and 2007
select round(avg(qlenident),1), round(avg(pident),1), 
  round(avg(gapopen),1), round(avg(mismatch),1), count(seqfile)
from seqdat_yr
where begin_yr between 2004 and 2007;