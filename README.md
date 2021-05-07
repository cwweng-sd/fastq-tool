# fastq-tool

This perl script could be used to parse the FASTQ file for various respects of processing.

Date: 2016/10/14<br>
Upgrade: 2018/1/23<br>
Author: Jeff<br>
<br>
Usage: fastq_tool.pl [OPTIONS] < in.fastq<br>
<br>
Options:<br>
-h, -?, --help<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This help message.<br>
-cr, --count-read<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Count total reads and bases from FASTQ file.<br>
-c,  --clip-primer<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fetch reads that contained forward primer and reverse primer sequences on read 5'-end and read 3'-end, respectively, and then clip primers.<br>
-fl, --fetch-long-read [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Filter out the length of read shorter than the user-defined.<br>
-fs, --fetch-short-read [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Filter out the length of read longer than the user-defined.<br>
-fa, --fastq-to-fasta<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Convert FASTQ format to FASTA format.<br>
-fp, --fetch-paired-read [reverse_R2.fastq]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Compare the FASTQ R1/R2 files and then output the paired reads into individual forward and reverse files. The FASTQ file of reverse reads need be named the reverse_R2.fastq! The ID of each read must have the original illumina format!<br>
-m,  --merge-pair-read [mismatch] [reverse_R2.fastq]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Compare the paired-reads and then output forward reads that present similar context of DNA with reverse reads. The length of forward and reverse reads need be same. The FASTQ file of reverse reads need be named the reverse_R2.fastq! The ID of each read must have the original illumina format!<br>
-rc, --reverse-complement<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Reverse complement of reads.<br>
-rl, --re-label-read [read_label]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Re-label reads to user defined label as prefix of label. The suffix of label will use continuous numeric format, such as 1, 2, 3 and so on.<br>
-tf, --trim-read-five-prime [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Trim reads from 5'-end according to the user-defined length.<br>
-tt, --trim-read-three-prime [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Trim reads from 3'-end according to the user-defined length.<br>
-ff, --fetch-read-five-prime [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fetch subsequence from 5'-end of reads according to the user-defined length.<br>
-ft, --fetch-read-three-prime [length]<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Fetch subsequence from 3'-end of reads according to the user-defined length.<br>
<br>
If the flag of -c is used, please provide the primer file including forward and reserve sequences.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Please refer the example file, i.e. primer.txt.<br>
