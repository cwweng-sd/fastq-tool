# This script used to parse the FASTQ file for various respects of processing.
#
# Date:2016/10/14
# Upgrade: 2018/1/23
# Author:C.W.Weng (Jeff)
 
#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

my $opts=parse_params();
do_pileup_to_fastq($opts);

exit;

#----------------

sub error{
    my(@msg)=@_;
    if(scalar@msg){croak(@msg);}
    die
        "Usage: fastq_tool.pl [OPTIONS] < in.fastq\n",
        "\n",
	"Options:\n",
        "   -h,  -?, --help                                          This help message.\n",
        "\n",
        "   -cr, --count-read                                        Count total reads and bases from FASTQ file.\n",
        "\n",
        "   -c,  --clip-primer                                       Fetch reads that contained forward primer and reverse primer sequences\n",
        "                                                            on read 5'-end and read 3'-end, respectively, and then clip primers.\n",
        "\n",
        "   -fl, --fetch-long-read [length]                          Filter out the length of read shorter than the user-defined.\n",
        "\n",
        "   -fs, --fetch-short-read [length]                         Filter out the length of read longer than the user-defined.\n",
        "\n",
        "   -fa, --fastq-to-fasta                                    Convert FASTQ format to FASTA format.\n",
        "\n",
        "   -fp, --fetch-paired-read [reverse_R2.fastq]              Compare the FASTQ R1/R2 files and then output the paired reads\n",
        "                                                            into individual forward and reverse files.\n",
        "                                                            The FASTQ file of reverse reads need be named the reverse_R2.fastq!\n",
        "                                                            The ID of each read must have the original illumina format!\n",
        "\n",
        "   -m,  --merge-pair-read [mismatch] [reverse_R2.fastq]     Compare the paired-reads and then output forward reads that\n",
        "                                                            present similar context of DNA with reverse reads.\n",
        "                                                            The length of forward and reverse reads need be same.\n",
        "                                                            The FASTQ file of reverse reads need be named the reverse_R2.fastq!\n",
        "                                                            The ID of each read must have the original illumina format!\n",
        "\n",
        "   -rc, --reverse-complement                                Reverse complement of reads.\n",
        "\n",
        "   -rl, --re-label-read [read_label]                        Re-label reads to user defined label as prefix of label.\n",
        "                                                            The suffix of label will use continuous numeric format, such as 1, 2, 3 and so on.\n",
        "\n",
        "   -tf, --trim-read-five-prime [length]                     Trim reads from 5'-end according to the user-defined length.\n",
        "\n",
        "   -tt, --trim-read-three-prime [length]                    Trim reads from 3'-end according to the user-defined length.\n",
        "\n",
        "   -ff, --fetch-read-five-prime [length]                    Fetch subsequence from 5'-end of reads according to the user-defined length.\n",
	"\n",
	"   -ft, --fetch-read-three-prime [length]                   Fetch subsequence from 3'-end of reads according to the user-defined length.\n",
	"\n",
	"If the flag of -c is used, please provide the primer file including forward and reserve sequences.\n",
        "Please refer the example file, i.e. primer.txt.\n",
        "\n";
}

sub parse_params{
    my %opts=();
    $opts{count_read}=0;
    $opts{clip_primer}=0;
    $opts{fetch_long_read}=0;
    $opts{fetch_short_read}=0;
    $opts{merge_pair_read}=0;
    $opts{fastq_to_fasta}=0;
    $opts{fetch_paired_read}=0;
    $opts{reverse_complement}=0;
    $opts{trim_read_five_prime}=0;
    $opts{trim_read_three_prime}=0;
    $opts{fetch_read_five_prime}=0;
    $opts{fetch_read_three_prime}=0;
    $opts{re_label_read}=0;

    $opts{fh_in}=*STDIN;

    while(my $arg=shift(@ARGV)){
        if($arg=~m/\W+\w*/){
            if($arg eq '-cr' || $arg eq '--count-read') { $opts{count_read}=1; next; }
            if($arg eq '-c'  || $arg eq '--clip-primer') { $opts{clip_primer}=1; next; }
            if($arg eq '-fl' || $arg eq '--fetch-long-read') { $opts{fetch_long_read}=1; next; }
            if($arg eq '-fs' || $arg eq '--fetch-short-read') { $opts{fetch_short_read}=1; next; }
            if($arg eq '-m'  || $arg eq '--merge-pair-read') { $opts{merge_pair_read}=1; next; }
            if($arg eq '-fa' || $arg eq '--fastq-to-fasta') { $opts{fastq_to_fasta}=1; next; }
            if($arg eq '-fp' || $arg eq '--fetch-paired-read') { $opts{fetch_paired_read}=1; next; }
            if($arg eq '-rc' || $arg eq '--reverse-complement') { $opts{reverse_complement}=1; next; }
            if($arg eq '-tf' || $arg eq '--trim-read-five-prime') { $opts{trim_read_five_prime}=1; next; }
            if($arg eq '-tt' || $arg eq '--trim-read-three-prime') { $opts{trim_read_three_prime}=1; next; }
            if($arg eq '-ff' || $arg eq '--fetch-read-five-prime') { $opts{fetch_read_five_prime}=1; next; }
	    if($arg eq '-ft' || $arg eq '--fetch-read-three-prime') { $opts{fetch_read_three_prime}=1; next; }
	    if($arg eq '-rl' || $arg eq '--re-label-read') { $opts{re_label_read}=1; next; }
            if($arg eq '-?'  || $arg eq '-h' || $arg eq '--help') { error(); }
        }
        if($arg eq "reverse_R2.fastq"){ $opts{reverse_reads_file}="reverse_R2.fastq"; next; }
        if($arg=~m/\d+/){
            if($opts{fetch_long_read}!=0 || $opts{fetch_short_read}!=0){
                if($arg<=0){ print "The length of read need larger than ZERO!\n"; }
                if($arg>0){ $opts{user_define_len}=$arg; }
                next;
            }elsif($opts{merge_pair_read}!=0){
                if($arg<0){ print "The number of mismatch need larger than or equal to ZERO!\n"; }
                if($arg>=0){ $opts{mismatch_base}=$arg; }
                next;
            }elsif($opts{trim_read_five_prime}!=0 || $opts{trim_read_three_prime}!=0){
                if($arg<0){ print "The defined length for trimming need larger than ZERO!\n"; }
                if($arg>=0){ $opts{user_define_len_trim}=$arg; }
                next;
            }elsif($opts{fetch_read_five_prime}!=0 || $opts{fetch_read_three_prime}!=0){
		if($arg<0){ print "The defined length of sub-sequence need larger than ZERO!\n"; }
		if($arg>=0){ $opts{user_define_len_capture}=$arg; }
		next;
	    }
        }
        if($arg=~m/\w+/ && $arg ne "reverse_R2.fastq"){ $opts{read_label_prefix}=$arg; next; }

        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }

    return \%opts;
}

sub do_pileup_to_fastq{
    my ($opts)=@_;
 
    my $count_read=$$opts{count_read} ? 1 : 0;
    my $clip_primer=$$opts{clip_primer} ? 1 : 0;
    my $fetch_long_read=$$opts{fetch_long_read} ? 1 : 0;
    my $fetch_short_read=$$opts{fetch_short_read} ? 1 : 0;
    my $merge_pair_read=$$opts{merge_pair_read} ? 1 : 0;
    my $fastq_to_fasta=$$opts{fastq_to_fasta} ? 1 : 0;
    my $fetch_paired_read=$$opts{fetch_paired_read} ? 1 : 0;
    my $reverse_complement=$$opts{reverse_complement} ? 1 : 0;
    my $trim_read_five_prime=$$opts{trim_read_five_prime} ? 1 : 0;
    my $trim_read_three_prime=$$opts{trim_read_three_prime} ? 1 : 0;
    my $fetch_read_five_prime=$$opts{fetch_read_five_prime} ? 1 : 0;
    my $fetch_read_three_prime=$$opts{fetch_read_three_prime} ? 1 : 0;
    my $re_label_read=$$opts{re_label_read} ? 1 : 0;

    if($count_read){ count_read(); }
    if($clip_primer){ clip_primer(); }
    if($fetch_long_read){ fetch_long_read(); }
    if($fetch_short_read){ fetch_short_read(); }
    if($merge_pair_read){ merge_pair_read(); }
    if($fastq_to_fasta){ fastq_to_fasta(); }
    if($fetch_paired_read){ fetch_paired_read(); }
    if($reverse_complement){ reverse_complement(); }
    if($trim_read_five_prime){ trim_read_five_prime(); }
    if($trim_read_three_prime){ trim_read_three_prime(); }
    if($fetch_read_five_prime){ fetch_read_five_prime(); }
    if($fetch_read_three_prime){ fetch_read_three_prime(); }
    if($re_label_read){ re_label_read(); }
}

sub count_read{
    open(OUTPUT,">count_read.txt")||die "Unable to output the summary of read count!\n";

    my $fh_in=$$opts{fh_in};
    my $read_base=0;
    my $total_bases=0;
    my $total_reads=0;
    while(<$fh_in>){
        $total_reads++;
        $_=<$fh_in>;
        chomp($_);
        my $read=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        $read_base=length($read);
        $total_bases=$total_bases+$read_base;
    }
    print OUTPUT "Total reads:$total_reads\tTotal bases:$total_bases\n";

    close OUTPUT; 
}

sub clip_primer{
    open(OUTPUT,">extracted_trim.fastq")||die "Unable to output the trimmed reads!\n";

    # process the primer file
    my @primers=parse_primers();
    my $forward_primer=$primers[0];
    my $revcomp_forward_primer=$primers[1];
    my $reverse_primer=$primers[2];
    my $revcomp_reverse_primer=$primers[3];
    my @forward_primer_ary=();
    my @revcomp_forward_primer_ary=();
    my @reverse_primer_ary=();
    my @revcomp_reverse_primer_ary=();
    if(length($forward_primer)==length($reverse_primer)){
        @forward_primer_ary=trim_primer($forward_primer);
        @revcomp_forward_primer_ary=trim_primer_revcomp($revcomp_forward_primer);
        @reverse_primer_ary=trim_primer($reverse_primer);
        @revcomp_reverse_primer_ary=trim_primer_revcomp($revcomp_reverse_primer);
        if($#forward_primer_ary==$#revcomp_forward_primer_ary && $#forward_primer_ary==$#reverse_primer_ary &&
           $#forward_primer_ary==$#revcomp_reverse_primer_ary){
            for(my $i=0;$i<=$#forward_primer_ary;$i++){
                print length($forward_primer_ary[$i])." mer\n";
                print "Forward primer: $forward_primer_ary[$i]\n";
                print "Revcomp of forward primer: $revcomp_forward_primer_ary[$i]\n";
                print "Reverse primer: $reverse_primer_ary[$i]\n";
                print "Revcomp of reverse primer: $revcomp_reverse_primer_ary[$i]\n";
            }
        }
    }
    print "\n";
    print "-----------------------------------------";
    print "\n";

    # fetch primer-contained reads
    my $fh_in=$$opts{fh_in};
    my $read_id;
    my $seq;
    my $seq_len;
    my $seq_quality;
    my $trim_seq;
    my $trim_seq_quality;
    my $trim_seq_len;
    my $trim_pair_seq;
    my $trim_pair_seq_quality;
    my $prime5seq_front_match_fprimer;
    my $prime5seq_end_match_fprimer;
    my $prime3seq_front_match_rprimer_revcomp;
    my $prime3seq_end_match_rprimer_revcomp;
    my $prime5seq_front_match_rprimer;
    my $prime5seq_end_match_rprimer;
    my $prime3seq_front_match_fprimer_revcomp;
    my $prime3seq_end_match_fprimer_revcomp;

    while(<$fh_in>){
        $read_id=$_;
        $_=<$fh_in>;
        chomp($_);
        $seq=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        chomp($_);
        $seq_quality=$_;
        $seq_len=length($seq);
        if($#forward_primer_ary==$#revcomp_forward_primer_ary && $#forward_primer_ary==$#reverse_primer_ary &&
           $#forward_primer_ary==$#revcomp_reverse_primer_ary){
            my @prime_5_seq_front=();
            my @prime_5_seq_end=();
            my @prime_3_seq_front=();
            my @prime_3_seq_end=();
            for(my $i=0;$i<=$#forward_primer_ary;$i++){
                # for @prime_5_seq_front, @prime_3_seq_end
                my $sub_seq_len_posi_front_start=$i;
                my $sub_seq_len_posi_front_end=length($forward_primer_ary[$i]);
                my $sub_seq_len_nega_end_start=-length($forward_primer_ary[0]);
                my $sub_seq_len_nega_end_end=length($forward_primer_ary[$i]);
                my $prime_5_seq_front=substr $seq,$sub_seq_len_posi_front_start,$sub_seq_len_posi_front_end;
                my $prime_3_seq_end=substr $seq,$sub_seq_len_nega_end_start,$sub_seq_len_nega_end_end;
                push @prime_5_seq_front,$prime_5_seq_front;
                push @prime_3_seq_end,$prime_3_seq_end;
                # for @prime_5_seq_end, @prime_3_seq_front
                my $sub_seq_len_posi_end_start=0;
                my $sub_seq_len_posi_end_end=length($forward_primer_ary[0])-$i;
                my $sub_seq_len_nega_front_start=-length($forward_primer_ary[$i]);
                my $prime_5_seq_end=substr $seq,$sub_seq_len_posi_end_start,$sub_seq_len_posi_end_end;
                my $prime_3_seq_front=substr $seq,$sub_seq_len_nega_front_start;
                push @prime_5_seq_end,$prime_5_seq_end;
                push @prime_3_seq_front,$prime_3_seq_front;
            }
            print "$read_id";
            if($#forward_primer_ary==$#prime_5_seq_front){
                $prime5seq_front_match_fprimer="false";
                print "prime5seq-front to forward: ";
                for(my $i=0;$i<=$#forward_primer_ary;$i++){
                    if($prime_5_seq_front[$i] eq $forward_primer_ary[$i]){
                        $prime5seq_front_match_fprimer="true";
                        print length($prime_5_seq_front[$i])."\t";
                        if(length($seq)==$seq_len){
                            $trim_seq=substr $seq,length($prime_5_seq_front[$i]);
                            $trim_seq_quality=substr $seq_quality,length($prime_5_seq_front[$i]);
                            $trim_seq_len=length($trim_seq);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#forward_primer_ary==$#prime_5_seq_end){
                $prime5seq_end_match_fprimer="false";
                print "prime5seq-end to forward: ";
                for(my $i=0;$i<=$#forward_primer_ary;$i++){
                    if($prime_5_seq_end[$i] eq $forward_primer_ary[$i]){
                        $prime5seq_end_match_fprimer="true";
                        print length($prime_5_seq_end[$i])."\t";
                        if(length($seq)==$seq_len){
                            $trim_seq=substr $seq,length($prime_5_seq_end[$i]);
                            $trim_seq_quality=substr $seq_quality,length($prime_5_seq_end[$i]);
                            $trim_seq_len=length($trim_seq);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#reverse_primer_ary==$#prime_5_seq_front){
                $prime5seq_front_match_rprimer="false";
                print "prime5seq-front to reverse: ";
                for(my $i=0;$i<=$#reverse_primer_ary;$i++){
                    if($prime_5_seq_front[$i] eq $reverse_primer_ary[$i]){
                        $prime5seq_front_match_rprimer="true";
                        print length($prime_5_seq_front[$i])."\t";
                        if(length($seq)==$seq_len){
                            $trim_seq=substr $seq,length($prime_5_seq_front[$i]);
                            $trim_seq_quality=substr $seq_quality,length($prime_5_seq_front[$i]);
                            $trim_seq_len=length($trim_seq);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#reverse_primer_ary==$#prime_5_seq_end){
                $prime5seq_end_match_rprimer="false";
                print "prime5seq-end to reverse: ";
                for(my $i=0;$i<=$#reverse_primer_ary;$i++){
                    if($prime_5_seq_end[$i] eq $reverse_primer_ary[$i]){
                        $prime5seq_end_match_rprimer="true";
                        print length($prime_5_seq_end[$i])."\t";
                        if(length($seq)==$seq_len){
                            $trim_seq=substr $seq,length($prime_5_seq_end[$i]);
                            $trim_seq_quality=substr $seq_quality,length($prime_5_seq_end[$i]);
                            $trim_seq_len=length($trim_seq);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#revcomp_forward_primer_ary==$#prime_3_seq_front){
                $prime3seq_front_match_fprimer_revcomp="false";
                print "prime3seq-front to revcomp forward: ";
                for(my $i=0;$i<=$#revcomp_forward_primer_ary;$i++){
                    if($prime_3_seq_front[$i] eq $revcomp_forward_primer_ary[$i]){
                        $prime3seq_front_match_fprimer_revcomp="true";
                        print length($prime_3_seq_front[$i])."\t";
                        if(length($trim_seq)==$trim_seq_len){
                            $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                            $trim_pair_seq_quality=substr $trim_seq_quality,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#revcomp_forward_primer_ary==$#prime_3_seq_end){
                $prime3seq_end_match_fprimer_revcomp="false";
                print "prime3seq-end to revcomp forward: ";
                for(my $i=0;$i<=$#revcomp_forward_primer_ary;$i++){
                    if($prime_3_seq_end[$i] eq $revcomp_forward_primer_ary[$i]){
                        $prime3seq_end_match_fprimer_revcomp="true";
                        print length($prime_3_seq_end[$i])."\t";
                        if(length($trim_seq)==$trim_seq_len){
                            $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                            $trim_pair_seq_quality=substr $trim_seq_quality,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#revcomp_reverse_primer_ary==$#prime_3_seq_front){
                $prime3seq_front_match_rprimer_revcomp="false";
                print "prime3seq-front to revcomp reverse: ";
                for(my $i=0;$i<=$#revcomp_reverse_primer_ary;$i++){
                    if($prime_3_seq_front[$i] eq $revcomp_reverse_primer_ary[$i]){
                        $prime3seq_front_match_rprimer_revcomp="true";
                        print length($prime_3_seq_front[$i])."\t";
                        if(length($trim_seq)==$trim_seq_len){
                            $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                            $trim_pair_seq_quality=substr $trim_seq_quality,0,length($trim_seq)-length($prime_3_seq_front[$i]);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if($#reverse_primer_ary==$#prime_3_seq_end){
                $prime3seq_end_match_rprimer_revcomp="false";
                print "prime3seq-end to revcomp reverse: ";
                for(my $i=0;$i<=$#revcomp_reverse_primer_ary;$i++){
                    if($prime_3_seq_end[$i] eq $revcomp_reverse_primer_ary[$i]){
                        $prime3seq_end_match_rprimer_revcomp="true";
                        print length($prime_3_seq_end[$i])."\t";
                        if(length($trim_seq)==$trim_seq_len){
                            $trim_pair_seq=substr $trim_seq,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                            $trim_pair_seq_quality=substr $trim_seq_quality,0,length($trim_seq)-length($prime_3_seq_end[$i]);
                        }
                        last;
                    }
                }
                print "\n";
            }
            if(($prime5seq_front_match_fprimer eq "true" && $prime3seq_end_match_rprimer_revcomp eq "true") ||
               ($prime5seq_front_match_rprimer eq "true" && $prime3seq_end_match_fprimer_revcomp eq "true") ||
               ($prime5seq_front_match_fprimer eq "true" && $prime3seq_front_match_rprimer_revcomp eq "true") ||
               ($prime5seq_front_match_rprimer eq "true" && $prime3seq_front_match_fprimer_revcomp eq "true") ||
               ($prime5seq_end_match_fprimer eq "true" && $prime3seq_front_match_rprimer_revcomp eq "true") ||
               ($prime5seq_end_match_rprimer eq "true" && $prime3seq_front_match_fprimer_revcomp eq "true") ||
               ($prime5seq_end_match_fprimer eq "true" && $prime3seq_end_match_rprimer_revcomp eq "true") ||
               ($prime5seq_end_match_rprimer eq "true" && $prime3seq_end_match_fprimer_revcomp eq "true")){
                print OUTPUT $read_id;
                print OUTPUT "$trim_pair_seq\n";
                print OUTPUT "+\n";
                print OUTPUT "$trim_pair_seq_quality\n";
                print "output\n";
            }
            print "\n";
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_long_read{
    open(OUTPUT,">defined_long_reads.fastq")||die "Unable to output the user-defined long reads!\n";

    my $fh_in=$$opts{fh_in};
    my $define_len=$$opts{user_define_len};
    while(<$fh_in>){
        my $read_id=$_;
        $_=<$fh_in>;
        chomp($_);
        my $seq=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        if(length($seq)>=$define_len){
            print OUTPUT "$read_id";
            print OUTPUT "$seq\n";
            print OUTPUT "+\n";
            print OUTPUT "$_";
        }else{
            next;
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_short_read{
    open(OUTPUT,">defined_short_reads.fastq")||die "Unable to output the user-defined short reads!\n";

    my $fh_in=$$opts{fh_in};
    my $define_len=$$opts{user_define_len};
    while(<$fh_in>){
        my $read_id=$_;
        $_=<$fh_in>;
        chomp($_);
        my $seq=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        if(length($seq)<=$define_len){
            print OUTPUT "$read_id";
            print OUTPUT "$seq\n";
            print OUTPUT "+\n";
            print OUTPUT "$_";
        }else{
            next;
        }
    }

    close $fh_in;
    close OUTPUT;
}

sub merge_pair_read{
    open(OUTPUT,">merge_pair_reads.fastq")||die "Unable to output the merged paired reads!\n";

    my $fh_in=$$opts{fh_in};
    my $reverse_file=$$opts{reverse_reads_file};
    open(FILE,$reverse_file)||die "The FASTQ file of reverse reads need be named the reverse_R2.fastq!\n";
    my @reverse_reads=<FILE>;
    close FILE;

    my $mismatch=$$opts{mismatch_base};

    while(<$fh_in>){
        my $read_id=$_;
        my @read_id_split1=();
        my @read_id_split2=();
        @read_id_split1=split ' ',$read_id;
        @read_id_split2=split ':',$read_id_split1[1];
        $_=<$fh_in>;
        chomp($_);
        my $seq=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        my $seq_quality=$_;

        my $reverse_read_id;
        my @reverse_read_id_split1;
        my @reverse_read_id_split2;
        my $reverse_seq;
        for(my $i=0;$i<=$#reverse_reads;$i++){
            my $index=$i+1;
            if($index%4==1){
                $reverse_read_id=$reverse_reads[$i];
                next;
            }
            if($index%4==2){
                chomp($reverse_reads[$i]);
                $reverse_seq=revcomp($reverse_reads[$i]);
                next;
            }
            if($index%4==3){ next; }
            if($index%4==0){
                @reverse_read_id_split1=();
                @reverse_read_id_split2=();
                @reverse_read_id_split1=split ' ',$reverse_read_id;
                @reverse_read_id_split2=split ':',$reverse_read_id_split1[1];
                if($read_id_split1[0] eq $reverse_read_id_split1[0]){
                    if($read_id_split2[0]==1 && $reverse_read_id_split2[0]==2){
                        if(length($seq)==length($reverse_seq)){
                            my $count_mismatch=0;
                            my $num_base=length($seq);
                            for(my $subseq_index=0;$subseq_index<$num_base;$subseq_index++){
                                my $base=substr $seq,$subseq_index,1;
                                my $reverse_base=substr $reverse_seq,$subseq_index,1;
                                if($base ne $reverse_base){
                                    $count_mismatch=$count_mismatch+1;
                                }
                            }
                            if($count_mismatch<=$mismatch){
                                print OUTPUT $read_id;
                                print OUTPUT "$seq\n";
                                print OUTPUT "+\n";
                                print OUTPUT "$seq_quality";
                                last;
                            }
                        }
                    }
                }
            }
        }
    }

    close $fh_in;
    close OUTPUT; 
}

sub fetch_paired_read{
    open(OUTPUT1,">forward_reads.fastq")||die "Unable to output the forward reads!\n";
    open(OUTPUT2,">reverse_reads.fastq")||die "Unable to output the reverse reads!\n";

    my $fh_in=$$opts{fh_in};
    my %forward_reads;
    while(<$fh_in>){
        my $read_id=$_;
        my @read_id_split=();
        @read_id_split=split ' ',$read_id;
        $_=<$fh_in>;
        my $seq=$_;
        $_=<$fh_in>;
        my $fastq_sign=$_;
        $_=<$fh_in>;
        my $seq_quality=$_;
        my $read_info=$read_id."\t".$seq."\t".$fastq_sign."\t"."$seq_quality";
        $forward_reads{$read_id_split[0]}=$read_info;
    }
    close $fh_in;

    my $reverse_file=$$opts{reverse_reads_file};
    open(FILE,$reverse_file)||die "The FASTQ file of reverse reads need be named the reverse_R2.fastq!\n";
    my %reverse_reads;
    while(<FILE>){
        my $read_id=$_;
        my @read_id_split=();
        @read_id_split=split ' ',$read_id;
        $_=<FILE>;
        my $seq=$_;
        $_=<FILE>;
        my $fastq_sign=$_;
        $_=<FILE>;
        my $seq_quality=$_;
        my $read_info=$read_id."\t".$seq."\t".$fastq_sign."\t"."$seq_quality";
        $reverse_reads{$read_id_split[0]}=$read_info;
    }
    close FILE;

    foreach my $key(keys %forward_reads){
        if(defined $forward_reads{$key} && defined $reverse_reads{$key}){
            my @forward_read=split '\t',$forward_reads{$key};
            foreach my $forward_data(@forward_read){
                print OUTPUT1 $forward_data;
            }
            my @reverse_read=split '\t',$reverse_reads{$key};
            foreach my $reverse_data(@reverse_read){
                print OUTPUT2 $reverse_data;
            }
        }
    }

    close OUTPUT1;
    close OUTPUT2;
}

sub fastq_to_fasta{
    open(OUTPUT,">seq.fa")||die "Unable to output the file of FASTA format!\n";

    my $fh_in=$$opts{fh_in};
    while(<$fh_in>){
        my $seq_desc=$_;
        $_=<$fh_in>;
        my $seq=$_;
        $_=<$fh_in>;
        $_=<$fh_in>;
        print OUTPUT ">$seq_desc";
        print OUTPUT "$seq";
    }

    close $fh_in;
    close OUTPUT;
}

sub reverse_complement{
    open(OUTPUT,">revcomp.fastq")||die "Unable to output the reverse-complement reads!\n";

    my $fh_in=$$opts{fh_in};
    while(<$fh_in>){
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
        my $revcomp_read=revcomp($_);
        print OUTPUT "$revcomp_read\n";
        $_=<$fh_in>;
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
        my $rev_read_qual=reverse($_);
        print OUTPUT "$rev_read_qual\n";
    }

    close $fh_in;
    close OUTPUT;
}

sub trim_read_five_prime{
    open(OUTPUT,">reads_trim_five_prime.fastq")||die "Unable to output the 5'-end trimmed reads!\n";

    my $fh_in=$$opts{fh_in};
    my $define_len_trim=$$opts{user_define_len_trim};
    while(<$fh_in>){
        print OUTPUT $_;
        $_=<$fh_in>;
        my $trim_read=substr $_,$define_len_trim;
        print OUTPUT $trim_read;
        $_=<$fh_in>;
        print OUTPUT $_;
        $_=<$fh_in>;
        my $trim_read_qual=substr $_,$define_len_trim;
        print OUTPUT $trim_read_qual;
    }

    close $fh_in;
    close OUTPUT;
}

sub trim_read_three_prime{
    open(OUTPUT,">reads_trim_three_prime.fastq")||die "Unable to output the 3'-end trimmed reads!\n";

    my $fh_in=$$opts{fh_in};
    my $define_len_trim=$$opts{user_define_len_trim};
    while(<$fh_in>){ 
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
        my $trim_read=substr $_,0,-$define_len_trim;
        print OUTPUT "$trim_read\n";
        $_=<$fh_in>; 
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
        my $trim_read_qual=substr $_,0,-$define_len_trim;
        print OUTPUT "$trim_read_qual\n";
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_read_five_prime{
    open(OUTPUT,">subseq_five_prime.fastq")||die "Unable to output sub-sequence from the directuon of 5'-end!\n";
    my $fh_in=$$opts{fh_in};
    my $define_len=$$opts{user_define_len_capture};
    while(<$fh_in>){
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
	my $subseq=substr $_,0,$define_len;
        print OUTPUT "$subseq\n";
        $_=<$fh_in>;
        print OUTPUT $_;
        $_=<$fh_in>;
        chomp($_);
	my $read_qual=substr $_,0,$define_len;
        print OUTPUT "$read_qual\n";
    }

    close $fh_in;
    close OUTPUT;
}

sub fetch_read_three_prime{
    open(OUTPUT,">subseq_three_prime.fastq")||die "Unable to output sub-sequence from the directuon of 3'-end!\n";
    my $fh_in=$$opts{fh_in};
    my $define_len=$$opts{user_define_len_capture};
    while(<$fh_in>){
        print OUTPUT $_;
	$_=<$fh_in>;
        chomp($_);
	my $subseq=substr $_,-$define_len;
	print OUTPUT "$subseq\n";
	$_=<$fh_in>;
	print OUTPUT $_;
	$_=<$fh_in>;
	chomp($_);
	my $read_qual=substr $_,-$define_len;
	print OUTPUT "$read_qual\n"; 
    }

    close $fh_in;
    close OUTPUT;
}

sub re_label_read{
    open(OUTPUT,">re_label_read.fastq")||die "Unable to output the re-labelled reads!\n";

    my $fh_in=$$opts{fh_in};
    my $read_prefix=$$opts{read_label_prefix};
    my $count_read=0; 
    while(<$fh_in>){
        $count_read+=1;
        my $read_label=$read_prefix.$count_read;
        print OUTPUT "$read_label\n";
        $_=<$fh_in>;
        print OUTPUT "$_";
        $_=<$fh_in>;
        print OUTPUT "$_";
        $_=<$fh_in>;
        print OUTPUT "$_";
    }

    close OUTPUT;
}

sub parse_primers{
    my @total_primers=();
    my $forward_primer;
    my $revcomp_forward_primer;
    my $reverse_primer;
    my $revcomp_reverse_primer;
    open(PRIMER_FILE,"primer.txt")||die "Please provide the primer file in the current directory.\nThe example primer file can be a reference.\n";
    my @primers=<PRIMER_FILE>;
    close PRIMER_FILE;
    foreach my $primer(@primers){
        chomp($primer);
        my $primer_direction=substr $primer,0,1;
        if($primer_direction eq 'F'){
            $primer=~s/F://g;
            $forward_primer=$primer;
            $revcomp_forward_primer=revcomp($forward_primer);
            next;
        }
        if($primer_direction eq 'R'){
            $primer=~s/R://g;
            $reverse_primer=$primer;
            $revcomp_reverse_primer=revcomp($reverse_primer);
            next;
        }
    }
    push @total_primers,$forward_primer;
    push @total_primers,$revcomp_forward_primer;
    push @total_primers,$reverse_primer;
    push @total_primers,$revcomp_reverse_primer;
    # retern the primers: forward primer, revcomp of forward primer, reverse primer, revcomp of reverse primer
    return @total_primers;
}

sub trim_primer{
    my $primer=$_[0];
    my $mer_len=abs(3-length($primer));
    my @primer_ary=();
    for(my $i=0;$i<=$mer_len;$i++){
        my $primer_seq=substr $primer,$i;
        push @primer_ary,$primer_seq;
    }
    return @primer_ary;
}

sub trim_primer_revcomp{
    my $primer=$_[0];
    my $mer_len=abs(3-length($primer));
    my @revcomp_primer_ary=();
    for(my $i=0;$i<=$mer_len;$i++){
        my $revcomp_primer_seq=substr $primer,0,length($primer)-$i;
        push @revcomp_primer_ary,$revcomp_primer_seq;
    }
    return @revcomp_primer_ary;
}

sub revcomp{
   my $seq=$_[0];
   $seq=~tr/ACGTacgt/TGCAtgca/;
   $seq=reverse($seq);
   return $seq;
}
