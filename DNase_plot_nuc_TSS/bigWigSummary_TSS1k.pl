#!/usr/bin/perl -w
use strict;

open IN,"/WPS1/lilin/Database/mm9_annotation.20150821/refGene.TSS_1000.bed" or die $!;
my $i;

while(<IN>){
	chomp;
	my @a = split /\s+/,$_;
	if ( $a[3] eq "+" ){
		my $tmp1= `/WPS1/huboqiang/software/UCSC/bigWigSummary $ARGV[0] $a[0] $a[1] $a[2] 100`;
		if ( $tmp1 ){
			print "$a[0]:$a[1]-$a[2]\t$tmp1";
		}
	}
	elsif ( $a[3] eq "-" ){
		my $tmp =`/WPS1/huboqiang/software/UCSC/bigWigSummary $ARGV[0] $a[0] $a[1] $a[2] 100`;
		if ( $tmp ){
			my @out=split /\s+/,$tmp;
			print "$a[0]:$a[1]-$a[2]";
			foreach $i(1..100){
				my $word=pop @out;
				print "\t$word";
			}print "\n";
		}
	}
}
close IN;
 
