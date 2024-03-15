#! /usr/bin/perl -w
###########################
# Author: Shawn Yost <seyost@ucsd.edu>
# 
############################

use warnings;
use strict;
use Getopt::Long;

sub convert($$);
sub get_allele_freq($$$);

my @gold_het_qscore;

my $h;
my $gold = "tmp";
my $sample = "tmp";
my $quartile = 0.1;
my $mmRate_column = 12;
my $ref_column = 2;
my $cons_column = 3;
my $seq_column = 8;
my $out_file = "tmp";


GetOptions ("q=s" => \$quartile,
            "o=s" => \$out_file,
            "sample=s" => \$sample,
            "h" => \$h,
            "gold=s" => \$gold,
            "mmr=i" => \$mmRate_column,
            "ref=i" => \$ref_column,
            "cons=i" => \$cons_column,
            "seq=s" => \$seq_column);

die "Author: Shawn Yost <seyost\@ucsd.edu> \n
This program filters variants based on their Q score (Alternate allele frequency - local mismatch rate). The program uses a set of gold standard SNPs’ Q scores to determine the minimum Q score for the inputted variants. The program requires 3 files. The first file is a  variant file with the local-mismatch rates calculated by filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl.  The second file is a gold standard variant file with the local-mismatch rates calculated by filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl. The third file is the output file were unfiltered variants will be printed to.\n\nSee manual for more details\n 
Run: perl filter_on_local_mismatch_Rates.pl [Options]\n
Options:\n
	-h {} = Print this help message\n
	-gold {FILE} = Gold standard set of variants. The file should contain the variants' local-mismatch rates calculated by filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl (REQUIRED)\n
	-sample {FILE} = Variant file to filter. The file should contain the variants' local-mismatch rates calculated by filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl (REQUIRED)\n
	-o {FILE} = Output file (REQUIRED)\n
	-mmr {INT} = 0-based column in the variant file (gold/sample) that stores the local mismatch rate [12]\n
	-seq {INT} = 0-based column in the variant file (gold/sample) that stores the 'read base' information in a samtools pileup format [8]\n
	-ref {INT} = 0-based column in the variant file (gold/sample) that stores the reference base information [2]\n
	-cons {INT} = 0-based column in the variant file (gold/sample) that stores the consensus base and/or the alternate-allele [3]\n
	-q {DOUBLE} = quartile to use on the gold standard variants to calculate the minimum Q score used to filter 'sample' variants [0.10]\n
\n" if($h or $sample eq "tmp" or $gold eq "tmp" or $out_file eq "tmp");







# Open the gold set of variants and calculate the alternate allele frequency.
open(GOLD, $gold);
while(my $line = <GOLD>){
  chomp $line;
  my @fields = split(/\t/, $line);
  my $ref = $fields[$ref_column];
  my $cons = $fields[$cons_column];
  my $seq = $fields[$seq_column];
  my $mmrate = $fields[$mmRate_column];
  $ref =~ tr/[a-z]/[A-Z]/;
  $cons =~ tr/[a-z]/[A-Z]/;
  
  my $alt = convert($ref, $cons);
  
  my ($alt_cov, $total_cov) = get_allele_freq($ref, $alt, $seq);
 
  my $alt_freq = $alt_cov/$total_cov;  
  my $q_score = $alt_freq - $mmrate;
    
  push(@gold_het_qscore, $q_score);
}
close(GOLD);

my @sorted_gold_het_qscore = sort {$a <=> $b} @gold_het_qscore;

#calculate the minimum Q score to filter the variants on.
my $gold_het_length = $#sorted_gold_het_qscore;
my $quartile_index = int($gold_het_length * $quartile);
my $min_qscore = $sorted_gold_het_qscore[$quartile_index];

open(OUT, ">$out_file");

#Read in and filter the 'sample' variants
open(SAMP, $sample);
while(my $line = <SAMP>){
  chomp $line;
  my @fields = split(/\t/, $line);
  my $ref = $fields[$ref_column];
  my $cons = $fields[$cons_column];
  my $seq = $fields[$seq_column];
  my $mmrate = $fields[$mmRate_column];
  
  $ref =~ tr/[a-z]/[A-Z]/;
  $cons =~ tr/[a-z]/[A-Z]/;
  
  
  my $alt = convert($ref, $cons);
  
  my ($alt_cov, $total_cov) = get_allele_freq($ref, $alt, $seq);
  
  my $alt_freq = $alt_cov/$total_cov;  
  my $q_score = $alt_freq - $mmrate;
  
  if($q_score >= $min_qscore){
    print OUT "$line\n";
  }else{
    
  }
  
}
close(SAMP);
close(OUT);








#############################################################
#  Input = The reference and consensus base call
#  Output = The alternate allele for the consensus base call
#  Purpose = To convert from IUPAC Necleic acid coes to A, T,
#							 C, and G
#############################################################
sub convert($$){
	my ($r, $call) = @_;
	
	if($call eq "A" or $call eq "T" or $call eq "C" or $call eq "G"){
		return $call;
	}elsif($call eq "R"){
		if($r eq "A"){
			return "G";
		}else{
			return "A";
		}	
	}elsif($call eq "Y"){
		if($r eq "C"){
			return "T";
		}else{
			return "C";
		}
	}elsif($call eq "S"){
		if($r eq "G"){
			return "C";
		}else{
			return "G";
		}
	}elsif($call eq "W"){
		if($r eq "A"){
			return "T";
		}else{
			return "A";
		}
	}elsif($call eq "K"){
		if($r eq "G"){
			return "T";
		}else{
			return "G";
		}
	}elsif($call eq "M"){
		if($r eq "A"){
			return "C";
		}else{
			return "A";
		}
	}else{
		return "N";
	} 
	
}


#############################################################
#  Input = Reference base, Alternate allele, Pileup format "read bases"
#  Output = Alternate allele coverage and total coverage
#  Purpose = Extracts the alternate allele coverage and read coverage
#							from the 'read base' format from a samtools pileup file
#############################################################
sub get_allele_freq($$$){
	my ($ref, $alt, $read_bases) = @_;
	my %SNPs = ('A',0,'T',0,'C',0,'G',0);

	if($read_bases =~ m/[\$\^\+-]/) {
		$read_bases =~ s/\^.//g; #removing the start of the read segement mark
		$read_bases =~ s/\$//g; #removing end of the read segment mark
		while ($read_bases =~ m/[\+-]{1}(\d+)/g) {
			my $indel_len = $1;
			$read_bases =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
		}
	}

	my @bases = split(//, $read_bases);
	for my $base ( 0 .. @bases - 1 ) {
		if( $bases[ $base ] =~ m/[ATGC]/i ){
			$SNPs{ uc( $bases[ $base ] ) } += 1;
		}elsif( $bases[ $base ] =~ m/[\.,]/ ){
			$SNPs{ uc( $ref ) } += 1;
		}  
	}
	
	my $alt_coverage = $SNPs{$alt};
	my $total_coverage = 0;
	foreach my $b (keys %SNPs){
		$total_coverage += $SNPs{$b};
	}
	
	return($alt_coverage, $total_coverage); 
}
