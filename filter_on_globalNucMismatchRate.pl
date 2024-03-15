#! /usr/bin/perl -w
###########################
# Author: Shawn Yost <seyost@ucsd.edu>
# 
############################

use warnings;
use strict;
use Getopt::Long;

sub split_variants($\%\%);
sub convert($$);
sub get_allele_freq($$$);


my $chrom_column = 0;
my $pos_column = 1;
my $ref_column = 2;
my $cons_column = 3;
my $seq_column = 8;
my $var_file;
my $h;
my $global_subs_file;
my $tmp_file;
my $alpha = 0.05;
my $total_possible_variant_locations = 0;
my $out_file;
my $PATH;

my %variants;
my %var_out;
my %global_substitution_rates;
my %pvals;

my @pvals_array;

GetOptions ("chr=i" => \$chrom_column,
            "pos=i" => \$pos_column,
            "v=s" => \$var_file,
            "h" => \$h,
            "cons=i" => \$cons_column,
            "ref=i" => \$ref_column,
            "seq=i" => \$seq_column,
            "gnmr=s" => \$global_subs_file,
            "t=s" => \$tmp_file,
            "a=s" => \$alpha,
            "tot=i" => \$total_possible_variant_locations,
            "o=s" => \$out_file,
            "p=s" => \$PATH);

die "Author: Shawn Yost <seyost\@ucsd.edu> \n
Run: perl filter_on_globalNucMismatchRate.pl [Options]\n
This program calculates a false discovery rate according to the Benjemin-Hotchberg method and filters variants below a specified alpha. There are three required input files and one required option. The variant file is in a SAMTools pileup consensus format (FILE TYPE 1). The second file is the global nucleotide mismatch rates as calculated by 'calculate_golobal_nucleotide_mismatch_rate.pl'. The third file is a temp file that will be deleted when the program finishes. The last required option is the total possible locations in the genome in which a variant can be called (-tot). This option specifies the 'm' parameter in the Benjemin-Hotchberg method (http://en.wikipedia.org/wiki/False_discovery_rate).\n\nFor more information see the manual (manual.pdf)\n\n 
Options:\n
	-h {} = Print out help message\n
	-gnmr {FILE} = Global nucleotide mismatch rate file (REQUIRED)\n
	-v {FILE} = Variant file (REQUIRED)\n
	-o {FILE} = Output file (REQUIRED)\n
	-t {FILE} = tmp file (REQUIRED)\n
	-chr {INT} = 0-based column in the variant file that stores the chromosome information [0]\n
	-seq {INT} = 0-based column in the variant file that stores the 'read base' information in a samtools pileup format [8]\n
	-ref {INT} = 0-based column in the variant file that stores the reference base information [2]\n
	-cons {INT} = 0-based column in the variant file that stores the consensus base and/or the alternate-allele [3]\n
	-pos {INT} = 0-based column in the variant file that stores the genomic coordinate information [0]\n
	-a {DOUBLE} = alpha to use in the Benjemin-Hotchberg method, also refers to the FDR [0.05]\n
	-tot {INT} = total possible locations in the genome in which you could possibly call a somatic variant (REQUIRED)\n
	-p {PATH} = The path to the 'calc_binom_pval.r' file.  This file is in the same folder as the filter_on_globalNucMismatchRate.pl program\n
\n" if($h or !$var_file or !$out_file or !$global_subs_file or $total_possible_variant_locations == 0 or !$PATH);


split_variants($var_file, %variants, %var_out);

open(GLOB, $global_subs_file);
while(my $line = <GLOB>){
	chomp $line;
	my ($subs_type, $gsr) = split(/\t/, $line);
	$global_substitution_rates{$subs_type} = $gsr;
}
close(GLOB);

open(TMP, ">$tmp_file");
foreach my $substitution_type (sort keys %global_substitution_rates){
	my $gsr = $global_substitution_rates{$substitution_type};
	foreach my $chr (sort keys %{$variants{$substitution_type}}){
		foreach my $pos (sort {$a <=> $b} keys %{$variants{$substitution_type}{$chr}}){
			my @fields = split(/\t/, $variants{$substitution_type}{$chr}{$pos});
			my $ref = $fields[$ref_column];
			my $cons = $fields[$cons_column];
			my $seq = $fields[$seq_column];
			$ref =~ tr/[a-z]/[A-Z]/;
			$cons =~ tr/[a-z]/[A-Z]/;
			my $alt = convert($ref, $cons);
			my ($alt_coverage, $total_coverage) = get_allele_freq($ref, $alt, $seq);
			print TMP "$chr\t$pos\t$alt_coverage\t$total_coverage\t$gsr\n";
		}
	}
}
close(TMP);

#fix this so they input path or something...
`R --slave --args $tmp_file < $PATH/calc_binom_pval.r`;

my $max_rank = 0;

open(TMPIN, "$tmp_file");
while(my $line = <TMPIN>){
	chomp $line;
	my ($chr, $pos, $alt_cov, $tot_cov, $gsr, $p_value) = split(/\s/, $line);
	$pvals{$p_value}{$chr}{$pos} = 1;
	$max_rank++;
	push(@pvals_array, $p_value);
}
close(TMPIN);

my @sorted_pvalues = sort {$a <=> $b} @pvals_array;

open(OUT, ">$out_file");

my $rank = 1;
LOOP: while($sorted_pvalues[$rank-1] <= ($rank * $alpha)/$total_possible_variant_locations){
	foreach my $chr (sort keys %{$pvals{$sorted_pvalues[$rank-1]}}){
		foreach my $pos (sort {$a <=> $b} keys %{$pvals{$sorted_pvalues[$rank-1]}{$chr}}){
			print OUT "$var_out{$chr}{$pos}\n";
			$rank++;
			last LOOP if($rank > $max_rank or $sorted_pvalues[$rank-1] > ($rank * $alpha)/$total_possible_variant_locations);
		}
	}
}
close(OUT);

`rm $tmp_file`;







sub split_variants($\%\%){
	my ($file, $hash, $hash2) = @_;
	open(VARS, $file);
	while(my $line = <VARS>){
		chomp $line;
		my @fields = split(/\t/, $line);
		my $chr = $fields[$chrom_column];
		my $pos = $fields[$pos_column];
		my $ref = $fields[$ref_column];
		my $cons = $fields[$cons_column];

		$ref =~ tr/[a-z]/[A-Z]/;
		$cons =~ tr/[a-z]/[A-Z]/;
		
		$$hash2{$chr}{$pos} = $line;
		
		if($cons =~ m/A|T|C|G/i){
			if($ref eq "C" and $cons eq "T"){
				$$hash{"CGTA"}{$chr}{$pos} = $line;
			}elsif($ref eq "G" and $cons eq "A"){
				$$hash{"CGTA"}{$chr}{$pos}=$line;
			}elsif($ref eq "C" and $cons eq "G"){
				$$hash{"CGGC"}{$chr}{$pos}=$line;
			}elsif($ref eq "G" and $cons eq "C"){
				$$hash{"CGGC"}{$chr}{$pos}=$line;
			}elsif($ref eq "C" and $cons eq "A"){
				$$hash{"CGAT"}{$chr}{$pos}=$line;
			}elsif($ref eq "G" and $cons eq "T"){
				$$hash{"CGAT"}{$chr}{$pos}=$line;
			}elsif($ref eq "A" and $cons eq "T"){
				$$hash{"ATTA"}{$chr}{$pos}=$line;
			}elsif($ref eq "T" and $cons eq "A"){
				$$hash{"ATTA"}{$chr}{$pos}=$line;
			}elsif($ref eq "A" and $cons eq "G"){
				$$hash{"ATGC"}{$chr}{$pos}=$line;
			}elsif($ref eq "T" and $cons eq "C"){
				$$hash{"ATGC"}{$chr}{$pos}=$line;
			}elsif($ref eq "A" and $cons eq "C"){
				$$hash{"ATCG"}{$chr}{$pos}=$line;
			}elsif($ref eq "T" and $cons eq "G"){
				$$hash{"ATCG"}{$chr}{$pos}=$line;
			}
		}else{
			if($cons eq "R"){
				if($ref eq "A"){
					$$hash{"ATGC"}{$chr}{$pos}=$line;
				}elsif($ref eq "G"){
					$$hash{"CGTA"}{$chr}{$pos}=$line;
				}
			}elsif($cons eq "Y"){
				if($ref eq "C"){
					$$hash{"CGTA"}{$chr}{$pos}=$line;
				}elsif($ref eq "T"){
					$$hash{"ATGC"}{$chr}{$pos}=$line;
				}
			}elsif($cons eq "M"){
				if($ref eq "A"){
					$$hash{"ATCG"}{$chr}{$pos}=$line;
				}elsif($ref eq "C"){
					$$hash{"CGAT"}{$chr}{$pos}=$line;
				}
			}elsif($cons eq "K"){
				if($ref eq "T"){
					$$hash{"ATCG"}{$chr}{$pos}=$line;
				}elsif($ref eq "G"){
					$$hash{"CGAT"}{$chr}{$pos}=$line;
				}
			}elsif($cons eq "W"){
				if($ref eq "A"){
					$$hash{"ATTA"}{$chr}{$pos}=$line;
				}elsif($ref eq "T"){
					$$hash{"ATTA"}{$chr}{$pos}=$line;
				}
			}elsif($cons eq "S"){
				if($ref eq "C"){
					$$hash{"CGGC"}{$chr}{$pos}=$line;
				}elsif($ref eq "G"){
					$$hash{"CGGC"}{$chr}{$pos}=$line;
				}
			}else{
				print "Skipped: $line\nDoes not support Tri-allelic bases\n";
			}	
		}
	}
	close(VARS);	
}



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

