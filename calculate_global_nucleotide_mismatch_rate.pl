#! /usr/bin/perl -w
###########################
# Author: Shawn Yost <seyost@ucsd.edu>
# 
############################

use warnings;
use strict;
use Getopt::Long;

sub open_vars($\%);
sub open_1000Genomes($\%);
sub open_dbsnp($\%);
sub open_chrom_sizes($\@\%);
sub calculate_random_sites(\%\%$$\@\%\%$);
sub open_pileup($\%\%\%);


my $dbsnp_file;
my $var_file;
my $num_of_sites = 100000;
my $pileup_file;
my $MIN_COV = 3;
my $MAX_COV = 100;
my $out_file;
my $genomes_file;
my $chrom_file;
my $genome_covered = 0;
my $chrom_column = 0;
my $ref_column = 2;
my $position_column = 1;
my $seq_column = 8;
my $cov_column = 7;
my $h;

GetOptions ("v=s" => \$var_file,
						"d=s" => \$dbsnp_file,
						"n=i" => \$num_of_sites,
						"p=s" => \$pileup_file,
						"min=i" => \$MIN_COV,
						"max=i" => \$MAX_COV,
						"o=s" => \$out_file,
						"g=s" => \$genomes_file,
						"s=s" => \$chrom_file,
						"genCov=i" => \$genome_covered,
						"chr=i" => \$chrom_column,
						"ref=i" => \$ref_column,
						"pos=i" => \$position_column,
						"cov=i" => \$cov_column,
						"seq=i" => \$seq_column,
						"h" => \$h);

die "AUTHOR: Shawn Yost <seyost\@ucsd.edu>\n
calculate_global_nucleotide_mismatch_rate.pl calculates the global nucleotide mismatch rates for the given sample. See the manual (manual.pdf) for the definition of the 'global nucleotide mismatch rates'. There are 3 required files as inputs and one required option. The first file is a tab-delaminated file with a list of sites that are non-reference (called variants). The second file is the pileup file for the whole genome of the sample in which the global nucleotide mismatch rates will be calculated for. The third file is the output file that the global nucleotide mismatch rates will be printed to. The last required option is the '-genCov' option that specifies the percent of the genome that has >= 'min' coverage and <= 'max' coverage. It is HIGHLY recommended that you input a dbSNP file and/or a snps from the 1000 genomes project (see manual for file details).\n\nFor more information see the manual\n
To run: calculate_global_nucleotide_mismatch_rate.pl [OPTIONS] -v VAR_FILE -p PILEUP_FILE -o OUT_FILE\n
OPTIONS:
  -h {} = Print help message\n
  -v {FILE} = A list of sites that are called variant (i.e. a list of sites to skip when calculating the global nucleotide substitution rate) (REQUIRED)\n
  -d {FILE} = dbSNP file in the UCSC genome browser format, see manual (RECOMMENDED)\n
  -g {FILE} = 1000 genomes file, see manual for format (RECOMMENDED)\n
  -p {FILE} = Inputted file in a samtools pileup format (REQUIRED)\n
  -o {FILE} = Output file for the global nucleotide substitution rates (REQUIRED)\n
  -s {FILE} = A file containing all chromosomes to examine and the total number of bases for that chromosome, see manual for details (RECOMMENDED) [hg18 chromosomes]  
  -n {INT} = The number of random As, Ts, Cs, and Gs to examine (OPTIOAL) [100000] \n
  -min {INT} = Minimum coverage for bases used to calculate the global nucleotide substitution rate (OPTIONAL) [3]\n
  -max {INT} = Maximum coverage for bases used to calculate the global nucleotide substitution rate (OPTIONAL) [100]\n
  -genCov {INT} = The percent of the genome with >= MIN_COVERAGE and <= MAX_COVERAGE; ranges from 1-100 (REQUIRED)\n
  -chr {INT} = 0-based coordinate of the column with the chromosome information (OPTIONAL) [0]\n
  -pos {INT} = 0-based coordinate of the column with the chromosome position information (OPTIONAL) [1]\n
  -ref {INT} = 0-based coordinate of the column with the reference information (OPTIONAL) [2]\n
  -cov {INT} = 0-based coordinate of the column with the coverage information (OPTIONAL) [7]\n
  -seq {INT} = 0-based coordinate of the column with the samtools pileup 'read base' information (OPTIONAL) [8]\n
\n" if($h or !$out_file or !$var_file or !$pileup_file or $genome_covered == 0);


#$genome_covered = $genome_covered/100;

my %skipped_sites;
my %ref_sites;
my %order;
my %pileup_info;
my %chromosome_sizes;
my @chromosome_array;
my %random_sites;

@chromosome_array = ("chr7","chr7","chr7","chr20","chr22","chr14","chr14","chrY","chr19","chr8","chr8","chr8","chr1","chr1","chr1","chr1","chr1","chr11","chr11","chr11","chr6","chr6","chr6","chr6","chr17","chr17","chr21","chr16","chr16","chr18","chr18","chr3","chr3","chr3","chr3","chr12","chr12","chr12","chr15","chr15","chrX","chrX","chrX","chr4","chr4","chr4","chr4","chr2","chr2","chr2","chr2","chr2","chr9","chr9","chr9","chr13","chr13","chr10","chr10","chr10","chr5","chr5","chr5","chr5") if(!$chrom_file);

%chromosome_sizes = (
  "chr1" => 247249719,
  "chr2" => 242951149,
  "chr3" => 199501827,
  "chr4" => 191273063,
  "chr5" => 180857866,
  "chr6" => 170899992,
  "chr7" => 158821424,
  "chr8" => 146274826,
  "chr9" => 140273252,
  "chr10" => 135374737,
  "chr11" => 134452384,
  "chr12" => 132349534,
  "chr13" => 114142980,
  "chr14" => 106368585,
  "chr15" => 100338915,
  "chr16" => 88827254,
  "chr17" => 78774742,
  "chr18" => 76117153,
  "chr19" => 63811651,
  "chr20" => 62435964,
  "chr21" => 46944323,
  "chr22" => 49691432,
  "chrX" => 154913754,
  "chrY" => 57772954,
) if(!$chrom_file);

my $total_size = open_chrom_sizes($chrom_file, @chromosome_array, %chromosome_sizes) if($chrom_file);
print "Done processing $chrom_file\n"  if($chrom_file);

open_vars($var_file, %skipped_sites);
print "DONE processing $var_file\n";

open_dbsnp($dbsnp_file, %skipped_sites) if($dbsnp_file);
print "Done processing $dbsnp_file\n"  if($dbsnp_file);

open_1000Genomes($genomes_file, %skipped_sites) if($genomes_file);
print "Done processing $genomes_file\n"  if($genomes_file);

calculate_random_sites(%skipped_sites, %order, $num_of_sites, $genome_covered, @chromosome_array, %chromosome_sizes, %random_sites, $total_size);
print "Done Calculating random sites\n";

open_pileup("$pileup_file", %random_sites, %pileup_info, %ref_sites);
print "Done processing $pileup_file\n";


my $a_t_total_coverage = 0;
my $c_g_total_coverage = 0;

my $cgta_total_mismatches = 0;
my $cggc_total_mismatches = 0;
my $cgat_total_mismatches = 0;
my $atta_total_mismatches = 0;
my $atgc_total_mismatches = 0;
my $atcg_total_mismatches = 0;

my %ref_counts;

$ref_counts{"A"} = 0;
$ref_counts{"T"} = 0;
$ref_counts{"C"} = 0;
$ref_counts{"G"} = 0;


ORDERLOOP: foreach my $rank (sort {$a <=> $b} keys %order){

	my ($chr, $pos) = split(/\t/, $order{$rank});
	if(exists($pileup_info{$chr}{$pos})){
		my $reference = $ref_sites{$chr}{$pos};
		
		next ORDERLOOP if($ref_counts{$reference} > $num_of_sites);
		
		$ref_counts{$reference}++;
		
		my $seq = $pileup_info{$chr}{$pos};

		my %SNPs = ('A',0,'T',0,'C',0,'G',0);
		if($seq =~ m/[\$\^\+-]/){
			$seq =~ s/\^.//g; #removing the start of the read segement mark
			$seq =~ s/\$//g; #removing end of the read segment mark
			while ($seq =~ m/[\+-]{1}(\d+)/g) {
				my $indel_len = $1;
				$seq =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
			}
		}
		$seq =~ tr/[a-z]/[A-Z]/;

		my @bases = split(//, $seq);
		foreach my $base (@bases){
			if($base =~ m/[ATGC]/){
				$SNPs{$base}++;
			}else{
				$SNPs{$reference}++;
			}
		}
		
		if($reference eq "A"){
			$a_t_total_coverage += ($SNPs{"A"} + $SNPs{"T"} + $SNPs{"G"} + $SNPs{"C"});
			$atta_total_mismatches += $SNPs{"T"};
			$atgc_total_mismatches += $SNPs{"G"};
			$atcg_total_mismatches += $SNPs{"C"};
			
		}elsif($reference eq "C"){
			$c_g_total_coverage += ($SNPs{"A"} + $SNPs{"T"} + $SNPs{"G"} + $SNPs{"C"});
			$cgta_total_mismatches += $SNPs{"T"};
			$cggc_total_mismatches += $SNPs{"G"};
			$cgat_total_mismatches += $SNPs{"A"};
			
		}elsif($reference eq "T"){
			$a_t_total_coverage += ($SNPs{"A"} + $SNPs{"T"} + $SNPs{"G"} + $SNPs{"C"});
			$atta_total_mismatches += $SNPs{"A"};
			$atgc_total_mismatches += $SNPs{"C"};
			$atcg_total_mismatches += $SNPs{"G"};
			
		}elsif($reference eq "G"){
			$c_g_total_coverage += ($SNPs{"A"} + $SNPs{"T"} + $SNPs{"G"} + $SNPs{"C"});
			$cgta_total_mismatches += $SNPs{"A"};
			$cggc_total_mismatches += $SNPs{"C"};
			$cgat_total_mismatches += $SNPs{"T"};
			
		}			
	}
	
	last ORDERLOOP if($ref_counts{"A"} > $num_of_sites and $ref_counts{"T"} > $num_of_sites and $ref_counts{"G"} > $num_of_sites and $ref_counts{"C"} > $num_of_sites);
	
}

print STDERR "Only found ".$ref_counts{"A"}." number of 'A' reference bases, instead of the specified $num_of_sites\n" if($ref_counts{"A"} < $num_of_sites);
print STDERR "Only found ".$ref_counts{"T"}." number of 'T' reference bases, instead of the specified $num_of_sites\n" if($ref_counts{"T"} < $num_of_sites);
print STDERR "Only found ".$ref_counts{"G"}." number of 'G' reference bases, instead of the specified $num_of_sites\n" if($ref_counts{"G"} < $num_of_sites);
print STDERR "Only found ".$ref_counts{"C"}." number of 'C' reference bases, instead of the specified $num_of_sites\n" if($ref_counts{"C"} < $num_of_sites);


my $cgta_glob_sub_rate = $cgta_total_mismatches/$c_g_total_coverage;
my $cggc_glob_sub_rate = $cggc_total_mismatches/$c_g_total_coverage;
my $cgat_glob_sub_rate = $cgat_total_mismatches/$c_g_total_coverage;
my $atta_glob_sub_rate = $atta_total_mismatches/$a_t_total_coverage;
my $atgc_glob_sub_rate = $atgc_total_mismatches/$a_t_total_coverage;
my $atcg_glob_sub_rate = $atcg_total_mismatches/$a_t_total_coverage;

open(OUT, ">$out_file");

print OUT "ATTA\t$atta_glob_sub_rate\n";
print OUT "ATGC\t$atgc_glob_sub_rate\n";
print OUT "ATCG\t$atcg_glob_sub_rate\n";
print OUT "CGTA\t$cgta_glob_sub_rate\n";
print OUT "CGGC\t$cggc_glob_sub_rate\n";
print OUT "CGAT\t$cgat_glob_sub_rate\n";

close(OUT);







sub open_vars($\%){
	my ($file, $skipped_hash) = @_;
	open(VARS, $file);
	while(my $line = <VARS>){
		chomp $line;
		my @fields = split(/\t/, $line);
		my $chr = $fields[$chrom_column];
		my $pos = $fields[$position_column];
		$$skipped_hash{$chr}{$pos} = 1;
	}
	close(VARS);
}

sub open_1000Genomes($\%){
	my ($file, $skipped_hash) = @_;
	open(GENOMES, $file);
	while(my $line = <GENOMES>){
		chomp $line;
		my @fields = split(/\t/, $line);
		my $chr = "chr".$fields[0];
		my $pos = $fields[1];
		$$skipped_hash{$chr}{$pos} = 1;
	}
	close(GENOMES);	
}

sub open_dbsnp($\%){
	my ($file, $skipped_hash) = @_;
	open(DBSNP, $file);
DBSNPLINE:	while(my $line = <DBSNP>){
		chomp $line;
		my @fields = split(/\t/, $line);
		next DBSNPLINE if($fields[11] ne "single");
		my $chr = $fields[1];
		my $pos = $fields[3];
		$$skipped_hash{$chr}{$pos} = 1;
	}
	close(DBSNP);
}

sub open_chrom_sizes($\@\%){
	my ($chrom_file, $chrom_array, $chrom_hash) = @_;
	my $min_size = 10000000000;
	my $total_size = 0;
	open(CHROM, $chrom_file);
	while(my $line = <CHROM>){
		chomp $line;
		my ($chr, $start, $end) = split(/\t/, $line);
		$$chrom_hash{$chr} = $end;
		$total_size += ($end - $start);
		if($end < $min_size){
			$min_size = $end;
		}
	}
	
	foreach my $c (keys %$chrom_hash){
		my $size = $$chrom_hash{$c};
		my $push_count = int($size/$min_size + 0.5);
		for(my $i = 1; $i <= $push_count; $i++){
			push(@$chrom_array, $c);
		}
	}
	close(CHROM);
	return ($total_size);
}

sub calculate_random_sites(\%\%$$\@\%\%$){
	my ($skipped_hash, $order_hash, $num_of_sites, $genome_covered, $chrom_array, $chrom_hash, $sites_hash, $total_size) = @_;
	
	my $total_num_of_random_sites = (($num_of_sites * 4)*(100/$genome_covered)*2);

FORLOOP:	for(my $i = 1; $i <= $total_num_of_random_sites; $i++){
		my $random_chrom = $$chrom_array[int(rand(scalar(@$chrom_array)))];
		my $random_position = int(rand($$chrom_hash{$random_chrom}))+1;
		
		if(exists($$skipped_hash{$random_chrom}{$random_position}) or exists($$sites_hash{$random_chrom}{$random_position})){
			$i--;
			next FORLOOP;
		}
		
		$$order_hash{$i} = $random_chrom."\t".$random_position;
		$$sites_hash{$random_chrom}{$random_position} = 1;
	}
	
}

sub open_pileup($\%\%\%){
	my ($pileup_file, $sites_hash, $pileup_info_hash, $ref_hash) = @_;
	open(PILEUP, $pileup_file);
PILEUPLINE:	while(my $line = <PILEUP>){
		chomp $line;
		my @fields = split(/\t/, $line);
		my $chr = $fields[$chrom_column];
		my $pos = $fields[$position_column];
		my $seq = $fields[$seq_column];
		my $cov = $fields[$cov_column];
		my $ref = $fields[$ref_column];

		next PILEUPLINE if($ref eq "*");
		next PILEUPLINE if($cov < $MIN_COV or $cov > $MAX_COV);
		next PILEUPLINE if(!($ref =~ m/a|t|c|g/i));
		
		if(exists($$sites_hash{$chr}{$pos})){
			$ref =~ tr/[a-z]/[A-Z]/;
			$$pileup_info_hash{$chr}{$pos} = $seq;
			$$ref_hash{$chr}{$pos} = $ref;
		}
	}
	close(PILEUP);
}



