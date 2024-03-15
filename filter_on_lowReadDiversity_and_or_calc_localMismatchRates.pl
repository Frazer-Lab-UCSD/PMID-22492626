#! /usr/bin/perl -w
###########################
# Author: Shawn Yost <seyost@ucsd.edu>
# 
############################

use warnings;
use strict;
use Getopt::Long;

sub get_reads($$\%);
sub convert($$);


my $h;
my $variant_file = "tmp";  #Assumes it is TAB seperated
my $BAM_file = "tmp";
my $TMP_file = "tmp";
my $OUT_PREFIX = "tmp";
my $chromosome = 0;
my $position = 1;
my $reference_base = 2;
my $consensus_base = 3;
my $remove_sameStart_only;
my $calculate_mismatch_rate_only;
my $window = 10;
my $minimum_start_locations = 3;


GetOptions("h" => \$h,
						"v=s" => \$variant_file,
						"b=s" => \$BAM_file,
						"t=s" => \$TMP_file,
						"o=s" => \$OUT_PREFIX,
						"chr=i" => \$chromosome,
						"pos=i" => \$position,
						"ref=i" => $reference_base,
						"cons=i" => $consensus_base,
						"window=i" => $window,
						"rso" => \$remove_sameStart_only,
						"cmro" => \$calculate_mismatch_rate_only,
						"minstart=i" => \$minimum_start_locations);

die "<AUTHOR: Shawn Yost <seyost\@ucsd.edu> \n
This program implements two different methods. The first method removes variants with low read diversity. This program removes variants with biased read diversity: Duplicate sequencing reads carrying an error can result in false positive calls. The program removes candidate variants supported by reads with less than 'minstart' different start positions. The second method this program implements is to calculate the local-mismatch rate for a given variant. The local mismatch rate (LMR) = m/(n+m), where (m) is the number of positions matching the reference and (n) the number of mismatched (excluding the candidate variant itself). The LMR is calculated for a given 'window' around the variant. The user can specify to run Both methods or only one method (default is to run both methods).\n\nSee manual for more help\n
Run: perl filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl [Options]\n
Options:\n
	-h {} = Print out help message\n
	-v {FILE} = Variant file (REQUIRED)\n
	-b {FILE} = Sorted-indexed BAM file that has been processed with SAMTools calmd -e (REQUIRED)\n
	-t {FILE} = Temporary file (REQUIRED)\n
	-o {FILE_PREFIX} = Output file prefix (REQUIRED)\n
	-chr {INT} = 0-based column in the variant file (v) that stores the chromosome information [0]\n
	-pos {INT} = 0-based column in the variant file (v) that stores the variant position information [1]\n
	-ref {INT} = 0-based column in the variant file (v) that stores the reference base information [2]\n
	-cons {INT} = 0-based column in the variant file (v) that stores the consensus base and/or the alternate-allele [3]\n
	-window {INT} = +/- Window size used to determine the misMatch rates of a variant [10]\n
	-minstart {INT} = Minimum number of reads with different starting alignment locations needed in order to accept a variant [3]\n
	-rso {} = Specifies to ONLY remove same-start-site variants (do NOT calculate misMatch rate)\n
	-cmro {} = Specifies to ONLY calculate the misMatch rate for all variants (do NOT filter variants)\n
\n" if($h or $variant_file eq "tmp" or $BAM_file eq "tmp" or $TMP_file eq "tmp" or $OUT_PREFIX eq "tmp");

open(PASS, ">$OUT_PREFIX\_sameStart\_mmRate.txt") if(!($calculate_mismatch_rate_only) and !($remove_sameStart_only));
open(FILT, ">$OUT_PREFIX\_filtered.txt") if(!($calculate_mismatch_rate_only));
open(MMR, ">$OUT_PREFIX\_mmRate.txt") if($calculate_mismatch_rate_only);
open(SSR, ">$OUT_PREFIX\_sameStart.txt") if($remove_sameStart_only);


open(VARIANTS, $variant_file);
READLINE: while(my $line = <VARIANTS>){
	chomp $line;
	my @fields = split(/\t/, $line);
	my $chr = $fields[$chromosome];
	my $pos = $fields[$position];
	my $ref = $fields[$reference_base];
	my $cons = $fields[$consensus_base];
	
	$ref =~ tr/[a-z]/[A-Z]/;
	if($ref ne "A" and $ref ne "T" and $ref ne "G" and $ref ne "C"){ #skip the variant if it is an indel
		print "SKIPPED: Indels are not supported\n$line\n";
		next READLINE;
	}

	$cons =~ tr/[a-z]/[A-Z]/;
	
	my $alt = convert($ref, $cons); #optain the alternate allele
	
	if($alt eq "N"){ #The script doesn't support tri-allelic bases
		print "SKIPPED: Currently this program does not support tri-allelic bases\n";
		print "$line\n";
		next READLINE;
	}

	#Make the window size
	my $end = $pos + $window;
	my $begin = $pos - $window;
	
	my %reads = ();
	get_reads($chr, $pos, %reads); #obtain the reads for the Window
	
	my %starts = (); #start sites for a given variant
	my $number_of_mismatches=0;
	my $number_of_matches=0;
	my $total=0;

	#loop through the reads
	foreach my $start (sort {$a <=> $b} keys %reads){
		my $alternate_alignment_position = $pos-$start;

			#loop through the reads with different cigar strings and alignments
READ:		foreach my $cig_align (keys %{$reads{$start}}){
			my $strand = $reads{$start}{$cig_align};
			my $cur_pos = $start;
			my ($cigar, $alignment) = split(/\t/, $cig_align);
			my @alignment_array = split(//, $alignment);
			my @cigar_counts = split(/\D+/, $cigar);
			my @cigar_types = split(/\d{1,2}/, $cigar);
			shift @cigar_types;

			my $current_alignment_position=0;
			my $cigar_string_count= 0;
			my $current_alignment_call;
			my $current_type;
			
			#Loop through the Cigar string to determin which sites are matches, deletions, or insertions
			for(my $j = 0; $j <= $#cigar_counts; $j++){
				#Loop through the alignment one base at a time using the CIGAR string
				for(my $i = 1; $i <= $cigar_counts[$j]; $i++){
					$current_type = $cigar_types[$j];
					
					#Skip deletions
					if($current_type eq "D"){
						$current_alignment_call = "D";
					}else{
						$current_alignment_call = $alignment_array[$current_alignment_position];
						$current_alignment_position++;
					}
					
###########################################
# Determine if the start site of the read is uniqu or already exists
###########################################
					if(!($calculate_mismatch_rate_only)){
						#If the current position is the alternate allele position
						if($alternate_alignment_position == $cigar_string_count){
							
							#If it is a match and the alternate allele matches the current alignment for that read then count a new start site
							if($current_type eq "M" and $alt eq $current_alignment_call){
								if($strand == 0){
									$starts{$start} = 1;
								}else{
									my $read_end = $start;
									for(my $k = 0; $k <= $#cigar_counts;$k++){
										if($cigar_types[$k] ne "I"){
											$read_end += $cigar_counts[$k];
										}
									}
									$starts{$read_end} = 1;
								}
							}
						}
					}

###########################################
# Calculate the mismatch rate
###########################################
					if(!($remove_sameStart_only)){
						#If the alignment position is within the specified window determine if it is a match or mismatch
						if(($cigar_string_count + $start) >= $begin and ($cigar_string_count + $start) <= $end){
							
							#If the alignment isn't an '=' and isn't the alternate allele for the given variant then count it as a mismatch
							if($current_alignment_call ne "=" and ($alternate_alignment_position != $cigar_string_count or ($alternate_alignment_position == $cigar_string_count and $current_alignment_call ne $alt))){
								$number_of_mismatches++;
							
								#else if it matches the reference count a match
							}elsif($current_alignment_call eq "="){
								$number_of_matches++;
							}
						}
					}
					
					
					
					$cigar_string_count++;
					next READ if(($cigar_string_count + $start) > $end);
				}
			}
		}
	}

	#determine if the variant has more then or equal to the minimum number of starting locations specified
	if(scalar(keys(%starts)) < $minimum_start_locations and (!($calculate_mismatch_rate_only))){
		print FILT "$line\n";
		next READLINE;
	}
	#print out and append the number of mismatches for the variant that passed the above filter.
	$total = $number_of_mismatches + $number_of_matches;
	my $mismatch_rate = $number_of_mismatches/$total;
	print PASS "$line\t$number_of_mismatches\t$number_of_matches\t$mismatch_rate\n" if(!($calculate_mismatch_rate_only) and !($remove_sameStart_only));
	print MMR "$line\t$number_of_mismatches\t$number_of_matches\t$mismatch_rate\n" if($calculate_mismatch_rate_only);
	print SSR "$line\n" if($remove_sameStart_only);
	
}
close(VARIANTS);

close(PASS) if(!($calculate_mismatch_rate_only) and !($remove_sameStart_only));
close(FILT) if(!($calculate_mismatch_rate_only));
close(MMR) if($calculate_mismatch_rate_only);
close(SSR) if($remove_sameStart_only);

`rm $TMP_file`;


#############################################################
#  Input = The chromosome and position of the variant and a 
#						HASH file to store the reads
#  Output = Fills in the HASH file with the read information
#  Purpose = To get the starting position, cigar string, 
#						alignment string, and strand of all reads aligning 
#						within +/- 10bp of the variant at a given 
#						chromosome and position.
#############################################################
sub get_reads($$\%){
	my ($chr, $pos, $hash) = @_;
	my $s = $pos - $window;
	my $end = $pos + $window;
	`samtools view $BAM_file $chr:$s-$end > $TMP_file`;
	open(TMP, "$TMP_file");
RLINE:	while(my $tmp_line = <TMP>){
		chomp $tmp_line;
		next RLINE if($tmp_line =~ m/^@/);
		my @fields = split(/\t/, $tmp_line);
		my $strand = $fields[1];
		my $start = $fields[3];
		my $cigar = $fields[5];
		my $alignment = $fields[9];

		$$hash{$start}{$cigar."\t".$alignment}=$strand;
	}
	close(TMP);
}

#############################################################
#  Input = The reference and consensus base call
#  Output = The alternate allele for the consensus base call
#  Purpose = To convert from IUPAC Necleic acid coes to A, T,
#							 C, and G
#############################################################
sub convert($$){
	my ($reference, $consensus_call) = @_; 
	
	if($consensus_call eq "A" or $consensus_call eq "T" or $consensus_call eq "C" or $consensus_call eq "G"){
		return $consensus_call;
	}elsif($consensus_call eq "R"){
		if($reference eq "A"){
			return "G";
		}else{
			return "A";
		}	
	}elsif($consensus_call eq "Y"){
		if($reference eq "C"){
			return "T";
		}else{
			return "C";
		}
	}elsif($consensus_call eq "S"){
		if($reference eq "G"){
			return "C";
		}else{
			return "G";
		}
	}elsif($consensus_call eq "W"){
		if($reference eq "A"){
			return "T";
		}else{
			return "A";
		}
	}elsif($consensus_call eq "K"){
		if($reference eq "G"){
			return "T";
		}else{
			return "G";
		}
	}elsif($consensus_call eq "M"){
		if($reference eq "A"){
			return "C";
		}else{
			return "A";
		}
	}else{
		return "N"; #Currently the program does not support bases that are tri-allelic
	} 
}

