#!/bin/sh

#########################################################
# AUTHOR: Shawn Yost <seyost@ucsd.edu>
#
# This is a shel script to test the 4 programs made to 
#     remove false positives caused by FFPE DNA damage.
#
# To run type "./test_scripts.sh"
# The script must be run in the directory containing the 
#     perl programs
#
# See the manual.pdf for more details on file formats and
#   help using the programs.
#
#########################################################

echo "#######  Calculate Global Nucleotide Mismatch Rate ########"
echo "#"
echo "# The test genome is chr1:1-1,000,100"
echo "# Only 33% of the test genome is covered there -genCov is 33"
echo "# We only want to use 10,000 random sites (-n) to calculate the GNMR"
echo "# I am imputing both a 1000 genomes and dbsnp file."
echo "# I am running: perl calculate_global_nucleotide_mismatch_rate.pl -v test_files/potential_somatic_variants_chr1_1_1000000.final -d test_files/test_dbsnp.txt -g test_files/test_1000Genomes.txt -p test_files/test_pileup_file.raw -s test_files/test_genome_size.txt -o test_files/test_global_nucletoide_mismatch_rates.txt -genCov 33 -n 10000 "
echo "#"
echo "#"
echo ""

perl calculate_global_nucleotide_mismatch_rate.pl -v test_files/potential_somatic_variants_chr1_1_1000000.final -d test_files/test_dbsnp.txt -g test_files/test_1000Genomes.txt -p test_files/test_pileup_file.raw -s test_files/test_genome_size.txt -o test_files/test_global_nucletoide_mismatch_rates.txt -genCov 33 -n 10000

echo ""
echo ""
echo "Using samtools calmd and indexing BAM"

samtools calmd -e test_files/test_bam.bam test_files/test_chr1_fasta.fa | samtools view -Sbh - > test_files/test_bam_calMD.bam
samtools index test_files/test_bam_calMD.bam

echo "####### Run first filter and calculate local mismatch rate ######"
echo "# "
echo "# Running: perl filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl -v test_files/potential_somatic_variants_chr1_1_1000000.final -b test_files/test_bam_calMD.bam -t tmp.txt -o test_files/potential_somatic_variants_chr1_1_1000000 "
echo "#" 
echo ""

perl filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl -v test_files/potential_somatic_variants_chr1_1_1000000.final -b test_files/test_bam_calMD.bam -t tmp.txt -o test_files/potential_somatic_variants_chr1_1_1000000


echo "############ Calculate local mismatch rate for Gold set of variants #####"
echo "#"
echo "# Running: perl filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl -v test_files/test_gold_set_variants.final -b test_files/test_bam_calMD.bam -t tmp.txt -o test_files/test_gold_set_variants -cmro"
echo "#"
echo "#"

perl filter_on_lowReadDiversity_and_or_calc_localMismatchRates.pl -v test_files/test_gold_set_variants.final -b test_files/test_bam_calMD.bam -t tmp.txt -o test_files/test_gold_set_variants -cmro


echo "############  Filter on LMR ############"
echo "# "
echo "# Running: perl filter_on_local_mismatch_Rates.pl -gold test_files/test_gold_set_variants_mmRate.txt -sample test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRate.txt -o test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered.txt "
echo "# "
echo ""

perl filter_on_local_mismatch_Rates.pl -gold test_files/test_gold_set_variants_mmRate.txt -sample test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRate.txt -o test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered.txt 


echo "############ Filter on GNMR  ###########"
echo "#"
echo "# Since there are only ~334,000 sites I can call a variant at"
echo "# I set the -tot option to 334,000.  Also note that for the "
echo "# '-p' option I set the path to the current directory since"
echo "# that is the directory that contains the calc_binom_pval.r file"
echo "#"
echo "# Running: perl filter_on_globalNucMismatchRate.pl -gnmr test_files/test_global_nucletoide_mismatch_rates.txt -v test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered.txt -t tmp.txt -o test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered_GNMR_filtered.txt -tot 334000 -p $(pwd)"
echo "#"
echo ""

perl filter_on_globalNucMismatchRate.pl -gnmr test_files/test_global_nucletoide_mismatch_rates.txt -v test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered.txt -t tmp.txt -o test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered_GNMR_filtered.txt -tot 334000 -p $(pwd)


echo "The final variant file should only contain a single variant"

more test_files/test_final_output.txt

echo ""
echo "If the line below does not match the line above then the script did not run properly"

more test_files/potential_somatic_variants_chr1_1_1000000_sameStart_mmRateFiltered_GNMR_filtered.txt

