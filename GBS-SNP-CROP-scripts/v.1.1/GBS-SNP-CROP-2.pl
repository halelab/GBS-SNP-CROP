#!/usr/bin/perl

###############################################################################################################################
# GBS-SNP-CROP, Step 2. For description, please see Melo et al. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.
###############################################################################################################################

##########################################################################################
# Requirement 1: Java 7 or high
# Requirement 2: Trimmomatic software (Bolger et al., 2014)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-2.pl -d <data type, PE = Paired-End or SR = Single-End> -fq <FASTQ file name seed> -t <number of threads> -ph <Phred quality score: 33 or 64>\n" 
."-l <Trimmomatic LEADING value> -sl <Trimmomatic SLIDINGWINDOW value> -tr <Trimmomatic TRAILING value> -m <Trimmomatic MINLEN value>.\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n";

my ($dataType,$fastq_seed,$threads,$phred,$leading,$sliding,$trailing,$minlen);

GetOptions(
'd=s' => \$dataType,      # string - "PE" or "SE"
'fq=s' => \$fastq_seed,   # string
't=s' => \$threads,       # numeric
'ph=s' => \$phred,        # numeric
'l=s' => \$leading,       # numeric
'sl=s' => \$sliding,      # numeric
'tr=s' => \$trailing,     # numeric
'm=s' => \$minlen,        # numeric
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 2, v.1.1\n#################################\n";

############################
# Trimming Paired-End data 
############################

if ($dataType eq "PE") {
	print "Trimming Paired-End reads...\n";

	print "\nConcatenating library $fastq_seed R1 reads...";
	my $outR1 = join ("","$fastq_seed","_R1parsed",".merged",".fastq",".gz");
	system ( "cat $fastq_seed*R1parsed.fastq.gz > $outR1" );
	print "DONE.";

	print "\nConcatenating library $fastq_seed R2 reads...";
	my $outR2 = join ("","$fastq_seed","_R2parsed",".merged",".fastq",".gz");
	system ( "cat $fastq_seed*R2parsed.fastq.gz > $outR2" );
	print "DONE.";

	print "\nTrimmomatic paired-end reads cleanning/trimming initiated...\n";

	my $trimoPER1OUT = join("","$fastq_seed","_PE_R1parsed",".fastq");
	my $trimoSER1OUT = join("","$fastq_seed","_SE_R1parsed",".fastq");
	my $trimoPER2OUT = join("","$fastq_seed","_PE_R2parsed",".fastq");
	my $trimoSER2OUT = join("","$fastq_seed","_SE_R2parsed",".fastq");

	system ( "java -jar /usr/local/bin/trimmomatic-0.33.jar PE -phred$phred -threads $threads $outR1 $outR2 $trimoPER1OUT $trimoSER1OUT $trimoPER2OUT $trimoSER2OUT LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );

	print "\nTrimmomatic paired-end reads cleaning/trimming completed.\n";
	
############################
# Trimming Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Trimming Single-End reads...\n";

	print "\nConcatenating library $fastq_seed reads...";
	my $outR1 = join ("","$fastq_seed","_R1parsed",".merged",".fastq",".gz");
	system ( "cat $fastq_seed*R1parsed.fastq.gz > $outR1" );
	print "DONE.";

	print "\nTrimmomatic single-end reads cleanning/trimming initiated...\n";

	my $trimoSER1OUT = join("","$fastq_seed","_SE_R1parsed",".fastq");
	
	system ( "java -jar /usr/local/bin/trimmomatic-0.33.jar SE -phred$phred -threads $threads $outR1 $trimoSER1OUT LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );

	print "\nTrimmomatic single-end reads cleaning/trimming completed.\n";
}

print "\n\nPlease cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
