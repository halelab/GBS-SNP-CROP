#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-2.pl help
###########################################################################################################
############################################################
# Requirement 1: Java 7 or higher - script tested with version 8 (update 221)
# Requirement 2: Trimmomatic (Bolger et al., 2014) - script tested with v0.39
# Requirement 3: cpan module Getopt::Long
############################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

######################
# The help function
######################
my $help = $ARGV[0];
my ($trimmomatic,$dataType,$fastq_seed,$threads,$phred,$adaptor,$leading,$sliding,$trailing,$minlen);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.1\n"
	."Step 2: Trim based on quality\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-2.pl [options]\n\n"
	."Options:\n"
	."-tm: Path to Trimmomatic jar file. String. Default: /usr/local/bin/trimmomatic-0.39.jar.\n"
	."-d: Data type. Either PE (Paired-End) or SE (Single-End). String. Required.\n"
	."-fq: FASTQ file name seed (convention: FileNameSeed_001.R1parsed.fq.gz). String. Required\n"
	."-t: Number of independent threads used. Numeric. Default: 10\n"
	."-ph: Trimmomatic Phred scale. Either 33 or 64. Numeric. Default: 33\n"
	."-ad: Trimmomatic ILLUMINACLIP parameter (FASTA file with adaptors:Seed mismatches:Palindrome clip threshold:Simple clip threshold). String. Default: TruSeq3-PE.fa:2:30:10\n"
	."-l: Trimmomatic LEADING parameter value. Numeric. Default: 30\n"
	."-sl: Trimmomatic SLIDINGWINDOW parameter (Window size:Required quality). Colon-separated numeric. Default: 4:30\n"
	."-tr: Trimmomatic TRAILING parameter value. Numeric. Default: 30\n"
	."-m: Trimmomatic MINLEN parameter value. Numeric. Default: 32\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameter values
#################################
$trimmomatic = '/usr/local/bin/trimmomatic-0.39.jar';
$threads = 10;	$phred = 33;
$adaptor = 'TruSeq3-PE.fa:2:30:10';
$leading = 30;	$sliding = '4:30';
$trailing = 30;	$minlen = 32;

GetOptions(
'tm=s' => \$trimmomatic,  # string
'd=s' => \$dataType,      # string
'fq=s' => \$fastq_seed,   # string
't=s' => \$threads,       # numeric
'ph=s' => \$phred,        # numeric
'ad=s' => \$adaptor,      # string 
'l=s' => \$leading,       # numeric
'sl=s' => \$sliding,      # numeric
'tr=s' => \$trailing,     # numeric
'm=s' => \$minlen,        # numeric
) or die "$H";

print "\n#################################\n# GBS-SNP-CROP, Step 2, v.4.1\n#################################\n";
my $sttime = time;

############################
# Trimming Paired-End data 
############################

if ($dataType eq "PE") {
	print "Trimming paired-end reads ...\n";

	print "\nConcatenating library of $fastq_seed R1 reads... ";
	my $PER1IN = join ("","$fastq_seed",".R1parsed.merged.fq.gz");
	system ( "cat $fastq_seed*R1parsed.fq.gz > $PER1IN" );
	print "DONE.";

	print "\nConcatenating library of $fastq_seed R2 reads... ";
	my $PER2IN = join ("","$fastq_seed",".R2parsed.merged.fq.gz");
	system ( "cat $fastq_seed*R2parsed.fq.gz > $PER2IN" );
	print "DONE.";

	print "\nInitiating Trimmomatic paired-end read trimming...\n\n";

	my $PER1OUT = join("","$fastq_seed",".PE.R1parsed.merged.trimmed.fq.gz");
	my $SER1OUT = join("","$fastq_seed",".SE.R1parsed.merged.trimmed.fq.gz");
	my $PER2OUT = join("","$fastq_seed",".PE.R2parsed.merged.trimmed.fq.gz");
	my $SER2OUT = join("","$fastq_seed",".SE.R2parsed.merged.trimmed.fq.gz");
	
	if ( $adaptor ne '0' ) {
		system ( "java -jar $trimmomatic PE -threads $threads -phred$phred $PER1IN $PER2IN $PER1OUT $SER1OUT $PER2OUT $SER2OUT ILLUMINACLIP:$adaptor SLIDINGWINDOW:$sliding LEADING:$leading TRAILING:$trailing MINLEN:$minlen" );
	} elsif ( $adaptor eq '0' ) {
		system ( "java -jar $trimmomatic PE -threads $threads -phred$phred $PER1IN $PER2IN $PER1OUT $SER1OUT $PER2OUT $SER2OUT SLIDINGWINDOW:$sliding LEADING:$leading TRAILING:$trailing MINLEN:$minlen" );
	}

	system ( "rm *parsed.fq.gz" );
	print "\nDONE.\n";
	
############################
# Trimming Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Trimming single-end reads ...\n";

	print "\nConcatenating library of $fastq_seed reads... ";
	my $SER1IN = join ("","$fastq_seed",".R1parsed.merged.fq.gz");
	system ( "cat $fastq_seed*R1parsed.fq.gz > $SER1IN" );
	print "DONE.";

	print "\nInitiating Trimmomatic single-end read trimming...\n\n";

	my $SER1OUT = join("","$fastq_seed",".SE.R1parsed.merged.trimmed.fq.gz");
	
	if ( $adaptor ne '0' ) {
		system ( "java -jar $trimmomatic SE -threads $threads -phred$phred $SER1IN $SER1OUT ILLUMINACLIP:$adaptor SLIDINGWINDOW:$sliding LEADING:$leading TRAILING:$trailing MINLEN:$minlen" );
	} elsif ( $adaptor eq '0' ) {
		system ( "java -jar $trimmomatic SE -threads $threads -phred$phred $SER1IN $SER1OUT SLIDINGWINDOW:$sliding LEADING:$leading TRAILING:$trailing MINLEN:$minlen" );
	}

	system ( "rm *parsed.fq.gz" );
	print "\nDONE.\n";
}
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;