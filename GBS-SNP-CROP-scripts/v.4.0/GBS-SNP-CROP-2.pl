#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-2.pl
###########################################################################################################
############################################################
# Requirement 1: Java 7 or high
# Requirement 2: Trimmomatic software (Bolger et al., 2014)
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
	."Version: 4.0\n"
	."Step 2: Trim based on quality\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-2.pl [options]\n\n"
	."Options:\n"
	."-tm: Path to Trimmomatic jar file. String. Default: /usr/local/bin/trimmomatic-0.33.jar.\n"
	."-d: Data type. Either PE (Paired-End) or SE (Single-End). String. Required.\n"
	."-fq: FASTQ file name seed. String. Required\n"
	."-t: Number of independent threads used. Numeric. Default: 10\n"
	."-ph: Trimmomatic Phred scale. Either 33 or 64. Numeric. Default: 33\n"
	."-ad: Trimmomatic ILLUMINACLIP adaptor. String. Default: TruSeq3-PE.fa:2:30:10\n"
	."-l: Trimmomatic LEADING parameter value. Numeric. Default: 30\n"
	."-sl: Trimmomatic SLIDINGWINDOW parameter values. Colon-separated numeric. Default: 4:30\n"
	."-tr: Trimmomatic TRAILING parameter value. Numeric. Default: 30\n"
	."-m: Trimmomatic MINLEN parameter value. Numeric. Default: 32\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameters values
#################################
$trimmomatic = '/usr/local/bin/trimmomatic-0.33.jar';
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

print "\n#################################\n# GBS-SNP-CROP, Step 2, v.4.0\n#################################\n";
my $sttime = time;

############################
# Trimming Paired-End data 
############################

if ($dataType eq "PE") {
	print "Trimming Paired-End reads ...\n";

	print "\nConcatenating library $fastq_seed R1 reads ... ";
	my $outR1 = join ("","$fastq_seed","_R1parsed",".merged",".fq",".gz");
	system ( "cat $fastq_seed*R1parsed.fq.gz > $outR1" );
	print "DONE.";

	print "\nConcatenating library $fastq_seed R2 reads ... ";
	my $outR2 = join ("","$fastq_seed","_R2parsed",".merged",".fq",".gz");
	system ( "cat $fastq_seed*R2parsed.fq.gz > $outR2" );
	print "DONE.";

	print "\nTrimmomatic paired-end reads cleanning/trimming initiated ...\n\n";

	my $trimoPER1OUT = join("","$fastq_seed","_PE_R1parsed",".fq.gz");
	my $trimoSER1OUT = join("","$fastq_seed","_SE_R1parsed",".fq.gz");
	my $trimoPER2OUT = join("","$fastq_seed","_PE_R2parsed",".fq.gz");
	my $trimoSER2OUT = join("","$fastq_seed","_SE_R2parsed",".fq.gz");
	
	if ($adaptor ne '0') {
		system ( "java -jar $trimmomatic PE -phred$phred -threads $threads $outR1 $outR2 $trimoPER1OUT $trimoSER1OUT $trimoPER2OUT $trimoSER2OUT MINLEN:$minlen ILLUMINACLIP:$adaptor LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );
	} elsif ($adaptor eq '0') {
		system ( "java -jar $trimmomatic PE -phred$phred -threads $threads $outR1 $outR2 $trimoPER1OUT $trimoSER1OUT $trimoPER2OUT $trimoSER2OUT MINLEN:$minlen LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );
	}

	print "\nDONE.\n";
	
############################
# Trimming Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Trimming Single-End reads ...\n";

	print "\nConcatenating library $fastq_seed reads ... ";
	my $outR1 = join ("","$fastq_seed","_R1parsed",".merged",".fq",".gz");
	system ( "cat $fastq_seed*R1parsed.fq.gz > $outR1" );
	print "DONE.";

	print "\nTrimmomatic single-end reads cleanning/trimming initiated...\n\n";

	my $trimoSER1OUT = join("","$fastq_seed","_SE_R1parsed",".fq.gz");
	
	if ($adaptor ne '0') {
		system ( "java -jar $trimmomatic SE -phred$phred -threads $threads $outR1 $trimoSER1OUT MINLEN:$minlen ILLUMINACLIP:$adaptor LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );
	} elsif ($adaptor eq '0') {
		system ( "java -jar $trimmomatic SE -phred$phred -threads $threads $outR1 $trimoSER1OUT MINLEN:$minlen LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );
	}

	print "\nDONE.\n";
}
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
