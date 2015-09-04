#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 2. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

##########################################################################################
# Requirement 1: Java 7 or high
# Requirement 2: Trimmomatic software (Bolger et al., 2014)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-2.pl -fq <FASTQ file name seed> -t <threads> -ph <Phred quality score: 33 or 64> 
-l <Trimmomatic LEADING value> -sl <Trimmomatic SLIDINGWINDOW value> -tr <Trimmomatic TRAILING value> -m <Trimmomatic MINLEN value>\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($fastq_seed,$threads,$phred,$leading,$sliding,$trailing,$minlen);

GetOptions(
'fq=s' => \$fastq_seed,       # string
't=s' => \$threads,			  # numeric 
'ph=s' => \$phred,			  # numeric
'l=s' => \$leading,			  # numeric
'sl=s' => \$sliding,			  # numeric
'tr=s' => \$trailing,		  # numeric
'm=s' => \$minlen,			  # numeric 
) or die "$Usage\n$Manual\n";

print "\n\nConcatenating library $fastq_seed R1 reads...";
my $outR1 = join ("","$fastq_seed","_R1parsed",".fastq");
system ( "cat $fastq_seed*R1parsed.fastq > $outR1" );
print "DONE.";

print "\n\nConcatenating library $fastq_seed R2 reads...";
my $outR2 = join ("","$fastq_seed","_R2parsed",".fastq");
system ( "cat $fastq_seed*R2parsed.fastq > $outR2" );
print "DONE.";

print "\n\nTrimmomatic paire-end reads cleanning start...\n";

my $trimoPER1OUT = join("","$fastq_seed","_PE_R1parsed",".fastq");
my $trimoSRR1OUT = join("","$fastq_seed","_SR_R1parsed",".fastq");
my $trimoPER2OUT = join("","$fastq_seed","_PE_R2parsed",".fastq");
my $trimoSRR2OUT = join("","$fastq_seed","_SR_R2parsed",".fastq");

system ( "java -jar /usr/local/bin/trimmomatic-0.33.jar PE -threads $threads -phred$phred $outR1 $outR2 $trimoPER1OUT $trimoSRR1OUT $trimoPER2OUT $trimoSRR2OUT LEADING:$leading SLIDINGWINDOW:$sliding TRAILING:$trailing MINLEN:$minlen" );

print "\nRead trimming completed.\n";

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;