#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 5. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

##########################################################################################
# Requirement 1: BWA aligner (Li & Durbin, 2009)
# Requirement 2: SAMTools (Li et al., 2009)
##########################################################################################

use strict;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-5.pl -b <barcode-ID file name>  -ref <reference FASTA file> 
-sam_f <SAMTools flags controlled by small f> -sam_F <SAMTools flags controlled by CAPS F> -t <threads>.\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($barcodesID_file,$Reference,$f,$F,$threads);

GetOptions(
'b=s' => \$barcodesID_file,     # file
'ref=s' => \$Reference,   		# file
'samf=s' => \$f,     	  		# numeric 
'samF=s' => \$F,          		# numeric
't=s' => \$threads,       		# numeric
) or die "$Usage\n$Manual\n";

my @files = ();

open IN, "$barcodesID_file" or die "Can't find barcode_ID file\n";
	
while(<IN>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	
	push @files, $TaxaNames;
}

chomp (@files);

##################
# 1. BWA procedures
##################

# Index
print "\n\nIndexing reference fasta file...\n";
system ( "bwa index -a bwtsw $Reference" );
print "DONE.\n";


# BWA-mem mapping
foreach my $file (@files) {
		my $input_R1 = join (".", "$file","R1","fastq");
        my $input_R2 = join (".", "$file","R2","fastq");
        my $BWA_out = join(".","$file","sam");
		print "\nMapping paired $input_R1 $input_R2 FASTQ files to $Reference...\n";
		system ( "bwa mem -t $threads -M $Reference $input_R1 $input_R2 > $BWA_out" );
	
}
print "\n\nBWA-mem mapping was done!\n\n";


########################
# 2. SAMTools procedures
########################

# SAM to BAM
foreach my $file (@files) {
        my $input_sam = join (".", "$file","sam");
        my $view_out = join(".","$file","bam");
		print "\nProcessing $input_sam file...";
		if ($F > 0 && $f > 0) {
		system ( "samtools view -b -q30 -F$F -f$f $input_sam > $view_out" );
		} elsif ($F > 0 && $f == 0) {
		system ( "samtools view -b -q30 -F$F $input_sam > $view_out" );
		} elsif ($f > 0 && $F == 0) {
		system ( "samtools view -b -q30 -f$f $input_sam > $view_out" );
		} else {
		print "Unable to proceeed; please re-check your SAMTools flags values...";
		}	
}

print "\nAll SAM files were converted into binary files.\n";


# Sorting BAM files
foreach my $file (@files) {
        my $input_bam = join (".", "$file","bam");
        my $sort_out = join(".","$file","sorted");
		print "\nSorting $input_bam file...";
		system ( "samtools sort $input_bam $sort_out" );
}
print "\nAll BAM files were sorted.\n";


# Index sorted BAM files
foreach my $file (@files) {
        my $input_sorted = join (".","$file","sorted","bam");
		print "\nIndexing $input_sorted file...";
		system ( "samtools index $input_sorted" );
}
print "\nAll sorted BAM files were indexed.\n";


# Index reference FASTA file
print "\nIndexing the reference genome fasta file...";
system ( "samtools faidx $Reference" );
print "DONE.\n";


# Mpileup SNPs discovery
foreach my $file (@files) {
	my $input = join (".", "$file","sorted","bam");
	my $mpileup = join (".", "$file","mpileup");
	print "\nProducing mpileup file from $file ...\n";
	system ("samtools mpileup -Q30 -B -C 50 -f $Reference $input > $mpileup");
}
print "\n\nMpileup files were successfully created for each genotype.\n";

sub main {
   	my $dir = "alignments";
  	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
}   
main();

system ( "mv *bam* ./alignments" );
system ( "rm *.sam" );

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;