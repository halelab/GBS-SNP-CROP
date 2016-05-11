#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 5. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
##########################################################################################################################

##########################################################################################
# Requirement 1: BWA aligner (Li & Durbin, 2009)
# Requirement 2: SAMTools (Li et al., 2009)
##########################################################################################

use strict;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-5.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file name>  -ref <reference FASTA file> -Q <Phred score> -q <mapping quality score>\n"
." -f <SAMtools -f flag> -F <SAMtools _F flag> -t <threads> -Opt <any additional desired SAMtools options>.\n";
my $Manual = "Please see the User Manual on the GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($dataType,$barcodesID_file,$Reference,$phred_Q,$map_q,$f,$F,$threads,$sam_add);

GetOptions(
'd=s' => \$dataType,          # string - "PE" or "SE"
'b=s' => \$barcodesID_file,     # file
'ref=s' => \$Reference,         # file
'Q=s' => \$phred_Q,             # numeric
'q=s' => \$map_q,               # numeric
'f=s' => \$f,                # numeric 
'F=s' => \$F,                # numeric
't=s' => \$threads,             # numeric
'Opt=s' => \$sam_add,           # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 5, v.2.0\n#################################\n";

my @files = ();

open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
while(<$BAR>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	push @files, $TaxaNames;
}
close $BAR;
chomp (@files);

# 1. BWA procedures

# 1.1 Index
print "\nIndexing reference FASTA file...\n";
system ( "bwa index -a bwtsw $Reference" );
print "DONE.\n";
	
############################
# Aligning Paired-End data 
############################

if ($dataType eq "PE") {

# 1.2 BWA-mem mapping
	foreach my $file (@files) {
		my $input_R1 = join (".", "$file","R1","fastq");
		my $input_R2 = join (".", "$file","R2","fastq");
		my $BWA_out = join(".","$file","sam");
		print "\nMapping paired $input_R1 $input_R2 files to $Reference...\n";
		system ( "bwa mem -t $threads -M $Reference $input_R1 $input_R2 > $BWA_out" );
	}
	print "\n\nBWA-mem mapping was completed!\n\n";

############################
# Aligning Single-End data 
############################

} elsif ($dataType eq "SE") {

# 1.2 BWA-mem mapping
	foreach my $file (@files) {
		my $input_R1 = join (".", "$file","R1","fastq");
		my $BWA_out = join(".","$file","sam");
		print "\nMapping single $input_R1 file to $Reference...\n";
		system ( "bwa mem -t $threads -M $Reference $input_R1 > $BWA_out" );
	}
	print "\n\nBWA-mem mapping was completed!\n\n";
}

# 2. SAMTools procedures

# 2.1 SAM to BAM
foreach my $file (@files) {
	my $input_sam = join (".", "$file","sam");
	my $view_out = join(".","$file","bam");
	print "\nProcessing $input_sam file...";

	if ($F > 0 && $f > 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q$phred_Q -f$f -F$F $sam_add $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q$phred_Q -F$F $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q$phred_Q -f$f $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F > 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q$phred_Q -f$f -F$F $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q$phred_Q -F$F $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q$phred_Q -f$f $input_sam > $view_out" );
	} else {
		print "Unable to proceeed; please re-check the syntax of all declared SAMtools flags and options...";
	}
}
print "\nAll SAM files were converted into binary (BAM) files.\n";

# 2.2 Sorting BAM files
foreach my $file (@files) {
	my $input_bam = join (".", "$file","bam");
	my $sort_out = join(".","$file","sorted");
	print "\nSorting $input_bam file...";
	system ( "samtools sort $input_bam $sort_out" );
}
print "\nAll BAM files were sorted.\n";

# 2.3 Index sorted BAM files
foreach my $file (@files) {
	my $input_sorted = join (".","$file","sorted","bam");
	print "\nIndexing $input_sorted file...";
	system ( "samtools index $input_sorted" );
}
print "\nAll sorted BAM files were indexed.\n";

# Index reference FASTA file
print "\nIndexing the reference genome FASTA file...";
system ( "samtools faidx $Reference" );
print "DONE.\n";

# Mpileup SNPs discovery
foreach my $file (@files) {
	my $input = join (".", "$file","sorted","bam");
	my $mpileup = join (".", "$file","mpileup");
	print "\nProducing mpileup file from $file ...\n";
	system ("samtools mpileup -Q$phred_Q -q$map_q -B -C 50 -f $Reference $input > $mpileup");
}
print "\n\nAn mpileup file was successfully created for each genotype.\n";

sub main {
	my $dir = "alignments";
	unless(-e $dir, or mkdir $dir) {die "Directory $dir does not exist and cannot be created.\n";}
}
main();

system ( "mv *bam* ./alignments" );
system ( "rm *.sam" );

print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
