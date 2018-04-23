#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 5. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y
##########################################################################################################################

##########################################################################################
# Requirement 1: BWA aligner (Li & Durbin, 2009)
# Requirement 2: SAMTools (Li et al., 2009)
##########################################################################################

use strict;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);
use Parallel::ForkManager;

my $Usage = "Usage: perl GBS-SNP-CROP-5.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file name>  -ref <reference FASTA file> -Q <Phred score> -q <mapping quality score>\n"
." -f <SAMTools -f flag> -F <SAMTools _F flag> -t <threads> -Opt <any additional desired SAMTools options>.\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($dataType,$barcodesID_file,$Reference,$phred_Q,$map_q,$f,$F,$threads,$sam_add);

GetOptions(
'd=s' => \$dataType,          	# string - "PE" or "SE"
'b=s' => \$barcodesID_file,     # file
'ref=s' => \$Reference,         # file
'Q=s' => \$phred_Q,             # numeric
'q=s' => \$map_q,               # numeric
'f=s' => \$f,               	# numeric 
'F=s' => \$F,                	# numeric
't=s' => \$threads,             # numeric
'Opt=s' => \$sam_add,           # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 5, v.3.0\n#################################\n";
my $pm = new Parallel::ForkManager($threads);
my $sttime = time;

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
print "\nIndexing reference FASTA file ...\n";
system ( "bwa index -a bwtsw $Reference" );
print "DONE.\n\n";
	
############################
# Aligning Paired-End data 
############################

if ($dataType eq "PE") {

# 1.2 BWA-mem mapping
	foreach my $file (@files) {
		my $input_R1 = join (".", "$file","R1","fq","gz");
  	    my $input_R2 = join (".", "$file","R2","fq","gz");
        my $BWA_out = join(".","$file","sam");
		print "Mapping paired $input_R1 $input_R2 files to $Reference ...\n";
		system ( "bwa mem -t $threads -M $Reference $input_R1 $input_R2 > $BWA_out" );
	}
	print "DONE.\n";

############################
# Aligning Single-End data 
############################

} elsif ($dataType eq "SE") {

# 1.2 BWA-mem mapping
	foreach my $file (@files) {
		my $input_R1 = join (".", "$file","R1","fq","gz");
  	    my $BWA_out = join(".","$file","sam");
		print "Mapping single $input_R1 file to $Reference ...\n";
		system ( "bwa mem -t $threads -M $Reference $input_R1 > $BWA_out" );
	}
	print "DONE.\n";
}

# 2. SAMTools procedures

# 2.1 SAM to BAM
print "\nProcessing the SAM files ...";
foreach my $file (@files) {
	my $pid = $pm->start and next;
	my $input_sam = join (".", "$file","sam");
    my $view_out = join(".","$file","bam");
#	print "\nProcessing $input_sam file ...";

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
		print "Unable to proceeed; please re-check the syntax of all declared SAMTools flags and options...";
	}
	$pm->finish;
}
$pm->wait_all_children;
print "\nDONE.\n";

# 2.2 Sorting BAM files
print "\nSorting the BAM files ...";
foreach my $file (@files) {
	my $pid = $pm->start and next;
	my $input_bam = join (".", "$file","bam");
    my $sort_out = join(".","$file","sorted.bam");
#	print "\nSorting $input_bam file ...";
	system ( "samtools sort $input_bam -o $sort_out" );
	$pm->finish;
}
$pm->wait_all_children;
print "\nDONE.\n";

# 2.3 Index sorted BAM files
print "\nIndexing the sorted BAM files ...";
foreach my $file (@files) {
	my $pid = $pm->start and next;
	my $input_sorted = join (".","$file","sorted","bam");
#	print "\nIndexing $input_sorted file ...";
	system ( "samtools index $input_sorted" );
	$pm->finish;
}
$pm->wait_all_children;
print "\nDONE.\n";

# Index reference FASTA file
print "\nIndexing the reference genome FASTA file ...";
system ( "samtools faidx $Reference" );
print "\nDONE.\n\n";

# Mpileup SNPs discovery
print "Producing the mpileup files ...\n";
foreach my $file (@files) {
	my $pid = $pm->start and next;
	my $input = join (".", "$file","sorted","bam");
	my $mpileup = join (".", "$file","mpileup");
#	print "Producing mpileup file from $file ...\n";
	system ("samtools mpileup -Q$phred_Q -q$map_q -B -C 50 -f $Reference $input > $mpileup");
	$pm->finish;
}
$pm->wait_all_children;
print "DONE.\n\n";

sub main {
   	my $dir = "alignments";
  	unless(-e $dir, or mkdir $dir) {die "Directory $dir does not exist and cannot be created.\n";}
}   
main();

system ( "mv *bam* ./alignments" );
system ( "rm *.sam" );
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
