#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-5.pl help
###########################################################################################################
##########################################################################################
# Requirement 1: BWA aligner (Li & Durbin, 2009)
# Requirement 2: SAMTools (Li et al., 2009)
##########################################################################################

use strict;
no warnings 'uninitialized';
use Getopt::Long qw(GetOptions);
use Parallel::ForkManager;

######################
# The help function
######################
my $help = $ARGV[0];
my ($bwa,$samtools,$dataType,$barcodesID_file,$Reference,$phred_Q,$map_q,$F,$f,$threads,$sam_add);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.0\n"
	."Step 5: Align with BWA-mem and process with SAMtools\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-5.pl [options]\n\n"
	."Options:\n"
	."-bw: Path to BWA executable file. String. Default: /usr/local/bin/bwa\n"
	."-st: Path to SAMTools executable file. String. Default: /usr/local/bin/samtools\n"
	."-d: Data type. Either PE (Paired-End) or SE (Single-End). String. Required.\n"
	."-b: BarcodeID file. File. Required.\n"
	."-ref: Reference FASTA file, either Mock Reference or true reference. File. Default: GSC.MR.Genome.fa\n"
	."-Q: Phred score base call quality. Numeric. Default: 30\n"
	."-q: Alignment quality. Numeric. Default: 30\n"
	."-F: SAMtools flags controlled by CAPS F. Numeric. Default: 2308\n"
	."-f: SAMtools flags controlled by small f. Numeric. Default: 0\n"
	."-t: Number of independent threads used. Numeric. Default: 10\n"
	."-Opt: If desired, any additional options for SAMtools view. String within “quotes”. Default: 0 (nothing)\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameters values
#################################
$bwa = '/usr/local/bin/bwa';		$samtools = '/usr/local/bin/samtools';
$Reference = 'GSC.MR.Genome.fa';	$phred_Q = 30;
$map_q = 30;				$F = 2308;
$f = 2;					$threads = 10;	
$sam_add = 0;

GetOptions(
'bw=s' => \$bwa,	          	# string
'st=s' => \$samtools,          	# string
'd=s' => \$dataType,          	# string
'b=s' => \$barcodesID_file,     # file
'ref=s' => \$Reference,         # file
'Q=s' => \$phred_Q,             # numeric
'q=s' => \$map_q,               # numeric
'F=s' => \$F,                	# numeric
'f=s' => \$f,               	# numeric 
't=s' => \$threads,             # numeric
'Opt=s' => \$sam_add,           # string
) or die "$H\n";

#########################
# Starting GBS-SNP-CROP
#########################
print "\n#################################\n# GBS-SNP-CROP, Step 5, v.4.0\n#################################\n";
my $pm = new Parallel::ForkManager($threads);
my $sttime = time;

# Creating a directory
my $dir = "alignments"; 
unless(-e $dir, or mkdir $dir) {die "Directory $dir cannot be created.\n";}

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

	if ($F > 0 && $f > 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q $map_q -f $f -F $F $sam_add $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q $map_q -F $F $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add ne '0') ) {
		system ( "samtools view -b -q $map_q -f $f $sam_add $input_sam > $view_out" );
	} elsif ($f > 0 && $F > 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q $map_q -f $f -F $F $input_sam > $view_out" );
	} elsif ($F > 0 && $f == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q $map_q -F $F $input_sam > $view_out" );
	} elsif ($f > 0 && $F == 0 && ($sam_add eq '0') ) {
		system ( "samtools view -b -q $map_q -f $f $input_sam > $view_out" );
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
	system ("samtools mpileup -Q $phred_Q -q $map_q -B -C 50 -f $Reference $input > $mpileup");
	$pm->finish;
}
$pm->wait_all_children;
print "DONE.\n\n";

system ( "mv *bam* *.sam ./alignments" );
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
