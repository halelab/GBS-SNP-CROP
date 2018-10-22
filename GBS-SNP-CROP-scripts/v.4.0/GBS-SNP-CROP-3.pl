#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-3.pl
###########################################################################################################

use warnings;
use Getopt::Long qw(GetOptions);

######################
# The help function
######################
my $help = $ARGV[0];
my ($dataType,$barcodesID_file,$fastq_seed);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.0\n"
	."Step 3: Demultiplex\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-3.pl [options]\n\n"
	."Options:\n"
	."-d: Data type. Either PE (Paired-End) or SE (Single-End). String. Required.\n"
	."-b: BarcodeID file. File. Required.\n"
	."-fq: FASTQ file name seed. String. Required\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameters values
#################################
GetOptions(
'd=s' => \$dataType,            # string - "PE" or "SE"
'b=s' => \$barcodesID_file,     # file
'fq=s' => \$fastq_seed,         # string
) or die "$H";

#########################
# Starting GBS-SNP-CROP
#########################
print "\n#################################\n# GBS-SNP-CROP, Step 3, v.4.0\n#################################\n";
my $sttime = time;

# Creating directory
my $dir = "demultiplexed";
unless(-e $dir, or mkdir $dir) {die "Directory $dir cannot be created.\n";}

##################################
# Demultiplexing Paired-End data 
##################################

if ($dataType eq "PE") {
	print "Demultiplexing Paired-End reads ...\n";

	my $input1 = join ("","$fastq_seed","_PE_R1parsed",".fq.gz");
	print "\nCreating genotype-specific FASTQ files from $input1 file ...\n";

	my %barcode_hash;
	
	open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
	while(<$BAR>) {
		my $barcodesID = $_;
		chomp $barcodesID;
		my @barcode_ID = split("\t", $barcodesID);

		if ( $barcode_hash{$barcode_ID[0]} ) {
			die "Redundant barcodes in barcode-ID file!";
		} else {
			$barcode_hash{$barcode_ID[0]} = $barcode_ID[1];
		}
	}
	close $BAR;

	open my $IN1, '-|', 'gzip', '-dc', $input1 or die "Can't open FASTQ file: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R1","fq");
		open $key, " | gzip > ./demultiplexed/$filename.gz" or die "Can't open $_ file\n";
	}

	my @read1 = ();
	my $i1 = 1;
	while(<$IN1>) {
		if ($i1 % 4 != 0) {
			push @read1, $_;
			$i1++;
		} else {
			push @read1, $_;
			chomp (@read1);
			if ( $read1[0] =~ /^(@.*:N:0:)(\w{0,10})$/ && $barcode_hash{$2} ) {
				print $2 "$read1[0]\n$read1[1]\n$read1[2]\n$read1[3]\n";
				@read1 = ();
				$i1++;
				next;
			} else {
				@read1 = ();
				$i1++;
				next;
			}
		}
	}
	close $IN1;

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}", "fq");
		close $key;
	}
	print "DONE.\n";

	my $input2 = join ("","$fastq_seed","_PE_R2parsed",".fq.gz");
	print "\nCreating genotype-specific FASTQ files from $input2 file ...\n";

	open my $IN2, '-|', 'gzip', '-dc', $input2 or die "Can't open file $input2: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R2","fq");
		open $key, " | gzip > ./demultiplexed/$filename.gz" or die "Can't open $_ file\n";
	}

	my @read2 = ();
	my $i2 = 1;
	while(<$IN2>) {
		if ($i2 % 4 != 0) {
			push @read2, $_;
			$i2++;
		} else {
			push @read2, $_;
			chomp (@read2);
			if ( $read2[0] =~ /^(@.*:N:0:)(\w{0,10})$/ && $barcode_hash{$2} ) {
				print $2 "$read2[0]\n$read2[1]\n$read2[2]\n$read2[3]\n";
				@read2 = ();
				$i2++;
				next;
			} else {
				@read2 = ();
				$i2++;
				next;
			}
		}
	}
	close $IN2;

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}", "fq");
		close $key;
	}
	print "DONE.\n";
	
##################################
# Demultiplexing Single-End data 
##################################

} elsif ($dataType eq "SE") {
	print "Demultiplexing Single-End reads ...\n";
	
	my $input1 = join ("","$fastq_seed","_SE_R1parsed",".fq.gz");
	print "\nCreating genotype-specific FASTQ files from $input1 file ...\n";

	my %barcode_hash;
	
	open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
	while(<$BAR>) {
		my $barcodesID = $_;
		chomp $barcodesID;
		my @barcode_ID = split("\t", $barcodesID);

		if ( $barcode_hash{$barcode_ID[0]} ) {
			die "Redundant barcodes in barcode-ID file!";
		} else {
			$barcode_hash{$barcode_ID[0]} = $barcode_ID[1];
		}
	}
	close $BAR;

	open my $IN1, '-|', 'gzip', '-dc', $input1 or die "Can't open FASTQ file: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R1","fq");
		open $key, " | gzip > ./demultiplexed/$filename.gz" or die "Can't open $filename file\n";;
	}

	my @read1 = ();
	my $i1 = 1;
	while(<$IN1>) {
		if ($i1 % 4 != 0) {
			push @read1, $_;
			$i1++;
		} else {
			push @read1, $_;
			chomp (@read1);
			if ( $read1[0] =~ /^(@.*:N:0:)(\w{0,10})$/ && $barcode_hash{$2} ) {
				print $2 "$read1[0]\n$read1[1]\n$read1[2]\n$read1[3]\n";
				@read1 = ();
				$i1++;
				next;
			} else {
				@read1 = ();
				$i1++;
				next;
			}
		}
	}
	close $IN1;

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}", "fq");
		close $key;
	}
	print "DONE.\n";
}

print "\nPlease, see the 'demultiplexed' directory.\n";
print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
