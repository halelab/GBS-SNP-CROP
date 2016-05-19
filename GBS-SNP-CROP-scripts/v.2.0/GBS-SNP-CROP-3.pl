#!/usr/bin/perl

#############################################################################################################################
# GBS-SNP-CROP, Step 3. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
#############################################################################################################################

use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-3.pl -d <data type, PE = Paired-End or SE = Single-End> -b <barcode-ID file name> -fq <FASTQ file name seed>\n";
my $Manual = "Please see the User Manual on the GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n";

my ($dataType,$barcodesID_file,$fastq_seed);

GetOptions(
'd=s' => \$dataType,            # string - "PE" or "SE"
'b=s' => \$barcodesID_file,     # file
'fq=s' => \$fastq_seed,         # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 3, v.2.0\n#################################\n";

sub main {
	my $dir = "demultiplexed";
	unless(-e $dir, or mkdir $dir) {die "Directory $dir does not exist and cannot be created.\n";}
}
main();
##################################
# Demultiplexing Paired-End data 
##################################

if ($dataType eq "PE") {
	print "Demultiplexing Paired-End reads...\n";

	my $input1 = join ("","$fastq_seed","_PE_R1parsed",".fastq");
	print "\nCreating genotype-specific FASTQ files from $input1 file...\n";

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

	open my $IN1, "<", "$input1" or die "Can't open FASTQ file: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R1","fastq");
		open $key, ">", "./demultiplexed/$filename" or die "Can't open $filename!";
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
		my $filename = join(".", "$barcode_hash{$key}", "fastq");
		close $key;
	}
	print "DONE!\n";

	my $input2 = join ("","$fastq_seed","_PE_R2parsed",".fastq");
	print "\nCreating genotype-specific FASTQ files from $input2 file...\n";

	open my $IN2, "<", "$input2" or die "Can't open FASTQ file: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R2","fastq");
		open $key, ">", "./demultiplexed/$filename" or die "Can't open $filename!";
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
		my $filename = join(".", "$barcode_hash{$key}", "fastq");
		close $key;
	}
	print "DONE!\n";

##################################
# Demultiplexing Single-End data 
##################################

} elsif ($dataType eq "SE") {
	print "Demultiplexing Single-End reads...\n";

	my $input1 = join ("","$fastq_seed","_SE_R1parsed",".fastq");
	print "\nCreating genotype-specific FASTQ files from $input1 file...\n";

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

	open my $IN1, "<", "$input1" or die "Can't open FASTQ file: $!\n";

	foreach my $key (keys %barcode_hash) {
		my $filename = join(".", "$barcode_hash{$key}","R1","fastq");
		open $key, ">", "./demultiplexed/$filename" or die "Can't open $filename!";
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
		my $filename = join(".", "$barcode_hash{$key}", "fastq");
		close $key;
	}
	print "DONE!\n";
}

print "Please, see the 'demultiplexed' directory.\n";

print "\n\nPlease cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
