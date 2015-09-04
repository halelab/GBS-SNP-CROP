#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 3. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-3.pl -b <barcode-ID file name> -fqN <FASTQ file name> 
-r <the type of reads being demultiplexed, R1 or R2>\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($barcodesID_file,$fastq_file,$read_pair);

GetOptions(
'b=s' => \$barcodesID_file,     # file
'fqN=s' => \$fastq_file,     	# string
'r=s' => \$read_pair, 		  	# string 
) or die "$Usage\n$Manual\n";

my %barcode_hash;
	
open IN, "$barcodesID_file" or die "Can't find barcode_ID file\n";
	
while(<IN>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode_ID = split("\t", $barcodesID);

	if ( $barcode_hash{$barcode_ID[0]} ) {
		die "Redundant barcodes in barcode-ID file!";
	} else {
		$barcode_hash{$barcode_ID[0]} = $barcode_ID[1];
	}
}
	
close IN;

sub main {
	my $dir = "demultiplexed";
   	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
 }
main();

open IN, "$fastq_file" or die "Can't open FASTQ file: $!\n";

foreach my $key (keys %barcode_hash) {
	my $filename = join(".", "$barcode_hash{$key}","$read_pair","fastq");
	open $key, ">", "./demultiplexed/$filename";
}

my @read = ();
my $i = 1;
while(<IN>) {
	if ($i % 4 != 0) {
		push @read, $_;
		$i++;
	} else {
		push @read, $_;
		chomp (@read);
		if ( $read[0] =~ /^(@.*:N:0:)(\w{0,10})$/ && $barcode_hash{$2} ) {
			print $2 "$read[0]\n$read[1]\n$read[2]\n$read[3]\n";
			@read = ();
			$i++;
			next;
		} else {
			@read = ();
			$i++;
			next;
		}
	}
}

foreach my $key (keys %barcode_hash) {
	my $filename = join(".", "$barcode_hash{$key}", "fastq");
	close $key;
}
print "\nGenotypes specific FASTQ files were successfully extracted from $fastq_file.\nPlease, see the 'demultiplexed' directory.\n";

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;