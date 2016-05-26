#!/usr/bin/perl

##########################################################################################################################
# GBS-SNP-CROP, Step 1. For description, please see Melo et al. (2016) BMC Bioinformatics 17:29. DOI 10.1186/s12859-016-0879-y.
##########################################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Zlib;

my $Usage = "Usage: perl GBS-SNP-CROP-1.pl -b <barcode-ID file name> -fq <FASTQ file name seed>\n" 
."-s <start number of FASTQ file> -e <end number of FASTQ file> -enz1 <enzyme 1 restriction site residue sequence>\n" 
."-enz2 <enzyme 2 restriction site residue sequence>.\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics 17:29. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($barcodesID_file,$fastq_seed,$fastq_start_num,$fastq_end_num,$enzyme1,$enzyme2);

GetOptions(
'b=s' => \$barcodesID_file,   # file
'fq=s' => \$fastq_seed,       # string
's=s' => \$fastq_start_num,   # numeric 
'e=s' => \$fastq_end_num,     # numeric
'enz1=s' => \$enzyme1,        # string
'enz2=s' => \$enzyme2,        # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 1, v.1.0\n#################################\n";

my $RC_enz1 = reverse $enzyme1;
$RC_enz1 =~ tr/acgtACGT/tgcaTGCA/;
my $RC_enz2 = reverse $enzyme2;
$RC_enz2 =~ tr/acgtACGT/tgcaTGCA/;

sub main {
   	my $dir1 = "summaries"; 
   	my $dir2 = "distribs"; 
   	my $dir3 = "parsed";
   	my $dir4 = "singles"; 
   	unless(-e $dir1, or mkdir $dir1) {die "Directory $dir1 just exist.\n";}
   	unless(-e $dir2, or mkdir $dir2) {die "Directory $dir2 just exist.\n";}
   	unless(-e $dir3, or mkdir $dir3) {die "Directory $dir3 just exist.\n";}
	unless(-e $dir4, or mkdir $dir4) {die "Directory $dir4 just exist.\n";}
}

main();

for ( my $file_index = $fastq_start_num; $file_index <= $fastq_end_num; $file_index++ ) {

	my $R1file = "";
	my $R2file = "";
	my $fileseed = "";
	
	if ( length ($file_index) == 1 ) {
		$R1file = join ("", "$fastq_seed","_R1_00", "$file_index", ".fastq.gz");
		$R2file = join ("", "$fastq_seed","_R2_00", "$file_index", ".fastq.gz");
		$fileseed = join ("", "$fastq_seed","_00", "$file_index");
	} elsif ( length ($file_index) == 2 ) {
		$R1file = join ("", "$fastq_seed","_R1_0", "$file_index", ".fastq.gz");
		$R2file = join ("", "$fastq_seed","_R2_0", "$file_index", ".fastq.gz");
		$fileseed = join ("", "$fastq_seed","_0", "$file_index");
	} elsif ( length ($file_index) == 3 ) {
		$R1file = join ("", "$fastq_seed","_R1_", "$file_index", ".fastq.gz");
		$R2file = join ("", "$fastq_seed","_R2_", "$file_index", ".fastq.gz");
		$fileseed = join ("", "$fastq_seed","_", "$file_index");
	} else {
	die "Unable to identify input fastq.gz files.\n";
	}

	print "\nProcessing files $R1file and $R2file\n";
	print "Loading data...";
	
	my @R1read;
	my @R1reads;
	
	open my $IN1, '-|', 'gzip', '-dc', $R1file or die "Can't open file $R1file: $!\n";
	my $i = 1;
	while(<$IN1>) {
		if ($i % 4 != 0) {
			push @R1read, $_;
			$i++;
		} else {
			push @R1read, $_;
			chomp (@R1read);
			push @R1reads, [ @R1read ];
			@R1read = ();
			$i++;
		}
	}
	close $IN1;
		
	my @R2read;
	my @R2reads;
	
	open my $IN2, '-|', 'gzip', '-dc', $R2file or die "Can't open file $R2file: $!\n";
	my $j = 1;
	while(<$IN2>) {
		if ($j % 4 != 0) {
			push @R2read, $_;
			$j++;
		} else {
			push @R2read, $_;
			chomp (@R2read);
			push @R2reads, [ @R2read ];
			@R2read = ();
			$j++;
		}
	}
	close $IN2;
		
	print "\nConsolidating R1 and R2 reads...";
		
	my $size = scalar @R1reads - 1;
	for (my $k = 0; $k <= $size; $k++) {
		
		if ( $R1reads[$k][0] =~ /^(@.*:N:0:)\w{0,50}$/ ) {
			$R1reads[$k][0] = $1;
		} else {
			$R1reads[$k][0] = $R1reads[$k][0];
		}
		
		$R2reads[$k][0] = "";
		$R1reads[$k][2] = $R1reads[$k][3];
		if ( $R2reads[$k][3] =~ /^#+$/ ) {
			$R1reads[$k][3] = "";
			$R2reads[$k][1] = "";
			$R1reads[$k][4] = "";
			$R2reads[$k][3] = "";
		} else {
		$R1reads[$k][3] = $R2reads[$k][1];
		$R2reads[$k][1] = "";
		$R1reads[$k][4] = $R2reads[$k][3];
		$R2reads[$k][3] = "";
		}
	}
		
	@R2reads = ();
	
	print "\nExtracting barcodes from R1 reads (exact matches or one mismatch)...";
		
	my @barcodes = ();
	my $barcode_string = "";
	
	open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
	while ( <$BAR> ) {
		my $barcodeID = $BAR;
		chomp $barcodeID;
		my @barcode = split("\t", $barcodeID);
		my $barcode_list = $barcode[0];
		my $TaxaNames = $barcode[1];
		push @barcodes, $barcode_list;
		$barcode_string = "$barcode_string $_";
	}
	close $BAR;
	chomp (@barcodes);
		
	$barcode_string = "$barcode_string ";
	$barcode_string =~ s/\n/ /g;
	
	my $total_raw_reads = $size + 1;
	my $no_R1_RE_site_tally = 0;
	my $no_R1_barcode_tally = 0;
	
	for (my $k = 0; $k <= $size; $k++) {
		my $R1seq = $R1reads[$k][1];
		if ( $R1seq =~ /^(\w{5,10})$enzyme1(\w+)$/ ) {
			my $inline_index = $1;
			$R1reads[$k][1] = $2;
			my $trim_length = length ($1) + 4;
			$R1reads[$k][2] = substr ( $R1reads[$k][2], $trim_length, $trim_length + length ($2) );
			
			if ( $barcode_string =~ /\s$inline_index\s/ ) {
				$R1reads[$k][0] = "$R1reads[$k][0]$inline_index";
				$R1reads[$k][5] = "$inline_index";
				next;
			} else {
				$R1reads[$k][5] = "";
			}
			
			my $barcode_match = $inline_index;
			my $hd_status = 0;
			
			foreach my $barcode (@barcodes) {
				if ( ( length ($inline_index) == length ($barcode) ) && hd ( $inline_index, $barcode ) == 1 ) {
					$barcode_match = $barcode;
					$hd_status++;
					next;
				} else {
					next;
				}
			}
						
			if ( $hd_status == 1 ) {
				$R1reads[$k][0] = "$R1reads[$k][0]$barcode_match";
				$R1reads[$k][5] = "$barcode_match";
			    next;
			} else {
				$no_R1_barcode_tally++;
				$R1reads[$k][0] = "";
				$R1reads[$k][1] = "";
				$R1reads[$k][2] = "";
				$R1reads[$k][3] = "";
				$R1reads[$k][4] = "";
				$R1reads[$k][5] = "";
				next;
			}
			
		} else {
			$no_R1_RE_site_tally++;
			$R1reads[$k][0] = "";
			$R1reads[$k][1] = "";
			$R1reads[$k][2] = "";
			$R1reads[$k][3] = "";
			$R1reads[$k][4] = "";
			$R1reads[$k][5] = "";
			next;
		}
	}
	
	print "\nTrimming R1 reads based on presence of $enzyme2 site, barcodes, and Illumina tag...";
	
	for (my $k = 0; $k <= $size; $k++) {
		if ( $R1reads[$k][0] eq "" ) {
			next;
		} else {
			my $R1seq = $R1reads[$k][1];
			if ( $R1seq =~ /^(\w+)$RC_enz2[AN][GN][AN][TN][CN][GN][GN][AN][AN]\w+$/ ) {
				my $R1_read = $1;
				my $R1_read_length = length ( $R1_read );
				$R1reads[$k][1] = $R1_read;
				$R1reads[$k][2] = substr ( $R1reads[$k][2], 0, $R1_read_length );
			} else {
				next;
			}
		}
	}
	
	print "\nTrimming R2 reads based on presence of $enzyme1 site, barcodes, and Illumina tag...";
	
	for (my $k = 0; $k <= $size; $k++) {
		if ( $R1reads[$k][3] eq "" ) {
			next;
		} else {
			my $R2seq = $R1reads[$k][3];
			if ( $R2seq =~ /^(\w+)$RC_enz1(\w{5,10})[AN][GN][AN][TN]/ ) {
				my $R2_read = $1;
				my $inline_index = scalar reverse $2;
				$inline_index =~ tr/actgACTG/tgacTGAC/;
				my $index_hd = hd ( $inline_index, $R1reads[$k][5] ); 
				
				if ( $index_hd <= 1 ) {
					$R1reads[$k][3] = $R2_read;
					my $R2_read_length = length ($R1reads[$k][3]);
					$R1reads[$k][4] = substr ( $R1reads[$k][4], 0, $R2_read_length );
					next;
				} else {
					next;
				}
			} else {
				next;
			}
		}
	}
	
	print "\nTrimming/culling R2 reads based on quality/length...";
	
	my $R2_min_length_tally = 0;
	my $dN_R2_tally = 0;
	
	for (my $k = 0; $k <= $size; $k++) {
		if ( ( $R1reads[$k][3] eq "" ) && ( $R1reads[$k][4] eq "" ) ) {
			next;
		} else {
			if ( $R1reads[$k][4] =~ /(#{4,})$/ ) {
				my $bad_bases = $1;
				my $trimmed_read_length = length ($R1reads[$k][4]) - length ($bad_bases);
				$R1reads[$k][3] = substr ( $R1reads[$k][3], 0, $trimmed_read_length );
				$R1reads[$k][4] = substr ( $R1reads[$k][4], 0, $trimmed_read_length );
			}
	
			my @R2_seq = split ("", $R1reads[$k][3]);
			my $N_count = 0;
			foreach my $base (@R2_seq) {
				if ( $base eq "N" ) {
					$N_count++;
				}
			}
			
			if ( $N_count >= (length ($R1reads[$k][3]) / 2 )) {
				$dN_R2_tally++;
				$R1reads[$k][3] = "";
				$R1reads[$k][4] = "";
				next;
			} else {
				next;
			}
		}
	}
	
	print "\nWriting results to output files...";
	
	my $R1singles_file = join (".","$fileseed","R1singles.fastq");
	my $R1out_file = join (".","$fileseed","R1parsed.fastq");
	my $R2out_file = join (".","$fileseed","R2parsed.fastq");
	my $summary_file = join (".","$fileseed","summary");
	my $read_dist_file = join (".","$fileseed","dist");
	
	my $R1_singleton_tally = 0;
	my $usable_pairs_tally = 0;

	open my $R1_singles, " | gzip > ./singles/$R1singles_file.gz" or die "Can't open $R1singles_file\n";
	open my $R1_OUT, " | gzip > ./parsed/$R1out_file.gz" or die "Can't open $R1out_file\n";
	open my $R2_OUT, " | gzip > ./parsed/$R2out_file.gz" or die "Can't open $R2out_file\n";

	my %R1_singles_lengths;
	my %R1_lengths;
	my %R2_lengths;
	my %R1R2_diff_lengths;
	
	for (my $k = 0; $k <= $size; $k++) {
		if ( $R1reads[$k][5] eq "" ) {
			next;
		} elsif ( ( $R1reads[$k][1] ne "" ) && ( $R1reads[$k][3] eq "" ) ) {
			$R1_singleton_tally++;
			print $R1_singles "$R1reads[$k][0]\n";
			print $R1_singles "$R1reads[$k][1]\n";
			print $R1_singles "+\n";
			print $R1_singles "$R1reads[$k][2]\n";
			
			my $R1_len = length ($R1reads[$k][1]);
			if ( exists $R1_singles_lengths{$R1_len} ) {
				$R1_singles_lengths{$R1_len}++;
			} else {
				$R1_singles_lengths{$R1_len} = 1;
			}
			next;
			
		} else {
			$usable_pairs_tally++;
			print $R1_OUT "$R1reads[$k][0]\n";
			print $R1_OUT "$R1reads[$k][1]\n";
			print $R1_OUT "+\n";
			print $R1_OUT "$R1reads[$k][2]\n";
	
			my $R1_len = length ($R1reads[$k][1]);
			if ( exists $R1_lengths{$R1_len} ) {
				$R1_lengths{$R1_len}++;
			} else {
				$R1_lengths{$R1_len} = 1;
			}
			
			my $R2_header = $R1reads[$k][0];
			$R2_header =~ s/ 1:/ 2:/;
			print $R2_OUT "$R2_header\n";
			print $R2_OUT "$R1reads[$k][3]\n";
			print $R2_OUT "+\n";
			print $R2_OUT "$R1reads[$k][4]\n";
	
			my $R2_len = length ($R1reads[$k][3]);
			if ( exists $R2_lengths{$R2_len} ) {
				$R2_lengths{$R2_len}++;
			} else {
				$R2_lengths{$R2_len} = 1;
			}
	
			my $diff = $R1_len - $R2_len;
			if ( exists $R1R2_diff_lengths{$diff} ) {
				$R1R2_diff_lengths{$diff}++;
			} else {
				$R1R2_diff_lengths{$diff} = 1;
			}
			next;
		}
	}
	close $R1_singles;
	close $R1_OUT;
	close $R2_OUT;

	open my $summary_OUT, ">", "./summaries/$summary_file" or die "Can't open $summary_file\n";

	print $summary_OUT "SUMMARY:\n";
	print $summary_OUT "Total raw read pairs = $total_raw_reads\n";
	print $summary_OUT "R1 reads with no identifiable restriction site = $no_R1_RE_site_tally\n";
	print $summary_OUT "R1 reads with no identifiable barcode = $no_R1_barcode_tally\n";
	print $summary_OUT "Usable R1 singletons = $R1_singleton_tally\n\n";
	my $usable_pairs_percent = ($usable_pairs_tally / $total_raw_reads) * 100;
	print $summary_OUT "Percentage of usable paired reads = $usable_pairs_percent\n";
	my $R1_singleton_percent = ( $R1_singleton_tally / $total_raw_reads ) * 100;
	print $summary_OUT "Percentage of usable, unpaired R1 reads = $R1_singleton_percent\n";
	close $summary_OUT;
	
	open my $dist_OUT, '>', "./distribs/$read_dist_file" or die "Can't open $read_dist_file\n";

	print $dist_OUT "Distribution of pairable R1 read lengths:\n";
	foreach my $key ( sort ( keys %R1_lengths ) ) {
		print $dist_OUT "R1", "\t", "$key", "\t", "$R1_lengths{$key}", "\n";
	}
	print $dist_OUT "\n\n";
	
	print $dist_OUT "Distribution of pairable R2 read lengths:\n";
	foreach my $key ( sort ( keys %R2_lengths ) ) {
		print $dist_OUT "R2", "\t", "$key", "\t", "$R2_lengths{$key}", "\n";
	}
	print $dist_OUT "\n\n";
	
	print $dist_OUT "Distribution of paired R1-R2 read length differences:\n";
	foreach my $key ( sort ( keys %R1R2_diff_lengths ) ) {
		print $dist_OUT "R1-R2", "\t", "$key", "\t", "$R1R2_diff_lengths{$key}", "\n";
	}
	print $dist_OUT "\n\n";
	
	print $dist_OUT "Distribution of singleton R1 read lengths:\n";
	foreach my $key ( sort ( keys %R1_singles_lengths ) ) {
		print $dist_OUT "R1_singles", "\t", "$key", "\t", "$R1_singles_lengths{$key}", "\n";
	}
	print $dist_OUT "\n\n";
	close $dist_OUT;
	
	print "\nData can be found in:\n$R1out_file\n$R2out_file\n$R1singles_file\n$read_dist_file\n$summary_file\n";

}

print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics 17:29. DOI 10.1186/s12859-016-0879-y.\n\n";

### SUB-ROUTINES ###

sub hd {
	my @var = $_;
	return length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );
}

exit;
