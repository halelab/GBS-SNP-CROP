#!/usr/bin/perl

###############################################################################################################################
# GBS-SNP-CROP, Step 1. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
###############################################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Zlib;

my $Usage = "Usage: perl GBS-SNP-CROP-1.pl -d <data type, PE = Paired-End or SE = Single-End> -b <barcode-ID file name> -fq <FASTQ file name seed>\n" 
."-s <start number of FASTQ file> -e <end number of FASTQ file> -enz1 <enzyme 1 restriction site residue sequence>\n" 
."-enz2 <enzyme 2 restriction site residue sequence>.\n";
my $Manual = "Please see the User Manual on the GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y.\n"; 

my ($dataType,$barcodesID_file,$fastq_seed,$fastq_start_num,$fastq_end_num,$enzyme1,$enzyme2);

GetOptions(
'd=s' => \$dataType,          # string - "PE" or "SE"
'b=s' => \$barcodesID_file,   # file
'fq=s' => \$fastq_seed,       # string
's=s' => \$fastq_start_num,   # numeric 
'e=s' => \$fastq_end_num,     # numeric
'enz1=s' => \$enzyme1,        # string
'enz2=s' => \$enzyme2,        # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 1, v.2.0\n#################################\n";

my $RC_enz1 = reverse $enzyme1;
$RC_enz1 =~ tr/acgtACGT/tgcaTGCA/;
my $RC_enz2 = reverse $enzyme2;
$RC_enz2 =~ tr/acgtACGT/tgcaTGCA/;

sub main {
	my $dir1 = "summaries";
	my $dir2 = "singles"; 
	my $dir3 = "distribs"; 
	my $dir4 = "parsed";
	unless(-e $dir1, or mkdir $dir1) {die "Directory $dir1 cannot be created.\n";}
	unless(-e $dir2, or mkdir $dir2) {die "Directory $dir2 cannot be created.\n";}
	unless(-e $dir3, or mkdir $dir3) {die "Directory $dir3 cannot be created.\n";}
	unless(-e $dir4, or mkdir $dir4) {die "Directory $dir4 cannot be created.\n";}
}
main();

############################
# Parsing Paired-End data 
############################

if ($dataType eq "PE") {
	print "Parsing Paired-End reads...\n";

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
		
		# Setting barcode ID, input and output files
		my @barcodes = ();
		my $barcode_string = "";

		open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
		while ( <$BAR> ) {
			my $barcodeID = $BAR;
			chomp $barcodeID;
			my @barcode = split("\t", $barcodeID);
			my $barcode_list = $barcode[0];
			push @barcodes, $barcode_list;
			$barcode_string = "$barcode_string $_";
		}
		close $BAR;
		chomp (@barcodes);

		$barcode_string = "$barcode_string ";
		$barcode_string =~ s/\n/ /g;

		my $total_raw_reads = 0;
		my $no_R1_RE_site_tally = 0;
		my $no_R1_barcode_tally = 0;
		
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

		# Start parsing
		print "\nProcessing files $R1file and $R2file\n";

		my @R1read;
		my @R2read;

		open my $IN1, '-|', 'gzip', '-dc', $R1file or die "Can't open file $R1file: $!\n";
		open my $IN2, '-|', 'gzip', '-dc', $R2file or die "Can't open file $R2file: $!\n";

		my $i = 1;
		my $progress = 0;

		while (!eof($IN1) and !eof($IN2)) {

			$progress++;
			if ($progress % 100 == 0) {
				print "\033[JStatus: ${progress} reads processed..."."\033[G";
			}

			my $file1 = <$IN1>;
			my $file2 = <$IN2>;
			if ($i % 4 != 0) {
				push @R1read, $file1;
				push @R2read, $file2;
				$i++;
			} else {
				push @R1read, $file1;
				push @R2read, $file2;
				chomp (@R1read);
				chomp (@R2read);
				$total_raw_reads++;
				$i++;

				# Consolidating R1 and R2 reads
				if ( $R1read[0] =~ /^(@.*:N:0:)\w{0,50}$/ ) {
					$R1read[0] = $1;
				} else {
					$R1read[0] = $R1read[0];
				}
				$R2read[0] = "";

				if ( $R2read[2] =~ /^#+$/ ) {
					$R2read[1] = "";
					$R2read[2] = "";
				} else {
					$R2read[1] = $R2read[1];
					$R2read[2] = $R2read[2];
				}

				# Extracting barcodes from R1 reads (exact matches or one mismatch)
				if ( $R1read[1] =~ /^(\w{3,10})$enzyme1(\w+)$/ ) { EXIT_IF: {
					my $inline_index = $1;
					$R1read[1] = $2;
					my $trim_length = length ($1) + length ($enzyme1);
					$R1read[3] = substr ( $R1read[3], $trim_length, $trim_length + length ($2) );

					if ( $barcode_string =~ /\s$inline_index\s/ ) {
						$R1read[0] = "$R1read[0]$inline_index";
						$R1read[4] = "$inline_index";
						last EXIT_IF;
					} else {
						$R1read[4] = "";
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
						$R1read[0] = "$R1read[0]$barcode_match";
						$R1read[4] = "$barcode_match";
					} else {
						$no_R1_barcode_tally++;
						$R1read[0] = "";
						$R1read[1] = "";
						$R1read[2] = "";
						$R1read[3] = "";
						$R1read[4] = "";
						$R2read[0] = "";
						$R2read[1] = "";
						$R2read[2] = "";
						$R2read[3] = "";
					}
				}} else {
					$no_R1_RE_site_tally++;
					$R1read[0] = "";
					$R1read[1] = "";
					$R1read[2] = "";
					$R1read[3] = "";
					$R1read[4] = "";
					$R2read[0] = "";
					$R2read[1] = "";
					$R2read[2] = "";
					$R2read[3] = "";
				}

				# Trimming R1 reads based on presence of $enzyme2 site residue, barcodes, and Illumina tag
				if ( $R1read[1] =~ /^(\w+)$RC_enz2[AN][GN][AN][TN][CN][GN][GN][AN][AN]\w+$/ ) {
					my $R1_read = $1;
					my $R1_read_length = length ( $R1_read );
					$R1read[1] = $R1_read;
					$R1read[3] = substr ( $R1read[3], 0, $R1_read_length );
				}

				# Trimming R2 reads based on presence of $enzyme1 site residue, barcodes, and Illumina tag
				if ( $R2read[1] =~ /^(\w+)$RC_enz1(\w{3,10})[AN][GN][AN][TN]/ ) {
					my $R2_read = $1;
					my $inline_index = scalar reverse $2;
					$inline_index =~ tr/actgACTG/tgacTGAC/;
					my $index_hd = hd ( $inline_index, $R1read[4] ); 
					if ( $index_hd <= 1 ) {
						$R2read[1] = $R2_read;
						my $R2_read_length = length ($R2read[1]);
						$R2read[3] = substr ( $R2read[3], 0, $R2_read_length );
					} 
				}

				# Trimming/culling R2 reads based on quality/length
				my $R2_min_length_tally = 0;
				my $dN_R2_tally = 0;
				if ( ( $R2read[1] ne "" ) && ( $R2read[2] ne "" ) ) {
					if ( $R2read[2] =~ /(#{4,})$/ ) {
						my $bad_bases = $1;
						my $trimmed_read_length = length ($R2read[3]) - length ($bad_bases);
						$R2read[1] = substr ( $R2read[1], 0, $trimmed_read_length );
						$R2read[3] = substr ( $R2read[3], 0, $trimmed_read_length );
					}
					my @R2_seq = split ("", $R2read[1]);
					my $N_count = 0;
					foreach my $base (@R2_seq) {
						if ( $base eq "N" ) {
							$N_count++;
						}
					}
					if ( $N_count >= (length ($R2read[1]) / 2 )) {
						$dN_R2_tally++;
						$R2read[1] = "";
						$R2read[2] = "";
					}
				}

				# Writing results to output files
				if ( ( $R1read[1] ne "" ) && ( $R2read[1] eq "" ) ) {
					$R1_singleton_tally++;
					print $R1_singles "$R1read[0]\n";
					print $R1_singles "$R1read[1]\n";
					print $R1_singles "+\n";
					print $R1_singles "$R1read[3]\n";
					my $R1_len = length ($R1read[1]);

					if ( exists $R1_singles_lengths{$R1_len} ) {
						$R1_singles_lengths{$R1_len}++;
					} else {
						$R1_singles_lengths{$R1_len} = 1;
					}

				} elsif (( $R1read[1] ne "" ) && ( $R2read[1] ne "" )) {
					$usable_pairs_tally++;
					print $R1_OUT "$R1read[0]\n";
					print $R1_OUT "$R1read[1]\n";
					print $R1_OUT "+\n";
					print $R1_OUT "$R1read[3]\n";

					my $R1_len = length ($R1read[1]);
					if ( exists $R1_lengths{$R1_len} ) {
						$R1_lengths{$R1_len}++;
					} else {
						$R1_lengths{$R1_len} = 1;
					}

					my $R2_header = $R1read[0];
					$R2_header =~ s/ 1:/ 2:/;
					print $R2_OUT "$R2_header\n";
					print $R2_OUT "$R2read[1]\n";
					print $R2_OUT "+\n";
					print $R2_OUT "$R2read[3]\n";

					my $R2_len = length ($R2read[1]);
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
				}
				@R1read = ();
				@R2read = ();
			}
		}
		close $IN1;
		close $IN2;
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


############################
# Parsing Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Parsing Single-End reads...\n";

	rmdir "singles";

	for ( my $file_index = $fastq_start_num; $file_index <= $fastq_end_num; $file_index++ ) {

		my $R1file = "";
		my $fileseed = "";

		if ( length ($file_index) == 1 ) {
			$R1file = join ("", "$fastq_seed","_R1_00", "$file_index", ".fastq.gz");
			$fileseed = join ("", "$fastq_seed","_00", "$file_index");
		} elsif ( length ($file_index) == 2 ) {
			$R1file = join ("", "$fastq_seed","_R1_0", "$file_index", ".fastq.gz");
			$fileseed = join ("", "$fastq_seed","_0", "$file_index");
		} elsif ( length ($file_index) == 3 ) {
			$R1file = join ("", "$fastq_seed","_R1_", "$file_index", ".fastq.gz");
			$fileseed = join ("", "$fastq_seed","_", "$file_index");
		} else {
			die "Unable to identify input fastq.gz files.\n";
		}

		# Setting barcode ID, input and output files
		my @barcodes = ();
		my $barcode_string = "";

		open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
		while ( <$BAR> ) {
			my $barcodeID = $BAR;
			chomp $barcodeID;
			my @barcode = split("\t", $barcodeID);
			my $barcode_list = $barcode[0];
			push @barcodes, $barcode_list;
			$barcode_string = "$barcode_string $_";
		}
		close $BAR;
		chomp (@barcodes);

		$barcode_string = "$barcode_string ";
		$barcode_string =~ s/\n/ /g;
		
		my $total_raw_reads = 0;
		my $no_R1_RE_site_tally = 0;
		my $no_R1_barcode_tally = 0;
		
		my $R1out_file = join (".","$fileseed","R1parsed.fastq");
		my $summary_file = join (".","$fileseed","summary");
		my $read_dist_file = join (".","$fileseed","dist");

		my $R1_tally = 0;
		my $usable_pairs_tally = 0;

		open my $R1_output, " | gzip > ./parsed/$R1out_file.gz" or die "Can't open $R1out_file\n";

		my %R1_singles_lengths;
		
		# Start parsing
		print "\nProcessing $R1file file\n";

		my @R1read;

		open my $IN1, '-|', 'gzip', '-dc', $R1file or die "Can't open file $R1file: $!\n";
		
		my $i = 1;
		my $progress = 0;
		
		while(<$IN1>) {
		
		$progress++;
		if ($progress % 100 == 0) {
			print "\033[JStatus: ${progress} reads processed..."."\033[G";
		}
			if ($i % 4 != 0) {
				push @R1read, $_;
				$i++;
			} else {
				push @R1read, $_;
				chomp (@R1read);
				$total_raw_reads++;
				$i++;

				# Consolidating reads. Checking header and barcode
				if ( $R1read[0] =~ /^(@.*:N:0:)\w{0,50}$/ ) {
					$R1read[0] = $1;
				} else {
					$R1read[0] = $R1read[0];
				}
				$R1read[2] = $R1read[3];

				# Extracting barcodes from reads (exact matches or one mismatch)
				if ( $R1read[1] =~ /^(\w{3,10})$enzyme1(\w+)$/ ) { EXIT_IF: {
					my $inline_index = $1;
					$R1read[1] = $2;
					my $trim_length = length ($1) + length ($enzyme1);
					$R1read[2] = substr ( $R1read[2], $trim_length, $trim_length + length ($2) );

					if ( $barcode_string =~ /\s$inline_index\s/ ) {
						$R1read[0] = "$R1read[0]$inline_index";
						$R1read[3] = "$inline_index";
						last EXIT_IF;
					} else {
						$R1read[3] = "";
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
						$R1read[0] = "$R1read[0]$barcode_match";
						$R1read[3] = "$barcode_match";
					} else {
						$no_R1_barcode_tally++;
						$R1read[0] = "";
						$R1read[1] = "";
						$R1read[2] = "";
						$R1read[3] = "";
					}
				}} else {
					$no_R1_RE_site_tally++;
					$R1read[0] = "";
					$R1read[1] = "";
					$R1read[2] = "";
					$R1read[3] = "";
				}

				# Trimming reads based on presence of $enzyme2 site residue, barcodes, and Illumina tag
				if ( $R1read[1] =~ /^(\w+)$RC_enz2[AN][GN][AN][TN][CN][GN][GN][AN][AN]\w+$/ ) {
					my $R1_read = $1;
					my $R1_read_length = length ( $R1_read );
					$R1read[1] = $R1_read;
					$R1read[2] = substr ( $R1read[2], 0, $R1_read_length );
				}

				# Writing singletons results to output files
				if ( $R1read[1] ne "" ) {
					$R1_tally++;
					print $R1_output "$R1read[0]\n";
					print $R1_output "$R1read[1]\n";
					print $R1_output "+\n";
					print $R1_output "$R1read[2]\n";

					my $R1_len = length ($R1read[1]);
					if ( exists $R1_singles_lengths{$R1_len} ) {
						$R1_singles_lengths{$R1_len}++;
					} else {
						$R1_singles_lengths{$R1_len} = 1;
					}
				}
				@R1read = ();
			}
		}
		close $IN1;
		close $R1_output;

		open my $summary_OUT, ">", "./summaries/$summary_file" or die "Can't open $summary_file\n";

		print $summary_OUT "SUMMARY:\n";
		print $summary_OUT "Total raw reads = $total_raw_reads\n";
		print $summary_OUT "Reads with no identifiable restriction site = $no_R1_RE_site_tally\n";
		print $summary_OUT "Reads with no identifiable barcode = $no_R1_barcode_tally\n";
		print $summary_OUT "Usable reads = $R1_tally\n\n";
		my $R1_singleton_percent = ( $R1_tally / $total_raw_reads ) * 100;
		print $summary_OUT "Percentage of usable reads = $R1_singleton_percent\n";
		close $summary_OUT;

		open my $dist_OUT, '>', "./distribs/$read_dist_file" or die "Can't open $read_dist_file\n";

		print $dist_OUT "Distribution of reads lengths:\n";
		foreach my $key ( sort ( keys %R1_singles_lengths ) ) {
			print $dist_OUT "R1_singles", "\t", "$key", "\t", "$R1_singles_lengths{$key}", "\n";
		}
		print $dist_OUT "\n\n";
		close $dist_OUT;

		print "\nData can be found in:\n$R1out_file\n$read_dist_file\n$summary_file\n";
	}
}

print "\n\nPlease cite: Melo et al. GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics (2016) 17:29 DOI 10.1186/s12859-016-0879-y..\n\n";

### SUB-ROUTINES ###

sub hd {
	my @var = $_;
	return length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] );
}

exit;
