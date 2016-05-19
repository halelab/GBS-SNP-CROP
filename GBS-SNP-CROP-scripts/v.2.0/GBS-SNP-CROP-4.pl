#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 4. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
#########################################################################################################################

##########################################################################################
# Requirement 1: PEAR bioinformatics tool (Zhang et al., 2014)
# Requirement 2: USEARCH bioinformatics tool (Edgar, 2010)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw/shuffle/;

my $Usage = "Usage: perl GBS-SNP-CROP-4.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file> -rl <Raw GBS read lengths> -pl <minimum length required after merging to retain read>\n"
."-p <p-value for PEAR (Zhang et al., 2014)> -id <nucleotide identity value required for USEARCH (Edgar, 2010) read clustering>\n"
." -t <number of threads dedicated to USEARCH clustering> -MR <Mock reference name>.\n";
my $Manual = "Please see the User Manual on the GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($dataType,$barcodesID_file,$raw_seq_length,$pear_length,$pvalue,$id,$threads,$MockRefName);

GetOptions(
'd=s' => \$dataType,          # string - "PE" or "SE"
'b=s' => \$barcodesID_file,   # file
'rl=s' => \$raw_seq_length,   # numeric
'pl=s' => \$pear_length,      # numeric
'p=s' => \$pvalue,            # numeric
'id=s' => \$id,               # numeric
't=s' => \$threads,           # numeric
'MR=s' => \$MockRefName,      # string
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 4, v.2.0\n#################################\n";

my @MR_taxa_files = ();

open my $BAR, "<", "$barcodesID_file" or die "Can't find barcode_ID file\n";
while ( <$BAR> ) {
	my $barcodeID = $_;
	chomp $barcodeID;
	my @barcode = split("\t", $barcodeID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	my $MR_geno = $barcode[2];

	if ($MR_geno eq "YES") {
		push @MR_taxa_files, $TaxaNames;
	}
}
close $BAR;
chomp (@MR_taxa_files);

sub main {
		my $dir = "fastaForRef"; 
		unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
	}
	main();

my $MR_Cluster = " ";
my $MR_Genome = " ";

############################
# Parsing Paired-End data 
############################

if ($dataType eq "PE") {
	print "Parsing Paired-End reads...\n\n";
	
	my $script_output = "Pear.log";
	open my $code_OUT, ">", "$script_output" or die "Can't open $script_output\n";
	print $code_OUT "### PEAR (Zhang et al., 2014) summary results:\n";
	
# 1. Use PEAR to merge the parsed R1 and R2 reads to create single reads, if possible

	foreach my $file (@MR_taxa_files) {
		my $R1input1 = join (".", "$file","R1","fastq");
		my $R2input2 = join (".", "$file","R2","fastq");
		my $Pear_out = join("","$file");
	
		print "Assembling paired $R1input1 and $R2input2 reads using PEAR...\n";
		print $code_OUT `pear -f $R1input1 -r $R2input2 -o $Pear_out -p $pvalue -n $pear_length -j $threads`;

		# Transform FASTQ PEAR assembly output file into FASTA format
		my $Ain =  join (".", "$file","assembled","fastq");
		my $Aout = join (".", "$file","assembled","fasta");
		open my $AIN, "<", "$Ain" or die "Can't open $Ain: $!\n";
		open my $AOUT, ">", "$Aout" or die "Can't load $Aout: $!\n";

		my @R1read;
		my @R1reads;
		
		my $i = 1;
		while(<$AIN>) {
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
		unlink $Ain;
		close $AIN;
	
		my $size = scalar @R1reads - 1;
		for (my $k = 0; $k <= $size; $k++) {
	
			$R1reads[$k][0] =~ s/@/>/g;
			$R1reads[$k][0] =~ s/ /:/g;
			print $AOUT "$R1reads[$k][0]\n$R1reads[$k][1]\n";
		}
		close $AOUT;
	}
	system ( "rm *.discarded.fastq" );

	print "Assembly completed for all mergeable read pairs.\n\n";


# 2. Stitch unassembled R1 and R2 reads together with an intermediate run of 20 high-quality A's


	print "Manually stitching together unassembled reads...\n\n";

	print $code_OUT "\n\n### Manually stitching together unassembled reads results:\n";

	foreach my $file (@MR_taxa_files) {
		my $R1file2 = join (".", "$file","unassembled","forward","fastq");
		my $R2file2 = join (".", "$file","unassembled","reverse","fastq");

		my @R1read;
		my @R1reads;

		open my $IN1, "<", "$R1file2" or die "Can't open $R1file2: $!\n";

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

		open my $IN2, "<", "$R2file2" or die "Can't open $R2file2: $!\n";

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

		my $size = scalar @R1reads - 1;

		my $stitched_file = join (".", "$file","stitched","fasta");

		open my $stitched_OUT, ">", "$stitched_file" or die "Can't open $stitched_file: $!\n";

		my $stitched_tally = 0;
		my $unstitched_tally = 0;
		my $seq_stitch = "A" x 20;

		for (my $k = 0; $k <= $size; $k++) {
			my $R1_length = length($R1reads[$k][1]);
			my $R2_length = length($R2reads[$k][1]);
			if ( $R1_length < ($raw_seq_length - 19) || $R2_length < ($raw_seq_length - 5) ) {
				$unstitched_tally++;
				next;
			} else {
				my $assembled_seq = join ("", "$R1reads[$k][1]", "$seq_stitch", "$R2reads[$k][1]");
				$R1reads[$k][0] =~ s/@/>/g;
				$R1reads[$k][0] =~ s/ /:/g;
				print $stitched_OUT "$R1reads[$k][0]\n";
				print $stitched_OUT "$assembled_seq\n";
				$stitched_tally++;
				next;
			}
		}
		close $stitched_OUT;

		my $unstitched_percentage = ( $unstitched_tally / ( $stitched_tally + $unstitched_tally ) ) * 100;

		print $code_OUT "\nFor pair of unassembled files $R1file2 and $R2file2:\n";
		print $code_OUT "Total number of stitched read pairs = $stitched_tally\n";
		print $code_OUT "Total number of unstitchable read pairs = $unstitched_tally\n";
		print $code_OUT "Percent of unassembled read pairs that were unstitchable = $unstitched_percentage\n\n";
	}
	system ( "rm *.unassembled*" );

# 3. Concatenate the merged and stitched reads for each genotype into a single FASTA file for use in building the Mock Reference 

	foreach my $file (@MR_taxa_files) {
		my $assembled = join (".", "$file","assembled","fasta");
		my $stitched = join (".", "$file","stitched","fasta");
		my $out = join(".","$file","AssembledStitched","fasta");

		print "Concatenating $assembled and $stitched files...\n";
		system ( "cat $assembled $stitched > $out" );
	}

# 4. Estimate the proportion of reads to be sampled for USEARCH clustering

	my $sum = 0;
	my @files_size = ();
	my $UsearchOUT;

	foreach my $file (@MR_taxa_files) {
		my $A_S_in = join(".","$file","AssembledStitched","fasta");
		$sum = ( $sum + (-s "$A_S_in") );
	}

	if ($sum <= 0) {
		print "GBS-SNP-CROP cannot access the sizes of the files selected to build the Mock Reference and so cannot calculate the proportion of those files to be used to randomly sample your data.\n"
		."Please refer to the GBS-SNP-CROP User Manual 'STAGE2: BUILD THE MOCK REFERENCE' for more details.\n\n";
		goto exit;
	}

	# Defining the value of sampling proportion based on the maximum memory allotment within USEARCH
	my $prop_sampling;
	my $S = ($sum * 5) + (8 * (65536)); #Based on USEARCH estimate
	if ($S <= 3500000000) {
		$prop_sampling = 1;
	} else {
		$prop_sampling = 3500000000 / $S;
		$prop_sampling = sprintf("%.3f", $prop_sampling);
	}

# 5. Use USEARCH to cluster reads

	# Without random sampling
	if ($prop_sampling >= 1) { # Sub-sampling is not required.
		print "\nYour Mock Reference will be assembled using all available reads from the designated genotypes.\n";

		my $UsearchIN = join(".","UsearchIN","fasta");
		system ( "cat *AssembledStitched.fasta > $UsearchIN" );

		# Step 5A: Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $UsearchIN -fastaout $out5A" );
	
		# Step 5B: Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");

		# Step 5C: Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );
	
		# Step 5D: Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		$UsearchOUT = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $UsearchOUT" );

		print "All sub-steps for clustering the population-level centroids were completed!\n";
		unlink $UsearchIN;

	# With random sampling
	} elsif ($prop_sampling < 1) {
		my $percent = $prop_sampling * 100;
		print "\nDue to the memory consumption limit of 4 Gb imposed by USEARCH (Zhang et al., 2014) and based on the sizes of your designated FASTA files, GBS-SNP-CROP will randomly sample $percent% of your reads for use in assembling the Mock Reference.\n\n";

		foreach my $file (@MR_taxa_files) {
			my $A_S_in = join(".","$file","AssembledStitched","fasta");
			my $sub = join (".","$file","AS","sub","fasta");

			print "Randomly sampling reads from $A_S_in file ...\n";

			# Loading data ...
			my @R1read;
			my @R1reads;

			open my $R1, "<", "$A_S_in" or die "Can't open $A_S_in file\n";
			my $i = 1;
			while(<$R1>) {
				if ($i % 2 != 0) {
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
			close $R1;

			# Creating an index array, the n sample and the sampled indexes
			my $total_reads = scalar @R1reads - 1;
			my @index = 0 .. $total_reads;
			my $n_sample = int ($prop_sampling * $total_reads);

			my @sample = (shuffle(@index))[0..($n_sample - 1)];
		
			# Grab sampled reads on outputs files
			open my $sub_out, ">", "$sub" or die "Can't load $sub file\n";

			foreach my $index (@sample) {
				print $sub_out "$R1reads[$index][0]\n$R1reads[$index][1]\n";
			}
			close $sub_out;
		}

		my $UsearchIN = join(".","UsearchIN","fasta");
		system ( "cat *AS.sub.fasta > $UsearchIN" );

		# Step 5A: Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $UsearchIN -fastaout $out5A" );

		# Step 5B: Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");

		# Step 5C: Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );

		# Step 5D: Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		$UsearchOUT = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $UsearchOUT" );

		print "All sub-steps for clustering population-level centroids were completed!\n";
		system ("rm *.AS.sub.fasta UsearchIN.fasta" );
	}

# 6. Building the mock reference, identifying poly-A coordinates for masking, and deleting all N-containing centroids

	$MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fasta");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fasta");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();

	open my $IN6, "<", "$UsearchOUT" or die "Unable to open $UsearchOUT file\n";
	{
		local $/ = ">";
		my @file = <$IN6>;
		shift @file;
		chomp @file;
		my $polyA = "A" x 20;
		foreach my $line(@file) {
			my ($name, $seq) = split /\n/, $line, 2;
			$seq = uc $seq;
			$seq =~ s/\n//g;
			if ($seq =~ "N") { 
				next;
			} else {
				my $seqG = join ("",$seq,$polyA); 	
				print $OUT1 ">$name\n$seq\n";
				print $OUT2 "$seqG";
				push @seqs, $seqG;
			}
		}
		my @bad_coords = ();
		my $mr = join '',@seqs;
		while ( $mr =~ m/(A{20,})/g ) { 
			my $end_pos = pos ($mr);
			my $length = length ($1);

			for (my $k = 1; $k <= $length; $k++) {
				my $mask_coord = $end_pos - $length + $k;
				print $OUT3 "$mask_coord\n";
			}
		}
	}
	close $IN6;

	system ( "mv *.AssembledStitched.fasta ./fastaForRef" );
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $UsearchOUT *assembled* *.stitched.fasta" );

############################
# Parsing Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Parsing Single-End reads...\n";

# 1. Transform FASTQ high quality reads into FASTA format
	foreach my $file (@MR_taxa_files) {
		my $R1input1 = join (".", "$file","R1","fastq");
		my $out = join (".", "$file","R1","fasta");
		open my $IN, "<", "$R1input1" or die "Can't open $R1input1\n";
		open my $OUT, ">", "$out" or die "Can't load $out\n";

		my @R1read;
		my @R1reads;

		my $i = 1;
		while(<$IN>) {
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
		close $IN;

		my $size = scalar @R1reads - 1;
		for (my $k = 0; $k <= $size; $k++) {

			$R1reads[$k][0] =~ s/@/>/g;
			$R1reads[$k][0] =~ s/ /:/g;
			print $OUT "$R1reads[$k][0]\n$R1reads[$k][1]\n";
		}
		close $OUT;
	}

# 2. Estimate the proportion of reads to be sampled for USEARCH clustering

	my $sum = 0;
	my @files_size = ();
	my $UsearchOUT;

	foreach my $file (@MR_taxa_files) {
		my $R1in = join(".","$file","R1","fasta");
		$sum = ( $sum + (-s "$R1in") );
	}

	if ($sum <= 0) {
		print "GBS-SNP-CROP cannot access the sizes of the files selected to build the Mock Reference and so cannot calculate the proportion of those files to be used to randomly sample your data.\n"
		."Please refer to the GBS-SNP-CROP User Manual 'STAGE2: BUILD THE MOCK REFERENCE' for more details.\n\n";
		goto exit;
	}

	# Defining the value of sampling proportion based on the maximum memory allotment within USEARCH
	my $prop_sampling;
	my $S = ($sum * 5) + (8 * (65536)); #Based on USEARCH estimate
	if ($S <= 3500000000) {
		$prop_sampling = 1;
	} else {
		$prop_sampling = 3500000000 / $S;
		$prop_sampling = sprintf("%.3f", $prop_sampling);
	}

# 3. Use USEARCH to cluster reads 

	# Without random sampling
	if ($prop_sampling >= 1) { # Sub-sampling is not required. 
		print "\nYour Mock Reference will be assembled using all available reads from the designated genotypes.\n";

		my $UsearchIN = join(".","UsearchIN","fasta");
		system ( "cat *.R1.fasta > $UsearchIN" );

		# Step 5A: Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $UsearchIN -fastaout $out5A" );

		# Step 5B: Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");

		# Step 5C: Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );

		# Step 5D: Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		$UsearchOUT = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $UsearchOUT" );

		print "All sub-steps for clustering population-level centroids were completed!\n";
		unlink $UsearchIN;

	# With random sampling
	} elsif ($prop_sampling < 1) {
		my $percent = $prop_sampling * 100;
		print "\nDue to the memory consumption limit of 4 Gb imposed by USEARCH (Zhang et al., 2014) and based on the sizes of your designated FASTA files, GBS-SNP-CROP will randomly sample $percent% of your reads for use in assembling the Mock Reference.\n\n";

		foreach my $file (@MR_taxa_files) {
			my $in = join(".","$file","R1","fasta");
			my $sub = join (".","$file","R1","sub","fasta");

			print "Randomly sampling reads from $in file ...\n";

			# Loading data ...
			my @R1read;
			my @R1reads;

			open my $R1, "<", "$in" or die "Can't open $in file\n";
			my $i = 1;
			while(<$R1>) {
				if ($i % 2 != 0) {
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
			close $R1;

			my $total_reads = scalar @R1reads - 1;
			my @index = 0 .. $total_reads;
			my $n_sample = int ($prop_sampling * $total_reads);

			my @sample = (shuffle(@index))[0..($n_sample - 1)];

			open my $sub_out, ">", "$sub" or die "Can't load $sub file\n";

			foreach my $index (@sample) {
				print $sub_out "$R1reads[$index][0]\n$R1reads[$index][1]\n";
			}
			close $sub_out;
		}

		my $UsearchIN = join(".","UsearchIN","fasta");
		system ( "cat *R1.sub.fasta > $UsearchIN" );

		# Step 5A: Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $UsearchIN -fastaout $out5A" );

		# Step 5B: Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");

		# Step 5C: Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );

		# Step 5D: Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		$UsearchOUT = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $UsearchOUT" );

		print "All sub-steps for clustering population-level centroids were completed!\n";
		system ("rm *.R1.sub.fasta UsearchIN.fasta" );
	}

# 4. Building the mock reference, identifying poly-A coordinates for masking, and deleting all N-containing centroids

	$MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fasta");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fasta");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();

	open my $IN4, "<", "$UsearchOUT" or die "Unable to open $UsearchOUT file\n";
	{
		local $/ = ">";
		my @file = <$IN4>;
		shift @file;
		chomp @file;
		my $polyA = "A" x 20;
		foreach my $line(@file) {
			my ($name, $seq) = split /\n/, $line, 2;
			$seq = uc $seq;
			$seq =~ s/\n//g;
			if ($seq =~ "N") { 
				next;
			} else {
				my $seqG = join ("",$seq,$polyA);
				print $OUT1 ">$name\n$seq\n";
				print $OUT2 "$seqG";
				push @seqs, $seqG;
			}
		}
		my @bad_coords = ();
		my $mr = join '',@seqs;
		while ( $mr =~ m/(A{20,})/g ) { 
			my $end_pos = pos ($mr);
			my $length = length ($1);

			for (my $k = 1; $k <= $length; $k++) {
				my $mask_coord = $end_pos - $length + $k;
				print $OUT3 "$mask_coord\n";
			}
		}
	}
	close $IN4;
	system ( "mv *.R1.fasta ./fastaForRef" );
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $UsearchOUT" );
}

print "\nCongratulations! Your '$MockRefName' Mock Reference genome was assembled using the following genotypes:\n@MR_taxa_files\n";

print "We recommend using '$MR_Genome' as the reference genome for mapping your parsed, high quality reads.";

print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
