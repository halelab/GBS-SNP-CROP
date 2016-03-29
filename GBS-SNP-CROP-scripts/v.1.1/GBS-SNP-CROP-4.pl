#!/usr/bin/perl

###########################################################################################################################
# GBS-SNP-CROP, Step 4. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y.
###########################################################################################################################

##########################################################################################
# Requirement 1: PEAR bioinformatics tool (Zhang et al., 2014)
# Requirement 2: USEARCH bioinformatics tool (Edgar, 2010)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-4.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file> -rl <Raw GBS read lengths> -pl <minimum length required after merging to retain read>\n"
."-p <p-value for PEAR (Zhang et al., 2014)> -id <nucleotide identity value required for USEARCH (Edgar, 2010) read clustering>\n"
." -t <number of threads dedicated to USEARCH clustering> -MR <Mock reference name>.\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

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

print "\n#################################\n# GBS-SNP-CROP, Step 4, v.1.1\n#################################\n";

my $MR_Genome = "";
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

		my $stitched_file = join (".", "$file","stitched","fastq");

		open my $stitched_OUT, ">", "$stitched_file" or die "Can't open $stitched_file: $!\n";

		my $stitched_tally = 0;
		my $unstitched_tally = 0;
		my $seq_stitch = "A" x 20;
		my $qual_stitch = "I" x 20;

		for (my $k = 0; $k <= $size; $k++) {
			my $R1_length = length($R1reads[$k][1]);
			my $R2_length = length($R2reads[$k][1]);
			if ( $R1_length < ($raw_seq_length - 19) || $R2_length < ($raw_seq_length - 5) ) {
				$unstitched_tally++;
				next;
			} else {
				my $assembled_seq = join ("", "$R1reads[$k][1]", "$seq_stitch", "$R2reads[$k][1]");
				my $assembled_qual = join ("", "$R1reads[$k][3]", "$qual_stitch", "$R2reads[$k][3]");
				print $stitched_OUT "$R1reads[$k][0]\n";
				print $stitched_OUT "$assembled_seq\n";
				print $stitched_OUT "+\n";
				print $stitched_OUT "$assembled_qual\n";
				$stitched_tally++;
				next;
			}
		}
		close $stitched_OUT;

		my $unstitched_percentage = ( $unstitched_tally / ( $stitched_tally + $unstitched_tally ) ) * 100;

		print $code_OUT "\nFor paired of unassembled files $R1file2 and $R2file2:\n";
		print $code_OUT "Total number of stitched read pairs = $stitched_tally\n";
		print $code_OUT "Total number of unstitchable read pairs = $unstitched_tally\n";
		print $code_OUT "Percent of unassembled read pairs that were unstitchable = $unstitched_percentage\n\n";
	}
	system ( "rm *.unassembled*" );

# 3. Concatenate the merged and stitched reads for each genotype into a single file for use in building the Mock Reference 

	foreach my $file (@MR_taxa_files) {
		my $assembled = join (".", "$file","assembled","fastq");
		my $stitched = join (".", "$file","stitched","fastq");
		my $out = join(".","$file","AssembledStitched","fastq");	

		print "Concatenating $assembled and $stitched files...\n";
		system ( "cat $assembled $stitched > $out" );
	}
	system ( "rm *.assembled.fastq *.stitched.fastq" );

# 4. Use USEARCH to cluster reads within each genotype 

# 4A. Sorting the FASTQ sequences in order of decreasing length (bp)
	print "\nSorting the FASTQ sequences in order of decreasing length...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4A = join (".", "$file","AssembledStitched","fastq");
		my $out4A = join(".","$file","sorted_by_length","fasta");
		print "Processing file $in4A...\n";
		system ( "usearch -sortbylength $in4A -fastaout $out4A" );
	}

# 4B. Finding initial clusters (centroids)
	print "\nFinding genotype-specific initial clusters (centroids)...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4B = join(".","$file","sorted_by_length","fasta");
		my $out4B = join(".","$file","clusters","fasta");
		print "Processing file $in4B...\n";
		system ( "usearch -cluster_fast $in4B -id $id -threads $threads -consout $out4B -sizeout" );
	}

# 4C. Sorting genotype-specific initial clusters (centroids) by depth 
	print "\nSorting genotype-specific initial clusters (centroids) by depth...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4C = join(".","$file","clusters","fasta");
		my $out4C = join(".","$file","sorted_by_size","fasta");
		print "Processing file $in4C...\n";
		system ( "usearch -sortbysize $in4C -fastaout $out4C" );
	}

# 4D. Reclustering the genotype-specific centroids
	print "\nRecluster the genotype-specific centroids...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4D = join(".","$file","sorted_by_size","fasta");
		my $out4D = join(".","$file","reclusters","fasta");
		print "Processing file $in4D...\n";
		system ( "usearch -cluster_fast $in4D -id $id -threads $threads -consout $out4D" );
	}

	print "\nAll sub-steps for clustering reads within genotype(s) were completed!\n";

# Define the Mock reference

	my $MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fasta");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fasta");
	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

# 5. Use USEARCH to cluster genotype-specific centroids across the population 

	if (scalar (@MR_taxa_files) > 1) {

		print "\n\nConcatenating the genotype-specific non-redundant centroids...";
		my $in5A = join(".",$MockRefName,"consensus","fasta");
		system ( "cat *reclusters.fasta > $in5A" );
	
# 5A. Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $in5A -fastaout $out5A" );
	
# 5B. Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");
	
# 5C. Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );
	
# 5D. Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		my $out5D = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $out5D" );
	
		print "All sub-steps for clustering population-level centroids were completed!\n";

# 6. Delete N-containing centroids and create the Mock Reference

		open my $IN5, "<", "$out5D" or die "Unable to open $out5D file\n";
		{
			local $/ = ">";
			my @file = <$IN5>;
			shift @file;
			chomp @file;
			foreach my $line(@file) {
				my ($name, $seq) = split /\n/, $line, 2;
				$seq =~ s/\n//g;
				$seq = uc $seq;
				if ($seq =~ "N") { 
					next;
				} else {
					print $OUT1 ">$name\n$seq\n";
					print $OUT2 "$seq";
				}
			}
		}
		close $IN5;

	} elsif (scalar (@MR_taxa_files) == 1) {

		foreach my $file (@MR_taxa_files) {
			my $out4D = join(".","$file","reclusters","fasta");
			open my $IN4, "<", "$out4D" or die "Unable to open $out4D file\n";
			{
				local $/ = ">";
				my @file = <$IN4>;
				shift @file;
				chomp @file;
				foreach my $line(@file) {
					my ($name, $seq) = split /\n/, $line, 2;
					$seq =~ s/\n//g;
					$seq = uc $seq;
					if ($seq =~ "N") { 
						next;
					} else {
						print $OUT1 ">$name\n$seq\n";
						print $OUT2 "$seq";
					}
				}
			}
			close $IN4;
		} 
	}

	sub main {
	my $dir = "fastqForRef"; 
	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
	}
	main();
	system ( "mv *.AssembledStitched.fastq ./fastqForRef" );
	system ( "rm  *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta *.reclusters.fasta *.consensus.fasta" );
	
############################
# Parsing Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Parsing Single-End reads...\n";

# 4. Use USEARCH to cluster reads within each genotype 

# 4A. Sorting the FASTQ sequences in order of decreasing length (bp)
	print "\nSorting the FASTQ sequences in order of decreasing length...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4A = join (".", "$file","R1","fastq");
		my $out4A = join(".","$file","sorted_by_length","fasta");
		print "Processing file $in4A...\n";
		system ( "usearch -sortbylength $in4A -fastaout $out4A" );
	}

# 4B. Finding initial clusters (centroids)
	print "\nFinding genotype-specific initial clusters (centroids)...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4B = join(".","$file","sorted_by_length","fasta");
		my $out4B = join(".","$file","clusters","fasta");
		print "Processing file $in4B...\n";
		system ( "usearch -cluster_fast $in4B -id $id -threads $threads -consout $out4B -sizeout" );
	}

# 4C. Sorting genotype-specific initial clusters (centroids) by depth 
	print "\nSorting genotype-specific initial clusters (centroids) by depth...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4C = join(".","$file","clusters","fasta");
		my $out4C = join(".","$file","sorted_by_size","fasta");
		print "Processing file $in4C...\n";
		system ( "usearch -sortbysize $in4C -fastaout $out4C" );
	}

# 4D. Reclustering the genotype-specific centroids
	print "\nRecluster the genotype-specific centroids...\n";
	foreach my $file (@MR_taxa_files) {
		my $in4D = join(".","$file","sorted_by_size","fasta");
		my $out4D = join(".","$file","reclusters","fasta");
		print "Processing file $in4D...\n";
		system ( "usearch -cluster_fast $in4D -id $id -threads $threads -consout $out4D" );
	}

	print "\nAll sub-steps for clustering reads within genotype(s) were completed!\n";

# Define Mock reference

	my $MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fasta");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fasta");
	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

# 5. Use USEARCH to cluster genotype-specific centroids across the population 

	if (scalar (@MR_taxa_files) > 1) {

		print "\n\nConcatenating the genotype-specific non-redundant centroids...";
		my $in5A = join(".",$MockRefName,"consensus","fasta");
		system ( "cat *reclusters.fasta > $in5A" );
	
# 5A. Sorting the full set of genotype-specific centroids in order of decreasing length
		print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
		my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
		system ( "usearch -sortbylength $in5A -fastaout $out5A" );
	
# 5B. Finding population-level initial clusters (centroids)
		print "\nFinding population-level initial clusters (centroids)...\n";
		my $out5B = join (".",$MockRefName,"clusters","fasta");
		system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");
	
# 5C. Sorting population-level initial clusters (centroids) by depth
		print "\nSorting population-level initial clusters (centroids) by depth...\n";
		my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
		system ( "usearch -sortbysize $out5B -fastaout $out5C" );
	
# 5D. Reclustering the population-level centroids
		print "\nReclustering the population-level centroids...\n";
		my $out5D = join (".",$MockRefName,"reclusters","fasta");
		system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $out5D" );
	
		print "All sub-steps for clustering population-level centroids were completed!\n";

# 6. Delete N-containing centroids and create the Mock Reference

		open my $IN5, "<", "$out5D" or die "Unable to open $out5D file\n";
		{
			local $/ = ">";
			my @file = <$IN5>;
			shift @file;
			chomp @file;
			foreach my $line(@file) {
				my ($name, $seq) = split /\n/, $line, 2;
				$seq =~ s/\n//g;
				$seq = uc $seq;
				if ($seq =~ "N") { 
					next;
				} else {
					print $OUT1 ">$name\n$seq\n";
					print $OUT2 "$seq";
				}
			}
		}
		close $IN5;

	} elsif (scalar (@MR_taxa_files) == 1) {

		foreach my $file (@MR_taxa_files) {
			my $out4D = join(".","$file","reclusters","fasta");
			open my $IN4, "<", "$out4D" or die "Unable to open $out4D file\n";
			{
				local $/ = ">";
				my @file = <$IN4>;
				shift @file;
				chomp @file;
				foreach my $line(@file) {
					my ($name, $seq) = split /\n/, $line, 2;
					$seq =~ s/\n//g;
					$seq = uc $seq;
					if ($seq =~ "N") { 
						next;
					} else {
						print $OUT1 ">$name\n$seq\n";
						print $OUT2 "$seq";
					}
				}
			}
			close $IN4;
		} 
	}
	system ( "rm  *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta *.reclusters.fasta *.consensus.fasta" );
}

print "\nCongratulations! Your '$MockRefName' Mock Reference genome was assembled using the following genotypes:\n@MR_taxa_files\n";
print "We recommend the using '$MR_Genome' as a reference genome for mapping your parsed, high quality reads.";
print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
