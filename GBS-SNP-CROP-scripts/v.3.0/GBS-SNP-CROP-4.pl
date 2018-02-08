#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 4. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y.
#########################################################################################################################

##########################################################################################
# Requirement 1: PEAR v0.9.6 (Zhang et al., 2014)
# Requirement 2: VSEARCH v2.6.2 (Rognes et al., 2016)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw/shuffle/;

my $Usage = "Usage: perl GBS-SNP-CROP-4.pl -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file> -rl <Raw GBS read lengths> -pl <minimum length required after merging to retain read>\n"
."-p <p-value for PEAR (Zhang et al., 2014)> -id <nucleotide identity value required for VSEARCH (Rognes et al., 2016) read clustering>\n"
." -t <number of threads dedicated to VSEARCH clustering> -MR <Mock reference name>.\n";
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

print "\n#################################\n# GBS-SNP-CROP, Step 4, v.3.0\n#################################\n";
my $sttime = time;

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
		my $R1input1 = join (".", "$file","R1","fq", "gz");
		my $R2input2 = join (".", "$file","R2","fq","gz");
		my $Pear_out = join("","$file");
	
		print "Assembling paired $R1input1 and $R2input2 reads using PEAR...\n";
		print $code_OUT `pear -f $R1input1 -r $R2input2 -o $Pear_out -p $pvalue -n $pear_length -j $threads`;

		# Transforming FASTQ PEAR assembly output file into FASTA
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
	print "DONE.\n\n";

# 2. Stitch unassembled R1 and R2 reads together with an intermediate run of 20 high-quality A's
	print "Manually stitching together unassembled reads...\n";
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

		print $code_OUT "\nFor paired of unassembled files $R1file2 and $R2file2:\n";
		print $code_OUT "Total number of stitched read pairs = $stitched_tally\n";
		print $code_OUT "Total number of unstitchable read pairs = $unstitched_tally\n";
		print $code_OUT "Percent of unassembled read pairs that were unstitchable = $unstitched_percentage\n\n";
	}
	print "DONE.\n\n";
	system ( "rm *.unassembled*" );

# 3. Concatenate the merged and stitched reads for each genotype into a single file for use in building the Mock Reference 
	foreach my $file (@MR_taxa_files) {
		my $assembled = join (".", "$file","assembled","fasta");
		my $stitched = join (".", "$file","stitched","fasta");
		my $out = join(".","$file","AssembledStitched","fa");

		print "Concatenating $assembled and $stitched files...\n";
		system ( "cat $assembled $stitched > $out" );
	}
	print "DONE.\n";

# 4. Use VSEARCH to cluster reads

	my $VsearchIN = join(".","VsearchIN","fa");
	system ( "cat *AssembledStitched.fa > $VsearchIN" );

	# Step 4A: Sorting the full set of genotype-specific centroids in order of decreasing length
	print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
	my $out4A = join (".",$MockRefName,"sorted_by_length","fasta");
	system ( "vsearch -sortbylength $VsearchIN -output $out4A" );
	
	# Step 4B: Finding population-level initial clusters (centroids)
	print "\nFinding population-level initial clusters (centroids)...\n";
	my $out4B = join (".",$MockRefName,"clusters","fasta");
	system ( "vsearch -cluster_fast $out4A -id $id -threads $threads -consout $out4B -sizeout");

	# Step 4C: Sorting population-level initial clusters (centroids) by depth
	print "\nSorting population-level initial clusters (centroids) by depth...\n";
	my $out4C = join (".",$MockRefName,"sorted_by_size","fasta");
	system ( "vsearch -sortbysize $out4B -output $out4C" );
	
	# Step 4D: Reclustering the population-level centroids
	print "\nReclustering the population-level centroids...\n";
	my $VsearchOUT = join (".",$MockRefName,"reclusters","fasta");
	system ( "vsearch -cluster_fast $out4C -id $id -threads $threads -consout $VsearchOUT" );

	print "\nAll sub-steps for clustering population-level centroids were completed!\n";
	unlink $VsearchIN;
	
# 5. Defining THE Mock reference, masked file and Delete N-containing centroids in order to create the Mock Reference

	$MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fa");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fa");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();

	open my $IN6, "<", "$VsearchOUT" or die "Unable to open $VsearchOUT file\n";
	{
		local $/ = ">";
		my @file = <$IN6>;
		shift @file;
		chomp @file;
		my $polyA = "A" x 20;
		my $i = 1;
		foreach my $line(@file) {
			my ($name, $seq) = split /\n/, $line, 2;
			$seq = uc $seq;
			$seq =~ s/\n//g;
			if ($seq =~ "N") { 
				next;
			} else {
				my $seqG = join ("",$seq,$polyA);
				my $NAME = join("","Cluster",$i); 	
				print $OUT1 ">$NAME\n$seq\n";
				print $OUT2 "$seqG";
				push @seqs, $seqG;
				$i++;
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

	system ( "mv *.AssembledStitched.fa ./fastaForRef" );
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $VsearchOUT *assembled* *.stitched.fasta" );

############################
# Parsing Single-End data 
############################

} elsif ($dataType eq "SE") {
	print "Parsing Single-End reads...\n";

# 1. Transforming FASTQ high quality reads into FASTA reads
	foreach my $file (@MR_taxa_files) {
		my $R1input1 = join (".", "$file","R1","fq","gz");
		my $out = join (".", "$file","R1","fa");
		open my $IN, '-|', 'gzip', '-dc', $R1input1 or die "Can't open file $R1input1: $!\n";
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

# 2. Use VSEARCH to cluster reads 

	my $VsearchIN = join(".","VsearchIN","fa");
	system ( "cat *.R1.fa > $VsearchIN" );

	# Step 2A: Sorting the full set of genotype-specific centroids in order of decreasing length
	print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
	my $out2A = join (".",$MockRefName,"sorted_by_length","fasta");
	system ( "vsearch -sortbylength $VsearchIN -output $out2A" );

	# Step 2B: Finding population-level initial clusters (centroids)
	print "\nFinding population-level initial clusters (centroids)...\n";
	my $out2B = join (".",$MockRefName,"clusters","fasta");
	system ( "vsearch -cluster_fast $out2A -id $id -threads $threads -consout $out2B -sizeout");

	# Step 2C: Sorting population-level initial clusters (centroids) by depth
	print "\nSorting population-level initial clusters (centroids) by depth...\n";
	my $out2C = join (".",$MockRefName,"sorted_by_size","fasta");
	system ( "vsearch -sortbysize $out2B -output $out2C" );

	# Step 2D: Reclustering the population-level centroids
	print "\nReclustering the population-level centroids...\n";
	my $VsearchOUT = join (".",$MockRefName,"reclusters","fasta");
	system ( "vsearch -cluster_fast $out2C -id $id -threads $threads -consout $VsearchOUT" );

	print "\nAll sub-steps for clustering population-level centroids were completed!\n";
	unlink $VsearchIN;

# 3. Defining THE Mock reference, masked file and Delete N-containing centroids in order to create the Mock Reference

	$MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fa");
	$MR_Genome = join (".","$MockRefName","MockRef_Genome","fa");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";

	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();

	open my $IN4, "<", "$VsearchOUT" or die "Unable to open $VsearchOUT file\n";
	{
		local $/ = ">";
		my @file = <$IN4>;
		shift @file;
		chomp @file;
		my $polyA = "A" x 20;
		my $i = 1;
		foreach my $line(@file) {
			my ($name, $seq) = split /\n/, $line, 2;
			$seq = uc $seq;
			$seq =~ s/\n//g;
			if ($seq =~ "N") { 
				next;
			} else {
				my $seqG = join ("",$seq,$polyA);
				my $NAME = join("","Cluster",$i); 
				print $OUT1 ">$NAME\n$seq\n";
				print $OUT2 "$seqG";
				push @seqs, $seqG;
				$i++;
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
	system ( "mv *.R1.fa ./fastaForRef" );
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $VsearchOUT" );
}

print "Your '$MockRefName' Mock Reference genome was assembled using the following genotypes:\n@MR_taxa_files\n";
print "We recommend the using '$MR_Genome' as a reference genome for mapping your parsed, high quality reads.";
print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
