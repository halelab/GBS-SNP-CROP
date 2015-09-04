#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 4. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

##########################################################################################
# Requirement 1: PEAR bioinformatics tool (Zhang et al., 2014)
# Requirement 2: USEARCH bioinformatics tool (Edgar, 2010)
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-4.pl -b <barcodesID file> -rl <Raw GBS read lengths (numeric).  Ex: 100 bp, 150 bp> 
-pl <minimum length required after merging to retain read > -p <p-value for PEAR (Zhang et al., 2014)> 
-id <nucleotide identity value required for USEARCH (Edgar, 2010) read clustering> -t <number of threads dedicated to USEARCH clustering>
-MR <Mock reference name>.\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($barcodesID_file,$raw_seq_length,$pear_length,$pvalue,$id,$threads,$MockRefName);

GetOptions(
'b=s' => \$barcodesID_file,   # file
'rl=s' => \$raw_seq_length,   # numeric
'pl=s' => \$pear_length,	  # numeric
'p=s' => \$pvalue, 			  # numeric
'id=s' => \$id,				  # numeric	 
't=s' => \$threads,		      # numeric
'MR=s' => \$MockRefName,      # string
) or die "$Usage\n$Manual\n";

##########################################################################################
# 1. Use PEAR to merge the parsed R1 and R2 reads to create single reads, if possible
##########################################################################################

my @MR_taxa_files = ();

open IN, "$barcodesID_file" or die "Can't find barcode_ID file\n";
	
while(<IN>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	my $MR_geno = $barcode[2];

	if ($MR_geno eq "YES") {
		push @MR_taxa_files, $TaxaNames;
	}
}

chomp (@MR_taxa_files);

foreach my $file (@MR_taxa_files) {
	my $R1input1 = join (".", "$file","R1","fastq");
	my $R2input2 = join (".", "$file","R2","fastq");
	my $out1 = join("","$file");
		
	print "\n\nAssembling paired $R1input1 and $R2input2 reads into single reads using PEAR...\n";
	system ( "pear -f $R1input1 -r $R2input2 -o $out1 -p $pvalue -n $pear_length -j $threads" );
}

close IN;

system ( "rm *.discarded.fastq" );

print "\nSome paired-end reads were successfully assembled into singletons.\n\n";

##########################################################################################
# 2. Stitch unassembled R1 and R2 reads together with an intermediate run of 20 high-quality A's
##########################################################################################

print "Manually stitching together the unassembled reads...\n\n";

foreach my $file (@MR_taxa_files) {
	my $R1file2 = join (".", "$file","unassembled","forward","fastq");
	my $R2file2 = join (".", "$file","unassembled","reverse","fastq");

	my @R1read;
	my @R1reads;
	
	open IN, "$R1file2" or die "Can't open $R1file2: $!\n";
	
	my $i = 1;
	while(<IN>) {
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
	
	my @R2read;
	my @R2reads;
	
	open IN, "$R2file2" or die "Can't open $R2file2: $!\n";
	
	my $j = 1;
	while(<IN>) {
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
	
	my $size = scalar @R1reads - 1;
	
	my $stitched_file = join (".", "$file","stitched","fastq");
		
	open stitched_OUT, ">", "$stitched_file";
		
	my $stitched_tally = 0;
	my $unstitched_tally = 0;
		
	for (my $k = 0; $k <= $size; $k++) {		
		my $R1_length = length($R1reads[$k][1]);
		my $R2_length = length($R2reads[$k][1]);		
		if ( $R1_length < ($raw_seq_length - 19) || $R2_length < ($raw_seq_length - 5) ) {
			$unstitched_tally++;
			next;
		} else {
			my $assembled_seq = join ("", "$R1reads[$k][1]", "AAAAAAAAAAAAAAAAAAAA", "$R2reads[$k][1]");
			my $assembled_qual = join ("", "$R1reads[$k][3]", "IIIIIIIIIIIIIIIIIIII", "$R2reads[$k][3]");
			print stitched_OUT "$R1reads[$k][0]\n";
			print stitched_OUT "$assembled_seq\n";
			print stitched_OUT "+\n";
			print stitched_OUT "$assembled_qual\n";
			$stitched_tally++;
			next;
		}
	}	
		
	close stitched_OUT;
			
	my $unstitched_percentage = ( $unstitched_tally / ( $stitched_tally + $unstitched_tally ) ) * 100;
	
  	print "For unassembled read pair $R1file2 and $R2file2:\n";	
	print "Total number of stitched read pairs = $stitched_tally\n";
	print "Total number of unstitchable read pairs = $unstitched_tally\n";
	print "Percent of read pairs that were unstitchable = $unstitched_percentage\n";
	print "DONE.\n\n\n";
	
}

close IN;

system ( "rm *.unassembled.*" );

##########################################################################################
# 3. Concatenate the merged and stitched reads for each genotype into a single file 
# for use in building the Mock Reference 
##########################################################################################

foreach my $file (@MR_taxa_files) {
	my $assembled = join (".", "$file","assembled","fastq");
	my $stitched = join (".", "$file","stitched","fastq");
	my $out = join(".","$file","AssembledStitched","fastq");
	
	print "Concatenating $assembled and $stitched files...";
	system ( "cat $assembled $stitched > $out" );
	print "DONE.\n";
}

#########################################################################################
# 4. Use USEARCH to cluster reads within each genotype 
#########################################################################################

# Step 4A: Sort the FASTQ sequences in order of decreasing length (bp)
print "\n\nSort the FASTQ sequences in order of decreasing length...\n";
foreach my $file (@MR_taxa_files) {
	my $in4A = join (".", "$file","AssembledStitched","fastq");
	my $out4A = join(".","$file","sorted_by_length","fasta");
	print "Processing file $in4A...\n";
	system ( "usearch -sortbylength $in4A -fastaout $out4A" );
	print "DONE.\n";
}

# Step 4B: Find initial clusters (centroids)
print "\nFind genotype-specific initial clusters (centroids)...\n";
foreach my $file (@MR_taxa_files) {
	my $in4B = join(".","$file","sorted_by_length","fasta");
	my $out4B = join(".","$file","clusters","fasta");
	print "Processing file $in4B...\n";
	system ( "usearch -cluster_fast $in4B -id $id -threads $threads -consout $out4B -sizeout" );
	print "DONE.\n";
}

# Step 4C: Sort reads by cluster size (# of reads in each cluster)
print "\nSort clusters (centroids) by depth...\n";
foreach my $file (@MR_taxa_files) {
	my $in4C = join(".","$file","clusters","fasta");
	my $out4C = join(".","$file","sorted_by_size","fasta");
	print "Processing file $in4C...\n";
	system ( "usearch -sortbysize $in4C -fastaout $out4C" );
	print "DONE.\n";
}

# Step 4D: Recluster the centroids
print "\nRecluster the genotype-specific centroids...\n";
# my $out4D = " ";
foreach my $file (@MR_taxa_files) {
	my $in4D = join(".","$file","sorted_by_size","fasta");
	my $out4D = join(".","$file","reclusters","fasta");
	print "Processing file $in4D...\n";
	system ( "usearch -cluster_fast $in4D -id $id -threads $threads -consout $out4D" );
	print "DONE.";
}

print "\n\nAll sub-steps for clustering reads within genotype(s) were completed!\n";

# Define Mock reference

my $MR_Cluster = join (".","$MockRefName","MockRef_Clusters","fasta");
my $MR_Genome = join (".","$MockRefName","MockRef_Genome","fasta");
open (OUT1, ">$MR_Cluster") || die "cant load file";
open (OUT2, ">$MR_Genome") || die "cant load file";

print OUT2 ">MockRefGenome\n";

#########################################################################################
# 5. Use USEARCH to cluster genotype-specific centroids across the population 
#########################################################################################

if (scalar (@MR_taxa_files) > 1) {

	print "\nConcatenating the genotype-specific non-redundant centroids...";
	my $in5A = join(".",$MockRefName,"consensus","fasta");
	system ( "cat *reclusters.fasta > $in5A" );
	print "DONE.\n";

	# Step 5A: Sort the FASTA sequences in order of decreasing length (bp)
	print "\nSort the non-redundant centroids in order of decreasing length...\n";
	my $out5A = join (".",$MockRefName,"sorted_by_length","fasta");
	system ( "usearch -sortbylength $in5A -fastaout $out5A" );
	print "DONE.\n";

	# Step 5B: Find initial centroids in population level
	print "\nFind population-level clusters (centroids)...\n";
	my $out5B = join (".",$MockRefName,"clusters","fasta");
	system ( "usearch -cluster_fast $out5A -id $id -threads $threads -consout $out5B -sizeout");
	print "DONE.\n";

	# Step 5C: Sort reads by cluster size (# of reads in each cluster)
	print "\nSort population-level clusters (centroids) by depth...\n";
	my $out5C = join (".",$MockRefName,"sorted_by_size","fasta");
	system ( "usearch -sortbysize $out5B -fastaout $out5C" );
	print "DONE.\n";

	# Step 5D: Recluster the centroids
	print "\nRecluster the population-level centroids...\n";
	my $out5D = join (".",$MockRefName,"reclusters","fasta");
	system ( "usearch -cluster_fast $out5C -id $id -threads $threads -consout $out5D" );
	print "DONE.\n";

	print "All sub-steps for clustering population-level centroids were completed!\n";

##########################################################################################
# 6. Delete N-containing centroids and create the Mock Reference
##########################################################################################

	open IN5, "$out5D" or die "Unable to open $out5D file\n";
	{						
		local $/ = ">";		
		my @file = <IN5>;
		shift @file;		
		chomp @file;
		foreach my $line(@file) {
			my ($name, $seq) = split /\n/, $line, 2; 	
			$seq =~ s/\n//g; 							
			if ($seq =~ "N") { 
				next;
			} else {
				print OUT1 ">$name\n$seq\n";
				print OUT2 "$seq";
			}
		}
	}
	close IN5;

} elsif (scalar (@MR_taxa_files) == 1) { 

	foreach my $file (@MR_taxa_files) {
		my $out4D = join(".","$file","reclusters","fasta");
		open IN4, "$out4D" or die "Unable to open $out4D file\n";
		{						
			local $/ = ">";		
			my @file = <IN4>;
			shift @file;		
			chomp @file;
			foreach my $line(@file) {
				my ($name, $seq) = split /\n/, $line, 2; 	
				$seq =~ s/\n//g; 							
				if ($seq =~ "N") { 
					next;
				} else {
					print OUT1 ">$name\n$seq\n";
					print OUT2 "$seq";
				}
			}
		}
		close IN4;
	} 
}
	
sub main {
   	my $dir = "fastqForRef"; 
   	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
}
main();

system ( "mv *.AssembledStitched.fastq ./fastqForRef" );
system ( "rm *.sorted_by_length.fasta" );
system ( "rm *.clusters.fasta" );
system ( "rm *.sorted_by_size.fasta" );
system ( "rm *.reclusters.fasta" );
system ( "rm *.consensus.fasta" );
system ( "rm *.assembled.fastq" );
system ( "rm *.stitched.fastq" );

print "\nGood job! Your '$MockRefName' Mock Reference genome was assembled using follow genotypes:\n@MR_taxa_files\n";

print "We recommend the use of '$MR_Genome' as a referenced genome for mapping your parsed, high quality paired-end reads.\n";

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;