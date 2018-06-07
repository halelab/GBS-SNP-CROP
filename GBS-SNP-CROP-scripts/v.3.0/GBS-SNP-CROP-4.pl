#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 4. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y.
#########################################################################################################################

#######################################################
# Requirement 1: PEAR v0.9.6 (Zhang et al., 2014)
# Requirement 2: VSEARCH v2.6.2 (Rognes et al., 2016)
#######################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw/shuffle/;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(natatime);

my $Usage = "Usage: perl GBS-SNP-CROP-4.pl -pr <PEAR path w/ exe file> -vs <Vsearch path w/ exe file> -d <data type, PE = Paired-End or SR = Single-End> -b <barcode-ID file> -rl <Raw GBS read lengths> -pl <minimum length required after merging to retain read>\n"
."-p <p-value for PEAR> -id <nucleotide identity value required for VSEARCH reads clustering> -t <number of threads dedicated to VSEARCH> -MR <Mock reference name>\n"
."-db <optional size of dereplication block> -rs <this flag activate remove singletons centroids>\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($pear,$vsearch,$dataType,$barcodesID_file,$raw_seq_length,$pear_length,$pvalue,$id,$threads,$MockRefName,$derep,$no_single);
$pear = '/usr/local/bin/pear';
$vsearch = '/usr/local/bin/vsearch';

GetOptions(
'pr=s' => \$pear, 	      # string - PEAR path w/ exe file
'vs=s' => \$vsearch,          # string - Vsearch path w/ exe file
'd=s' => \$dataType,          # string - "PE" or "SE"
'b=s' => \$barcodesID_file,   # file
'rl=s' => \$raw_seq_length,   # numeric
'pl=s' => \$pear_length,      # numeric
'p=s' => \$pvalue,            # numeric
'id=s' => \$id,               # numeric
't=s' => \$threads,           # numeric
'MR=s' => \$MockRefName,      # string
'db=i' => \$derep,	      # numeric
'rs' => \$no_single,	      # optional flag - it turn on the singletons centroids deletion
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 4, v.3.1\n#################################\n";
my $sttime = time;

# Reading the barcode ID file and set genotypes to build the MR
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

# If not specified, the default value for dereplication blocks is the number of genotypes set to build the MR
if (! defined $derep){
	$derep = scalar @MR_taxa_files;
}
# Creating a directory
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

		# Assuming well formatted fastq files, we can read four lines at a time. No need for chomping or storing the whole file in memory.
		while(!eof($AIN)) {
			my @R1read;
			$R1read[0] = readline($AIN); #fastq @header
			$R1read[1] = readline($AIN); #bases
			readline($AIN);              #+ (ignored)
			readline($AIN);              #quality (ignored)
			
			$R1read[0] =~ s/@/>/g;
			$R1read[0] =~ s/ /:/g;
			print $AOUT "$R1read[0]$R1read[1]";
		}
		unlink $Ain;
		close $AIN;
		close $AOUT;
		system ( "rm *.discarded.fastq" );
		
		# 2. Stitch unassembled R1 and R2 reads together with an intermediate run of 20 high-quality A's
		print $code_OUT "\n\n### Manually stitching together unassembled reads results:\n";

		my $R1file2 = join (".", "$file","unassembled","forward","fastq");
		my $R2file2 = join (".", "$file","unassembled","reverse","fastq");
		open my $IN1, "<", "$R1file2" or die "Can't open $R1file2: $!\n";
		open my $IN2, "<", "$R2file2" or die "Can't open $R2file2: $!\n";

		my @R1read;
		my @R1reads;
		
		while(!eof($IN1)) {
			my @R1read;
			$R1read[0] = readline($IN1); # fastq @header
			$R1read[1] = readline($IN1); # bases
			readline($IN1);              # + (ignored)
			readline($IN1);              # quality (ignored)

			chomp(@R1read);
			push @R1reads, [ @R1read ];
		}
		while(!eof($IN2)) {
			my @R2read;
			$R2read[0] = readline($IN2); #fastq @header
			$R2read[1] = readline($IN2); #bases
			readline($IN2);              #+ (ignored)
			readline($IN2);              #quality (ignored)

			chomp(@R2read);
			push @R2reads, [ @R2read ];
		}
		close $IN1;
		close $IN2;

		my $stitched_file = join (".", "$file","stitched","fasta");
		open my $stitched_OUT, ">", "$stitched_file" or die "Can't open $stitched_file: $!\n";

		my $stitched_tally = 0;
		my $unstitched_tally = 0;
		my $seq_stitch = "A" x 20;

		for (my $k = 0; $k <= scalar @R1reads - 1; $k++) {
			my $R1_length = length($R1reads[$k][1]);
			my $R2_length = length($R2reads[$k][1]);
			if ( $R1_length < ($raw_seq_length - 19) || $R2_length < ($raw_seq_length - 5) ) {
				$unstitched_tally++;
			} else {
				my $assembled_seq = join ("", "$R1reads[$k][1]", "$seq_stitch", "$R2reads[$k][1]");
				$R1reads[$k][0] =~ s/@/>/g;
				$R1reads[$k][0] =~ s/ /:/g;
				print $stitched_OUT "$R1reads[$k][0]\n";
				print $stitched_OUT "$assembled_seq\n";
				$stitched_tally++;
			}
		}
		close $stitched_OUT;

		my $unstitched_percentage = ( $unstitched_tally / ( $stitched_tally + $unstitched_tally ) ) * 100;
		print $code_OUT "\nFor paired of unassembled files $R1file2 and $R2file2:\n";
		print $code_OUT "Total number of stitched read pairs = $stitched_tally\n";
		print $code_OUT "Total number of unstitchable read pairs = $unstitched_tally\n";
		print $code_OUT "Percent of unassembled read pairs that were unstitchable = $unstitched_percentage\n\n";
		
		# Removing the unassembled files
		unlink $R1file2;
		unlink $R2file2;
		
		# 3. For each genotype, concatenate the merged and stitched reads into a single file for use in Vsearch (build the Mock Ref)
		# Also, remove redundancy among centroids entries (Nelson's dereplication).
		my $assembled = join (".", "$file","assembled","fasta");
		my $stitched = join (".", "$file","stitched","fasta");
		my $out = join(".","$file","AssembledStitched","fa");
		system ( "cat $assembled $stitched > $out" );
		unlink $assembled;
		unlink $stitched;
		
		# Dereplication process
		my $fasta_no_reps = join(".","$file","AssembledStitched","noreps", "fa");
		print $code_OUT `vsearch -derep_fulllength $out -sizeout -output $fasta_no_reps 2>&1`;
		# Replacing fasta with reps with new dereplicated file and compressing the result file
		system ( "mv $fasta_no_reps $out");
		system( "gzip $out");
	}
	print "DONE.\n";

	# 4. Use VSEARCH to cluster reads
	# Preparing the files and perform a population level dereplication
	my $VsearchIN = join(".","VsearchIN","fa");
	system("touch $VsearchIN");
	my $tmp = join(".","tmp","fa");

	# catting $derep plus previous results and dereplicating
	my $it = natatime($derep, @MR_taxa_files);
	while(my @files = $it->()){
		# list of files to be dereplicated in the current batch
		my $fasta_files = join('.AssembledStitched.fa.gz ', @files).'.AssembledStitched.fa.gz ';		
		# catting all together
		system("cp $VsearchIN $tmp");
		system("zcat $fasta_files >> $tmp");		
		# dereplication
		system ( "vsearch -derep_fulllength $tmp -sizein -sizeout -output $VsearchIN");
		# files are no longer necessary, we store them away
		system ("mv $fasta_files ./fastaForRef/");
	}
	
	# Removing singletons
	if ($no_single) {
		system("mv $VsearchIN $tmp");
		system("vsearch -fastx_filter $tmp -minsize 2 -fastaout $VsearchIN");
	}

	# Step 4A: Sorting the full set of genotype-specific centroids in order of decreasing length
	print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
	my $out4A = join (".",$MockRefName,"sorted_by_length","fasta");
	system ( "vsearch -sortbylength $VsearchIN -sizein -sizeout -output $out4A" );
	
	# Step 4B: Finding population-level initial clusters (centroids)
	print "\nFinding population-level initial clusters (centroids)...\n";
	my $out4B = join (".",$MockRefName,"clusters","fasta");
	system ( "vsearch -cluster_fast $out4A -sizein -sizeout -id $id -threads $threads -consout $out4B -sizeout");

	# Step 4C: Sorting population-level initial clusters (centroids) by depth
	print "\nSorting population-level initial clusters (centroids) by depth...\n";
	my $out4C = join (".",$MockRefName,"sorted_by_size","fasta");
	system ( "vsearch -sortbysize $out4B -sizein -sizeout -output $out4C" );
	
	# Step 4D: Reclustering the population-level centroids
	print "\nReclustering the population-level centroids...\n";
	my $VsearchOUT = join (".",$MockRefName,"reclusters","fasta");
	system ( "vsearch -cluster_fast $out4C -sizein -sizeout -id $id -threads $threads -consout $VsearchOUT" );

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
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $VsearchOUT $tmp" );

############################
# Parsing Single-End data 
############################
} elsif ($dataType eq "SE") {
	print "Parsing Single-End reads...\n\n";

	# 1. Transforming FASTQ high quality reads into FASTA reads
	foreach my $file (@MR_taxa_files) {
		my $R1input1 = join (".", "$file","R1","fq","gz");
		my $out = join (".", "$file","R1","fa");
		open my $IN, '-|', 'gzip', '-dc', $R1input1 or die "Can't open file $R1input1: $!\n";
		open my $OUT, ">", "$out" or die "Can't load $out\n";

		my @R1reads;
		while(! eof ($IN)) {
			$R1reads[0] = readline($IN); # fastq @header
			$R1reads[1] = readline($IN); # bases
			readline($IN);               # + (ignored)
			readline($IN);               # quality (ignored)

			# check on minimum read length (it contains the \n at the end, thus we use its length minus one
			if ((length($R1reads[1]) - 1) < $pear_length) {
				next;
			}
			# fastq -> fasta (no need for newlines, since we did not chomp)
			$R1reads[0] =~ s/@/>/g;
			$R1reads[0] =~ s/ /:/g;
			print $OUT "$R1reads[0]$R1reads[1]";
		}
		close $IN;
		close $OUT;
	}

	# 2. Use VSEARCH to to dereplicate each fasta file separatedly, then compressing and to cluster reads 
	# Dereplication
	foreach my $file (@MR_taxa_files) {
		my $in = join (".", "$file","R1","fa");
		my $no_reps = join (".", "$file","R1", "noreps","fa");
		system ( "vsearch -derep_fulllength $in -sizeout -minseqlength $pear_length -output $no_reps");
		system ( "mv $no_reps $in");
		system ( "gzip $in");
	}
	
	# Clustering reads
	my $VsearchIN = join(".","VsearchIN","fa");
	system ( "touch $VsearchIN" );
	my $tmp = join(".","tmp","fa");
	
	# catting $derep files plus previous results and aplly the population level dereplicating
	my $it = natatime($derep, @MR_taxa_files);
	while(my @files = $it->()){		
		my $fasta_files = join('.R1.fa.gz ', @files).".R1.fa.gz";
		system("cp $VsearchIN > $tmp");
		system("zcat $fasta_files >> $tmp");
		# dereplication
		system ( "vsearch -derep_fulllength $tmp -sizein -sizeout -minseqlength $pear_length -output $VsearchIN");
		system ("mv $fasta_files ./fastaForRef/");
	}
	
	# removing singletons
	if ($no_single){
		system("mv $VsearchIN $tmp");
		system("vsearch -fastx_filter $tmp -minsize 2 -fastaout $VsearchIN");
	}

	# Step 2A: Sorting the full set of genotype-specific centroids in order of decreasing length
	print "\n\nSorting the full set of genotype-specific centroids in order of decreasing length...\n";
	my $out2A = join (".",$MockRefName,"sorted_by_length","fasta");
	system ( "vsearch -sortbylength $VsearchIN -sizeout -output -output $out2A" );

	# Step 2B: Finding population-level initial clusters (centroids)
	print "\nFinding population-level initial clusters (centroids)...\n";
	my $out2B = join (".",$MockRefName,"clusters","fasta");
	system ( "vsearch -cluster_fast $out2A -sizeout -output -id $id -threads $threads -consout $out2B -sizeout");

	# Step 2C: Sorting population-level initial clusters (centroids) by depth
	print "\nSorting population-level initial clusters (centroids) by depth...\n";
	my $out2C = join (".",$MockRefName,"sorted_by_size","fasta");
	system ( "vsearch -sortbysize $out2B -sizeout -output -output $out2C" );

	# Step 2D: Reclustering the population-level centroids
	print "\nReclustering the population-level centroids...\n";
	my $VsearchOUT = join (".",$MockRefName,"reclusters","fasta");
	system ( "vsearch -cluster_fast $out2C -sizeout -output -id $id -threads $threads -consout $VsearchOUT" );

	print "\nAll sub-steps for clustering population-level centroids were completed!\n";
	unlink $VsearchIN;

	# 3. Defining The Mock Ref, masked file and Delete N-containing centroids in order to create the Mock Reference
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
	system ( "rm *.sorted_by_length.fasta *.clusters.fasta *.sorted_by_size.fasta $tmp $VsearchOUT" );
}

print "Your '$MockRefName' Mock Reference genome was assembled using the following genotypes:\n@MR_taxa_files\n";
print "We recommend the using '$MR_Genome' as a reference genome for mapping your parsed, high quality reads.";
print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
