#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-4.pl
###########################################################################################################
#########################################################
# Requirement 1: PEAR v0.9.6 (Zhang et al., 2014)
# Requirement 2: VSEARCH v2.6.2 (Rognes et al., 2016)
#########################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw/shuffle/;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(natatime);

######################
# The help function
######################
my $help = $ARGV[0];
my ($pear,$vsearch,$dataType,$barcodesID_file,$threads,$clustering,$raw_seq_length,$pear_length,$pvalue,$id,$derep,$minCl,$MockRefName);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.0\n"
	."Step 4: Cluster reads and assemble the Mock Reference\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-4.pl [options]\n\n"
	."Options:\n"
	."-pr: Path to Pear executable file. String. Default: /usr/local/bin/pear\n"
	."-vs: Path to Vsearch executable file. String. Default: /usr/local/bin/vsearch\n"
	."-d: Data type. Either PE (Paired-End) or SE (Single-End). String. Required.\n"
	."-b: BarcodeID file. File. Required.\n"
	."-t: Number of independent threads used. Numeric. Default: 10\n"
	."-cl: Vsearch clustering algorithm. Either consout or centroids. String. Default: consout\n"
	."-rl: Raw GBS read length. Numeric. Default: 150\n"
	."-pl: Minimum length required after merging to retain read. Numeric. Default: 32\n"
	."-p: p-value for PEAR. Numeric. Default: 0.01\n"
	."-id: Nucleotide identity value required for Vsearch read clustering. Numeric. Default: 0.93\n"
	."-db: Activate the Vsearch derreplication process in blocks. Numeric. Default: Population size\n"
	."-min: Minimum length (bp) for a Mock Reference cluster. Numeric. Default: 32\n"
	."-MR: Mock Reference name. String. Default: GSC.MockRef\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameters values
#################################
$pear = '/usr/local/bin/pear';
$vsearch = '/usr/local/bin/vsearch';
$threads = 10;			$clustering = 'consout';
$raw_seq_length = 150;	$pear_length = 32;
$pvalue = 0.01;			$id = 0.93;
$minCl = 32;			$MockRefName = 'GSC.MR';	

GetOptions(
'pr=s' => \$pear,             # string
'vs=s' => \$vsearch,          # string
'd=s' => \$dataType,          # string
'b=s' => \$barcodesID_file,   # file
't=s' => \$threads,           # numeric
'cl=s' => \$clustering,       # string
'rl=s' => \$raw_seq_length,   # numeric
'pl=s' => \$pear_length,      # numeric
'p=s' => \$pvalue,            # numeric
'id=s' => \$id,               # numeric
'db=s' => \$derep,            # numeric
'min=s' => \$minCl,           # numeric
'MR=s' => \$MockRefName,      # string
) or die "$H\n";

#########################
# Starting GBS-SNP-CROP
#########################
print "\n#################################\n# GBS-SNP-CROP, Step 4, v.4.0\n#################################\n";
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
my $dir = "FastaForRef"; 
unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}

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
	
		print $code_OUT "\n\n###Analyzing paired $R1input1 and $R2input2 reads...\n";
		print $code_OUT `$pear -f $R1input1 -r $R2input2 -o $Pear_out -p $pvalue -n $pear_length -j $threads`;
		
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
		
		my @R1reads;
		my @R2reads;

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
		my $seq_stitch = "T" x 20;

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
	
		# 4. VSEARCH 
		# Removing read redundancy: genotype level
		my $no_reps = join(".","$file","AssembledStitched","noreps", "fa");
		system ( "$vsearch -derep_fulllength $out -minseqlength $minCl -sizeout -output $no_reps" );
		system ( "mv $no_reps $out");
		system( "gzip $out");
	}
	print "DONE.\n\n";
	
	# 4. VSEARCH
	# Source fasta files are joined and dereplicated (but taking notes of numerosity). 
	# Dereplication is done (optionally) at blocks of files to save RAM (see $derep_blocks).
	# Preparing the files and perform a population level dereplication
	my $VsearchIN = join(".","VsearchIN","fa");
	system("touch $VsearchIN");
	my $tmp = join(".","tmp","fa");

	# catting $derep_blocks files (plus previous results) and dereplicating
	my $it = natatime($derep, @MR_taxa_files);
	while(my @files = $it->()){		
		my $fasta_files = join('.AssembledStitched.fa.gz ', @files).'.AssembledStitched.fa.gz ';		
		system("cp $VsearchIN $tmp");
		system("zcat $fasta_files >> $tmp");		
		system ( "$vsearch -derep_fulllength $tmp -sizein -sizeout -minseqlength $minCl -output $VsearchIN");
		system ("mv $fasta_files ./FastaForRef/");
	}

	# Finding population-level initial clusters (centroids)
	system ( "$vsearch -cluster_fast $VsearchIN -sizein -sizeout -id $id -threads $threads -$clustering VsearchOUT.fa");

	# 5. Defining THE Mock reference, masked file and Delete N-containing centroids in order to create the Mock Reference
	$MR_Cluster = join (".","$MockRefName","Clusters","fa");
	$MR_Genome = join (".","$MockRefName","Genome","fa");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";
	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();

	open my $IN6, "<", "VsearchOUT.fa" or die "Unable to open VsearchOUT.fa file\n";
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
			} 
			if (length $seq < $minCl) { # check the minimum length allowed for a cluster
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
				push @bad_coords, $mask_coord;
			}
		}
		while ( $mr =~ m/(T{18,})/g ) { 
			my $end_pos = pos ($mr);
			my $length = length ($1);

			for (my $k = 1; $k <= $length; $k++) {
				my $mask_coord = $end_pos - $length + $k;
				push @bad_coords, $mask_coord;	
			}
		}
		my @uniques = uniq(@bad_coords);
		foreach my $el (@uniques) {
			print $OUT3 "$el\n";
		}
	}
	close $IN6;
	system ( "rm tmp.fa VsearchIN.fa VsearchOUT.fa" );

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

			# fastq -> fasta (no need for newlines, since we did not chomp)
			$R1reads[0] =~ s/@/>/g;
			$R1reads[0] =~ s/ /:/g;
			print $OUT "$R1reads[0]$R1reads[1]";
		}
		close $IN;
		close $OUT;
	}
	
	# 2. VSEARCH 
	# Removing read redundancy: genotype level
	foreach my $file (@MR_taxa_files) {
		my $in = join (".", "$file","R1","fa");
		my $no_reps = join (".", "$file","R1", "noreps","fa");
		system ( "$vsearch -derep_fulllength $in -sizeout -minseqlength $minCl -output $no_reps");
		system ( "mv $no_reps $in");
		system ( "gzip $in");
	}
	
	# Source fasta files are joined and dereplicated (but taking notes of numerosity). 
	# Dereplication is done (optionally) at blocks of files to save RAM (see $derep_blocks).
	my $VsearchIN = join(".","VsearchIN","fa");
	system ("touch $VsearchIN");
	my $tmp = join(".","tmp","fa");
	
	# catting $derep files (plus previous results) and dereplicating
	my $it = natatime($derep, @MR_taxa_files);
	while(my @files = $it->()){	
		my $fasta_files = join('.R1.fa.gz ', @files).'.R1.fa.gz ';
		system("cp $VsearchIN $tmp");
		system("zcat $fasta_files >> $tmp");
		# dereplication
		system ( "$vsearch -derep_fulllength $tmp -sizein -sizeout -minseqlength $minCl -output $VsearchIN");
		system ("mv $fasta_files ./FastaForRef/");
	}	
	
	# Finding population-level initial clusters (centroids)
	system ( "$vsearch -cluster_fast $VsearchIN -sizein -sizeout -id $id -threads $threads -$clustering VsearchOUT.fa");

	# 3. Defining The Mock Ref, masked file and Delete N-containing centroids in order to create the Mock Reference
	$MR_Cluster = join (".","$MockRefName","Clusters","fa");
	$MR_Genome = join (".","$MockRefName","Genome","fa");
	my $masked_pos = "PosToMask.txt";

	open my $OUT1, ">", "$MR_Cluster" or die "Can't load file";
	open my $OUT2, ">", "$MR_Genome" or die "Can't load file";
	open my $OUT3, ">", "$masked_pos" or die "Can't load file";
	print $OUT2 ">MockRefGenome\n";

	my @seqs = ();
	open my $IN4, "<", "VsearchOUT.fa" or die "Unable to open VsearchOUT.fa file\n";
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
			} 
			if (length $seq < $minCl) { # check the minimum length allowed for a cluster
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
				push @bad_coords, $mask_coord;
			}
		}
		while ( $mr =~ m/(T{18,})/g ) { 
			my $end_pos = pos ($mr);
			my $length = length ($1);

			for (my $k = 1; $k <= $length; $k++) {
				my $mask_coord = $end_pos - $length + $k;
				push @bad_coords, $mask_coord;	
			}
		}
		my @uniques = uniq(@bad_coords);
		foreach my $el (@uniques) {
			print $OUT3 "$el\n";
		}
	}
	close $IN4;
	system ( "rm $VsearchIN VsearchOUT.fa $tmp" );	
}

print "\nYour '$MockRefName' Mock Reference genome was assembled using the following genotypes:\n@MR_taxa_files\n";
print "We recommend the using '$MR_Genome' as a reference genome for mapping your parsed, high quality reads.";
print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
