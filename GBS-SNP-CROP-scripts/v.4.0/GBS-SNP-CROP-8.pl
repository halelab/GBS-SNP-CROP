#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-8.pl
###########################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

######################
# The help function
######################
my $help = $ARGV[0];
my ($SNPmatrix,$output,$barcodesID_file,$tools);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.0\n"
	."Step 8: Creating input files for downstream analyses\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-8.pl [options]\n\n"
	."Options:\n"
	."-in: Final genotyping matrix input file. The output from step 7. File. Default: GSC.GenoMatrix.txt\n"
	."-out: Output label name without extension. String. Default: GSC\n"
	."-b: BarcodeID file. File. Required.\n"
	."-formats: The name(s) (R, Tassel, Plink, VCF, HetFreq) for which the GBS-SNP-CROP final genotyping matrix should be converted. If more than one format is desired, the names should be separated by commas without any space. String. Default: R,P,T,V,H\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}
#################################
# Setting the parameters values
#################################
$SNPmatrix = 'GSC.GenoMatrix.txt';	$output = 'GSC';
$tools =  'R,P,T,V,H';

GetOptions(
'in=s' => \$SNPmatrix,          # file
'out=s' => \$output,            # string
'b=s' => \$barcodesID_file,     # file
'formats=s' => \$tools,         # file (R, Tassel, Plink, VCF, HetFreq)
) or die "H\n";

#########################
# Starting GBS-SNP-CROP
#########################
print "\n#################################\n# GBS-SNP-CROP, Step 8, v.4.0\n#################################\n";
my $sttime = time;

# Defining outputs
my $R_out = join (".","$output","R","txt");
my $Tassel_out = join (".","$output","Tassel","hmp","txt");
my $Plink_out = join (".","$output","Plink","tped");
my $vcf_out = join (".","$output","vcf");
my $hetfreq_out = join (".","$output","HetFreq","txt");

###########################
# Creating R input file
###########################
if ($tools =~ "R" or $tools =~ "r") {
	print "\nConverting the genotyping matrix into an input file compatible with R software...";

	open my $IN1, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_R, ">", "$R_out" or die "Can't initialize $R_out output file";
	
	while(<$IN1>) {
		my @snp = split "\t", $_;
		chomp @snp;

		my $col = join ("\t",$snp[0],$snp[1]);
		my $primary = $snp[3];
		my $secondary = $snp[4];

		my @out = ();
		push @out, "$col";
	
		my $length = scalar(@snp) - 1;
	
		for ( my $i=12; $i<=$length; $i++ ) {
			my @geno1 = split /\|/, $snp[$i];
			my $geno2 = $geno1[0];

			if ($geno2 eq "-") {
				push @out, "NA";
				next;
			} elsif ($geno2 eq join("/",$primary,$primary)) {
				push @out, "0";
				next;
			} elsif ($geno2 eq join("/",$primary,$secondary)) {
				push @out, "0.5";
				next;
			} elsif ($geno2 eq join("/",$secondary,$primary)) {
				push @out, "0.5";
				next;
			} elsif ($geno2 eq join("/",$secondary,$secondary)) {
				push @out, "1";
				next;
			} else {
				push @out, "NA";
				next;
			}
		}
		my $line = join "\t", @out;
		print $OUT_R "$line\n";
	}

	close $IN1;
	close $OUT_R;
	print "\nDONE.\n";
	} else {
}
	
#####################################
# Creating Tassel HapMap input file
#####################################
if ($tools =~ "T" or $tools =~ "t") {
	print "\nConverting the genotyping matrix into an input file compatible with TASSEL GUI software...";

	my @TaxaNames = ();
	open my $BAR, "<", "$barcodesID_file" or die "Can't find $barcodesID_file file\n";
	while(<$BAR>) {
		my $barcodesID = $_;
		chomp $barcodesID;
		my @barcode = split("\t", $barcodesID);
		push @TaxaNames, $barcode[1];
	}
	s/_/./g for @TaxaNames;

	my @output;

	open my $IN2, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_T, ">", "$Tassel_out" or die "Can't initialize $Tassel_out output file";

	while(<$IN2>) {
		my @snp = split "\t", $_;
		chomp @snp;
		if ($snp[2] eq "Indel") {
			next;
		}
		my $col1 = join ("_",$snp[0],$snp[1]);
		my $alleles = join ("/",$snp[3],$snp[4]);

		my @out;

		my $length = scalar(@snp) - 1;

		for ( my $i=12; $i<=$length; $i++ ) {
			my @geno1 = split /\|/, $snp[$i];
			my $geno2 = $geno1[0];

			if ($geno2 eq "-") {
				push @out,"N";
				next;
			} elsif ($geno2 eq join("/","A","A")) {
				push @out,"A";
				next;
			} elsif ($geno2 eq join("/","C","C")) {
				push @out,"C";
				next;
			} elsif ($geno2 eq join("/","G","G")) {
				push @out,"G";
				next;
			} elsif ($geno2 eq join("/","T","T")) {
				push @out,"T";
				next;
			} elsif ( $geno2 eq join("/","A","C") or $geno2 eq join("/","C","A") ) {
				push @out,"M";
				next;
			} elsif ( $geno2 eq join("/","A","G") or $geno2 eq join("/","G","A") ) {
				push @out,"R";
				next;
			} elsif ( $geno2 eq join("/","A","T") or $geno2 eq join("/","T","A") ) {
				push @out,"W";
				next;
			} elsif ( $geno2 eq join("/","C","G") or $geno2 eq join("/","G","C") ) {
				push @out,"S";
				next;
			} elsif ( $geno2 eq join("/","C","T") or $geno2 eq join("/","T","C") ) {
				push @out,"Y";
				next;
			} elsif ( $geno2 eq join("/","G","T") or $geno2 eq join("/","T","G") ) {
				push @out,"K";
				next;
			} else {
				push @out,"N";
				next;
			}
		}
		my $line = join ("\t","$col1","$alleles","$snp[0]","$snp[1]","+","NA","GBS-SNP-CROP","GBS","RefName","Custom","QC+","@out");
		push @output,"$line\n";

	}
	print $OUT_T "#rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t@TaxaNames\n@output";

	close $IN2;
	close $OUT_T;
	print "\nDONE.\n";
	} else {
}
	
############################
# Creating PLINK input file 
############################
if ($tools =~ "P" or $tools =~ "p") {
	print "\nConverting the genotyping matrix into an input file compatible with Plink software...";

	open my $IN3, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_P, ">", "$Plink_out" or die "Can't initialize $Plink_out output file";

	while (<$IN3>) {
		chomp;
		my @input = split ("\t", $_);
		if ($input[2] eq "Indel") {
			next;
		}
		my $header = $input[0];
		my $position = $input[1];
		my $snp_identifier = join ("_",$header,$position);
		$header =~ s/chr0//;
		$header =~ s/chr//;
		my $Tpedfile = join("\t",$header,$snp_identifier,"0",$position);

		my @geno = splice @input, 12; 
		foreach (@geno) {
			my @geno1 = split /\|/, $_;
			s/-/0 0/ for @geno1;
			s/\// / for @geno1;
			$Tpedfile = join("\t",$Tpedfile,$geno1[0]);
		}
		print $OUT_P "$Tpedfile\n";
	}

	close $IN3;
	close $OUT_P;
	print "\nDONE.\n";
	} else {
}

#############################
# Creating VCF output format 
#############################
if ($tools =~ "V" or $tools =~ "v") {
	print "\nConverting the genotyping matrix into a VCF file format...";
	
	open my $IN4, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_V, ">", "$vcf_out" or die "Can't initialize $vcf_out output file";
	
	my @TaxaNames = ();
	open my $BAR, "<", "$barcodesID_file" or die "Can't find $barcodesID_file file\n";
	while(<$BAR>) {
		my $barcodesID = $_;
		chomp $barcodesID;
		my @barcode = split("\t", $barcodesID);
		push @TaxaNames, $barcode[1];
	}
	s/_/./g for @TaxaNames;
	my $Names = join "\t", @TaxaNames;
	
	# Printing output VCF header
	my $datestring = localtime();
	
	print $OUT_V "##fileformat=VCFv4.2\n"
	."##fileDate=$datestring\n"
	."##source=GBS-SNP-CROP\n"
	."##phasing=partial\n"
	."##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n"
	."##INFO=<ID=AF,Number=A,Type=Integer,Description=\"Allele Frequency\">\n"
	."##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
	."##INFO=<ID=AV,Number=1,Type=Integer,Description=\"Average Depth\">\n"
	."##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
	."##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	."##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n"
	."#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$Names\n";
	
	while (<$IN4>) {
		chomp;
		my @input = split ("\t", $_);
		my $identifier = $input[0];
		my $pos = $input[1];
		my $ref = $input[3];
		my $alt = $input[4];
		my $homo1 = $input[7];
		my $hetero = $input[8];
		my $homo2 = $input[9];
		
		# Alternative allele
		my $ref_p;
		my $alt_p;
		if ($input[2] eq "Indel") {
			if ($alt =~ "[+]") {
				$alt =~ s/([+]\d)//g; 
				$alt_p = "$ref"."$alt";
				$ref_p = $ref;
			} elsif ($alt =~ "[-]") {
				$alt =~ s/([-]\d)//g; 
				$ref_p = "$ref"."$alt";
				$alt_p = $ref;
			}
		} elsif ($input[2] eq "SNP") {
			$ref_p = $ref;
			$alt_p = $alt;
		}
		$alt_p =~ s/([+-]\d)//g;
		$alt =~ s/([+-]\d)//g;
		$ref_p =~ s/([+-]\d)//g;
		$ref =~ s/([+-]\d)//g;
			
		# Scored genotypes (NS)
		my $NS = $homo1 + $hetero + $homo2;
		# AV
		my $AV = sprintf("%.2f",$input[5]);


		# Allele Count (AC), Allele Freq (AF) and Depth (DP)
		my @alleles = ();
		my @depths = ();
		my @depths_alt = ();		

		for ( my $i=12; $i<=scalar(@input) - 1; $i++ ) {
			my @genos = split /\|/, $input[$i];
			
			my $x = $genos[0];
			$x =~ s/-/N\/N/g;
			my @x1 = split /\//, $x;
			$x1[0] =~ s/([+-]\d)//g;
			$x1[1] =~ s/([+-]\d)//g;
			push @alleles, $x1[0], $x1[1];

			my $y = $genos[1];
			$y =~ s/-/0\/0/g;
			my @y1 = split /\//, $y;
			push @depths, $y1[0], $y1[1];
			push @depths_alt, $y1[1];	
		}
		my $DP;
		$DP += $_ for @depths;
		
		# AC and AF
		my $Alt_cnt = 0;
		my $Ref_cnt = 0;
		
		foreach my $allele (@alleles) {
			#$allele =~ s/([+-]\d)//g;
			if ($allele eq $alt){
				$Alt_cnt++;
			} else {
				$Ref_cnt++;
			}
		}
		my $AC = $Alt_cnt;
		my $AF = sprintf("%.3f",($AC / ($AC + $Ref_cnt)));
		
		# Genotypes
		my @genos = ();
		for ( my $i=12; $i<=scalar(@input) - 1; $i++ ) {
			my @k = split /\|/, $input[$i];
			my @l = split /\//, $k[1];
			if ($k[0] eq "-") {
				$k[0] = "./.";
				$k[1] = ".,.";
				push @genos, join (":", $k[0],$k[1]);
				next;
			} else {
				$k[1] =~ s/\//,/g;
				$k[0] =~ s/([+-]\d)//g;
				$k[0] =~ s/$alt/1/g;
				$k[0] =~ s/$ref/0/g;
				push @genos, join (":", $k[0],$k[1]);
				next;
			}
		}
		my $Geno_Names = join "\t", @genos;
		my $line = join ("\t","$identifier","$pos",".","$ref_p","$alt_p",".","PASS","AC=$AC;AF=$AF;DP=$DP;AV=$AV;NS=$NS","GT:AD","$Geno_Names");
		print $OUT_V "$line\n";
	}

	close $IN4;
	close $OUT_V;
	print "\nDONE.\n";
	} else {
}

####################################################
# Estimating the depth ratio for heterozygotes loci 
####################################################
if ($tools =~ "H" or $tools =~ "h") {
	print "\nEstimating the depth ratio for all heterozygotes loci...";
	
	open my $IN5, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_H, ">", "$hetfreq_out" or die "Can't initialize $R_out output file";
	
	while(<$IN5>) {	
		my @snp = split "\t", $_;
		chomp @snp;

		my @out1 = ();
		my @out2 = ();

		my $col = join ("\t",$snp[0],$snp[1]);
		push @out1, "$col";
		push @out2, "$col";
	
		for ( my $i=12; $i<=scalar(@snp) - 1; $i++ ) {
			my @geno = split /\|/, $snp[$i];

			if ( ($geno[0] eq "-") or ($geno[0] eq "A/A") or ($geno[0] eq "C/C") or ($geno[0] eq "G/G") or ($geno[0] eq "T/T") ) {
				push @out1, "NA";
				push @out2, "NA";
				next;
		
			} else {
				my @depth = split /\//, $geno[1];
				if ($depth[0] >= $depth[1]) {
					my $freq1 = ($depth[1]/($depth[0] + $depth[1])); # lower / sum
					my $freq2 = ($depth[0]/($depth[0] + $depth[1])); # higher / sum
					push @out1, $freq1;
					push @out2, $freq2;
					next;
				} elsif ($depth[0] < $depth[1])  {
					my $freq1 = ($depth[0]/($depth[0] + $depth[1])); # lower / sum
					my $freq2 = ($depth[1]/($depth[0] + $depth[1])); # higher / sum
					push @out1, $freq1;
					push @out2, $freq2;
					next;
				}
			}
		}
		my $line1 = join "\t", @out1;
		my $line2 = join "\t", @out2;
		print $OUT_H "$line1\n$line2\n";
	}
	close $IN5;
	close $OUT_H;
	print "\nDONE.\n";
}

print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
