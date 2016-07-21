#!/usr/bin/perl

##################################################################################################################
# GBS-SNP-CROP, Step 8. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
##################################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-8.pl -in <input file: SNP genotyping matrix> -b <barcodesID file>\n"
."-formats <The name(s) of the software packages (R, Tassel, Plink) for which an input-compatible file format is desired>.\n";
my $Manual = "Please see the User Manual on the GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($SNPmatrix,$output,$barcodesID_file,$tools);

GetOptions(
'in=s' => \$SNPmatrix,          # file
'out=s' => \$output,            # file
'b=s' => \$barcodesID_file,     # file
'formats=s' => \$tools,         # file
) or die "Error in command line arguments.\n$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 8, v.2.0\n#################################\n";

# Defining outputs
my $R_out = join (".","$output","R_in","txt");
my $Tassel_out = join (".","$output","Tassel_in","hmp","txt");
my $Plink_out = join (".","$output","Plink_in","tped");

###########################
# Creating R input file
###########################
if ($tools =~ "R" or $tools =~ "r"){

	print "\nConverting the SNP genotyping matrix into an input file compatible with R software...";

	open my $IN1, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_R, ">", "$R_out" or die "Can't initialize $R_out output file";
	
	while(<$IN1>) {
		my @snp = split "\t", $_;
		chomp @snp;

		my $col = join ("\t",$snp[0],$snp[1]);
		my $primary = $snp[5];
		my $secondary = $snp[6];

		my @out = ();
		push @out, "$col";
	
		my $length = scalar(@snp) - 1;
	
		for ( my $i=10; $i<=$length; $i++ ) {
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
if ($tools =~ "T" or $tools =~ "t"){

	print "\nConverting the SNP genotyping matrix into an input file compatible with TASSEL GUI software...";

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

		my $col1 = join ("_",$snp[0],$snp[1]);
		my $alleles = join ("/",$snp[4],$snp[5]);

		my @out;

		my $length = scalar(@snp) - 1;

		for ( my $i=10; $i<=$length; $i++ ) {
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
if ($tools =~ "P" or $tools =~ "p"){

	print "\nConverting the SNP genotyping matrix into an input file compatible with Plink software...";

	open my $IN3, "<", "$SNPmatrix" or die "Can't load $SNPmatrix file";
	open my $OUT_P, ">", "$Plink_out" or die "Can't initialize $Plink_out output file";

	while (<$IN3>) {
		chomp;
		my @input = split ("\t", $_);

		my $header = $input[0];
		my $position = $input[1];
		my $snp_identifier = join ("_",$header,$position);
		$header =~ s/chr0//;
		$header =~ s/chr//;
		my $Tpedfile = join("\t",$header,$snp_identifier,"0",$position);

		my @geno = splice @input, 11; 
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

print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
