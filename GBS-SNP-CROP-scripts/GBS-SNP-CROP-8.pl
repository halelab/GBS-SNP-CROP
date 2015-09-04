#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP-8.pl
# For description, see Melo et al. (2015) DOI XXX
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-8.pl -in <input file: SNP genotyping matrix> -b <barcodesID file>
-formats <The name(s) of the software packages (R, Tassel, Plink) for which an input-compatible file format should be created>.\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($SNPmatrix,$barcodesID_file,$tools);

GetOptions(
'in=s' => \$SNPmatrix,    	# file
'b=s' => \$barcodesID_file,     # file
'formats=s' => \$tools,	# file
) or die "Error in command line arguments.\n$Usage\n$Manual\n";

#####################
# Defining outputs
#####################

my $R_out = join (".","GBS-SNP-CROP","R_in","txt");
my $Tassel_out = join (".","GBS-SNP-CROP","Tassel_in","hmp","txt");
my $Plink_out = join (".","GBS-SNP-CROP","PLINK_in","tped");

#############################
# Creating R input file
#############################
if ($tools =~ "R" or $tools =~ "r"){

	print "\nConverting the SNP genotyping matrix into an input file compatible with R software...";

	open (IN, "$SNPmatrix") || die "cant load file";
	open (OUT_R, ">$R_out") || die "cant load file";
	
	while(<IN>) {
		my @snp = split "\t", $_;
		chomp @snp;

		my $col = join ("\t",$snp[0],$snp[1]);
		my $primary = $snp[4];
		my $secondary = $snp[5];

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
			} elsif ($geno2 eq join("/",$secondary,$secondary)) {
				push @out, "1";
				next;
			} else {
				push @out, "NA";
				next;
			}
		}	
		my $line = join "\t", @out;
		print OUT_R "$line\n";
	}

	close IN;
	close OUT_R;
	print "DONE.\n";
	
	} else {
}

#######################################
# Creating Tassel HapMap input file
#######################################
if ($tools =~ "T" or $tools =~ "t"){
	
	print "\nConverting the SNP genotyping matrix into an input file compatible with TASSEL GUI software...";

	my @TaxaNames = ();
	open BAR, "$barcodesID_file" or die "Can't find $barcodesID_file file\n";
	while(<BAR>) {
		my $barcodesID = $_;
		chomp $barcodesID;
		my @barcode = split("\t", $barcodesID);
		push @TaxaNames, $barcode[1];
	}
	s/_/./g for @TaxaNames;

	my @output;
	
	open (IN, "$SNPmatrix") || die "cant load file";
	open (OUT_T, ">$Tassel_out") || die "cant load file";
	
	while(<IN>) {
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
		my $line = "$col1\t$alleles\t$snp[0]\t$snp[1]\t+\tNA\tGBS-SNP-CROP\tGBS\tRefName\tCustom\tQC+\t@out";
		push @output,"$line\n";

	}
	print OUT_T "#rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t@TaxaNames\n@output";

	close IN;
	close OUT_T;
	print "DONE.\n";
	
	} else {
}

#######################################
# Creating Tassel HapMap input file
#######################################
if ($tools =~ "P" or $tools =~ "p"){

	print "\nConverting the SNP genotyping matrix into an input file compatible with PLINK software...";
	
	open (IN, "$SNPmatrix") || die "cant load file";
	open (OUT_P, ">$Plink_out") || die "cant load file";
	
	while (<IN>) {
		chomp;
		my @input = split ("\t", $_);
		my $header = $input[0];
		my $position = $input[1];
		my $snp_identifier = join ("_",$header,$position);
		$header =~ s/chr0//;
		$header =~ s/chr//;
		my $Tpedfile = join("\t",$header,$snp_identifier,"0",$position);
	
		my @geno = splice @input, 10; 
		foreach (@geno) {
			my @geno1 = split /\|/, $_;
			s/-/0 0/ for @geno1;
			s/\// / for @geno1;	
			$Tpedfile = join("\t",$Tpedfile,$geno1[0]);
		}
		print OUT_P "$Tpedfile\n";
	}

	close IN;
	close OUT_P;
	print "DONE.\n";
	
	} else {
}

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;