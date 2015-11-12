#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 6. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-6.pl -b <barcodesID file> -out <output variant master matrix file>\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($barcodesID_file,$output_file);

GetOptions(
'b=s' => \$barcodesID_file,    # file
'out=s' => \$output_file,      # file
) or die "$Usage\n$Manual\n";

my @files = ();

open IN, "$barcodesID_file" or die "Can't find $barcodesID_file file\n";
	
while(<IN>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];
	
	push @files, $TaxaNames;
}

chomp (@files);

foreach my $file (@files) {
    my $mpileup_input = join (".", "$file","mpileup");
	my $count_out = join (".", "$file","count","txt");
			
	print "\n\nParsing $mpileup_input file...\n";		
	
	open (PILEUP, "$mpileup_input") || die "cant load $mpileup_input file";
	open (OUT1, ">$count_out") || die "cant load file";
		
	while (<PILEUP>){
		my @input1 = split("\t", $_);
		my $ref = $input1[2];
		my $algn = uc $input1[4];
		$algn =~ s/\,/\./g;
		$algn =~ s/\^.\././g;
	
	### Capturing Indels and put it out from mpileup string
		my @positions;
		my @sizes;
		my @Indels;
		while ( $algn =~ /\.{1}[\+|-]{1}(\d+)/g  ) { 
			push @positions, pos($algn);
			push @sizes, $1;
		}

		my $indices = scalar @positions - 1;

		for (my $k = $indices; $k >= 0; $k--) {
			my $indel = substr ( $algn, $positions[$k] - length($sizes[$k]) - 1, 1 + length($sizes[$k]) + $sizes[$k]);
			push @Indels, $indel;
			my $start = substr ( $algn, 0, $positions[$k] - length($sizes[$k]) - 2 );
			my $end = substr ( $algn, $positions[$k] + $sizes[$k] );
			$algn = join ("", "$start", "$end" );
		}

		# Count Indels frequency and store a hash with Indel type and count	
		my %Indel_cnt;
		$Indel_cnt{$_}++ foreach @Indels;
		
		my @IndelType = sort {$Indel_cnt{$b} <=> $Indel_cnt{$a}} keys %Indel_cnt;
		my @IndelCount = @Indel_cnt{@IndelType};
	
		### Work on a nucleotides specific string
		$algn =~ s/\$//g;   
		$algn =~ s/\*//g;
		$algn =~ s/\./$ref/g;
		
		my @bases = split(//, $algn);
		@bases = grep /\S/, @bases; # remove empty elements from array
		
		my $A=0;	my $C=0; 
		my $G=0;	my $T=0; 
		
        	for(my $x=0; $x<scalar(@bases);$x++) {
			if ($bases[$x] =~ /A/){
                   	$A++;
			}
            		if ($bases[$x] =~ /C/){
                	 $C++;
            		}
            		if ($bases[$x] =~ /G/){
                	    $G++;
            		}
            		if ($bases[$x] =~ /T/){
                	    $T++;
        		}
        	}
		if (scalar @IndelCount == 0) { 	
			print OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$A","$C","$G","$T","_","_","_","_"),"\n";
			next;
		} elsif (scalar @IndelCount == 1) {
			print OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$A","$C","$G","$T","$IndelCount[0]","$IndelType[0]","_","_"),"\n";
			next;
		} elsif (scalar @IndelCount > 1) {
			print OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$A","$C","$G","$T","$IndelCount[0]","$IndelType[0]","$IndelCount[1]","$IndelType[1]"),"\n";
			next;
		} else {
		next;
		}
	}
	
	close PILEUP;
	close OUT1;
	
	print "Filtering polymorphic sites from $count_out file...\n";
	
	my $countF_out = join (".", "$file","countF","txt");
	my $ref_file = join (".", "$file","ref","txt");
		
	open (COUNT, "$count_out") || die "Unable to identify $count_out\n";
	open (REF, ">$ref_file") || die "Unable to identify $ref_file\n";
	open (OUT2, ">$countF_out") || die "cant load file";
	
	my @output;
	while(<COUNT>) {
		my @input2 = split "\t", $_;
        chomp @input2;
        my $ref = $input2[2];
        my @geno = split(",", $input2[3]);
	    
	if ($ref eq "A" and $geno[0] == 0 and ($geno[1]!= 0 or $geno[2] != 0 or $geno[3] != 0) or
		$ref eq "A" and $geno[0] != 0 and ($geno[1]!= 0 or $geno[2] != 0 or $geno[3] != 0) or
     		$ref eq "A" and $geno[0] != 0 and ($geno[1]!= 0 or $geno[2] != 0 or $geno[3] != 0) or
	        $ref eq "C" and $geno[1] == 0 and ($geno[0]!= 0 or $geno[2] != 0 or $geno[3] != 0) or
    		$ref eq "C" and $geno[1] != 0 and ($geno[0]!= 0 or $geno[2] != 0 or $geno[3] != 0) or
       		$ref eq "G" and $geno[2] == 0 and ($geno[0]!= 0 or $geno[1] != 0 or $geno[3] != 0) or
         	$ref eq "G" and $geno[2] != 0 and ($geno[0]!= 0 or $geno[1] != 0 or $geno[3] != 0) or   	
        	$ref eq "T" and $geno[3] == 0 and ($geno[0]!= 0 or $geno[1] != 0 or $geno[2] != 0) or
    		$ref eq "T" and $geno[3] != 0 and ($geno[0]!= 0 or $geno[1] != 0 or $geno[2] != 0) ) {

       		@output = join ("\t", @input2);
    		print OUT2 "@output\n";
    		print REF join ("\t", $input2[0], $input2[1], $input2[2]),"\n";
    		}
	}
	print "DONE.\n";
	close COUNT;
	close OUT2;
}

print "\nCreating a polymorphic variant list...\n";

system ( "cat *.ref.txt | uniq > SNPsVerticalRef.txt" );

# Creating a count file list
system ( "ls *.count.txt > CountFileList.txt" );

print "DONE.\n";

print "\nCreating a variant depth count master matrix...\n";

my $posFile = "SNPsVerticalRef.txt";
my $countList = "CountFileList.txt";

open (POS, "$posFile") || die "cant load file $!";
open (LIST, "$countList") || die "cant load file $!";
open (DEST, ">$output_file") || die "cant load file $!";

my %posHash;

while (<POS>){
	chomp;
	my @input4 = split("\t", $_);
	my $chr_pos_ref = join("\t", "$input4[0]", "$input4[1]", "$input4[2]");
	if ( $posHash{$chr_pos_ref} ) {
		next;
	} else {
		$posHash{$chr_pos_ref} = $chr_pos_ref;
		next;
	}

}
close POS;

while ( my $fileName = <LIST> ) {
	chomp $fileName;
	open (GENO_BASE_COUNT_FILE, "$fileName") or die "can't load file $!";

	my %genoHash;

	while ( my $line = <GENO_BASE_COUNT_FILE>){
		my @input5 = split ("\t", $line);
		chomp @input5;
		my $chr_pos_ref = join ("\t", "$input5[0]", "$input5[1]", "$input5[2]");
		$genoHash{$chr_pos_ref} = $input5[3];
		@input5 = ();
	}

	foreach my $chr_pos_ref ( keys %posHash ){
		if ( $genoHash{$chr_pos_ref} ) {
			$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "$genoHash{$chr_pos_ref}");
		} else {
			$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_,_,_,_,_");
		}
	}

	%genoHash = ();
	close GENO_BASE_COUNT_FILE;
}

close LIST;

foreach my $key ( sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %posHash ) {
	print DEST "$posHash{$key}\n";
}

close DEST;

print "DONE.\n\n";

sub main {
   	my $dir = "mpileup"; 
   	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
   	}
main();

system ( "mv *.count.txt ./mpileup" );
system ( "mv *.mpileup ./mpileup" );
system ( "mv VariantVerticalRef.txt ./mpileup" );
system ( "mv CountFileList.txt ./mpileup" );
system ( "rm *.ref.txt" );
system ( "rm *.countF.txt" );

print "The variant count master matrix was successfully created.\nPlease, proceed with filtering and genotyping call on step 7.\n";

print "\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;
