#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 6. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-6.pl -b <barcodesID file> -out <output SNP master matrix file>\n";
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
		my $algn = $input1[4];
	
		my @bases = split(//, $algn);
		s/\./$ref/g for @bases;
		s/\,/$ref/g for @bases;
			
		my $aCount=0; my $ACount=0;
		my $gCount=0; my $GCount=0;
		my $cCount=0; my $CCount=0;
		my $tCount=0; my $TCount=0;

        my $x;
        for($x=0;$x<scalar(@bases);$x++) {
			if ($bases[$x] =~ /A/){
                   $ACount++;
			}
        	if ($bases[$x] =~ /a/){
                   $aCount++;
            }
            if ($bases[$x] =~ /C/){
                   $CCount++;
            }
            if ($bases[$x] =~ /c/){
                   $cCount++;
            }
            if ($bases[$x] =~ /G/){
                    $GCount++;
            }
            if ($bases[$x] =~ /g/){
                    $gCount++;
            }
            if ($bases[$x] =~ /T/){
                    $TCount++;
            }
            if ($bases[$x] =~ /t/){
                    $tCount++;
            }
        }
        my $A = $ACount + $aCount;
        my $C = $CCount + $cCount;
        my $G = $GCount + $gCount;
        my $T = $TCount + $tCount;
        	          
		print OUT1 join ("\t", $input1[0], $input1[1], $input1[2]),"\t", join(",","$A", "$C", "$G", "$T"),"\n";
		
		}
	
	close PILEUP;
	close OUT1;
	
	print "Filtering polymorphic SNPs from $count_out file...\n";
	
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

print "\nCreating a polymorphic SNP list...\n";

system ( "cat *.ref.txt | uniq > SNPsVerticalRef.txt" );

# Creating a count file list
system ( "ls *.count.txt > CountFileList.txt" );

print "DONE.\n";

print "\nCreating a SNP master matrix...\n";

my $posFile = "SNPsVerticalRef.txt";
my $countList = "CountFileList.txt";
#my $outFile = "SNPsDiscovered_Master_Matrix.txt";

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

while ( my $fileName = <LIST> )
{
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
			$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_,_,_");
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
system ( "mv SNPsVerticalRef.txt ./mpileup" );
system ( "mv CountFileList.txt ./mpileup" );
system ( "rm *.ref.txt" );
system ( "rm *.countF.txt" );

print "The SNP master matrix was successfully created.\nPlease, proceed with SNP filtering and genotyping step 7.\n";

print "\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;
