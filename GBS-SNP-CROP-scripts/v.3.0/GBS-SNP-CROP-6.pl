#!/usr/bin/perl

#########################################################################################################################
# GBS-SNP-CROP, Step 6. For description, see User Manual (https://github.com/halelab/GBS-SNP-CROP/wiki)
#########################################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw( min max );
use Parallel::ForkManager;

my $Usage = "Usage: perl GBS-SNP-CROP-6.pl -b <barcodesID file> -out <SNPs master matrix file> -indels <Use this flag to identify SNPs + Indels> -t <threads>\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($barcodesID_file,$output_file,$indels,$threads);

GetOptions(
'b=s' => \$barcodesID_file,		# file
'out=s' => \$output_file,		# file
'indels' => \$indels,	    	# string
't=s' => \$threads,				# numeric
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 6, v.3.0\n#################################\n\n";
my $pm = new Parallel::ForkManager($threads);
my $sttime = time;

my @files = ();
open my $IN, "<", "$barcodesID_file" or die "Can't find $barcodesID_file file\n";
while(<$IN>) {
	my $barcodesID = $_;
	chomp $barcodesID;
	my @barcode = split("\t", $barcodesID);
	my $barcode_list = $barcode[0];
	my $TaxaNames = $barcode[1];

	push @files, $TaxaNames;
}
close $IN;

open my $CountList, ">", "CountFileList.txt" or die "Unable to load CountFileList.tx file\n";

###################################
# Identifying both SNPs and indels 
###################################
if ($indels) {
	print "The GBS-SNP-CROP will identify and call both SNPs and indels from the alignment dataset.\n\n";
	print "Counting nucleotides and filtering monomorphic sites for all genotypes. Even using multi-threads this process can take a while ... be patient!\n";

	foreach my $file (@files) {
	    	my $mpileup_input = join (".", "$file","mpileup");
		my $count_out = join (".", "$file","count","txt");
		my $ref_file = join (".", "$file","ref","txt");
		print $CountList "$count_out\n";

		my $pid = $pm->start and next;
			
		open my $PILEUP, "<", "$mpileup_input" or die "Can't load $mpileup_input file";
		open my $OUT1, ">", "$count_out" or die "Can't initialized $count_out ouput file";
		open my $REF, ">", "$ref_file" or die "Unable to identify $ref_file\n";

		while (<$PILEUP>) {
	
			my @input1 = split("\t", $_);
			my $ref = $input1[2];
			my $algn = $input1[4];		
			$algn =~ s/\^.\././g;
			
			my @positions;
			my @sizes;
			my @Indels;
			while ( $algn =~ /[\.|,]{1}[\+|-]{1}(\d+)/g  ) { 
				push @positions, pos($algn);
				push @sizes, $1;
			}

			my $indices = scalar @positions - 1;

			for (my $k = $indices; $k >= 0; $k--) {
				my $indel = substr ( $algn, $positions[$k] - length($sizes[$k]) - 1, 1 + length($sizes[$k]) + $sizes[$k]);
				$indel = uc($indel);
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
	
			# Work on a nucleotides specific string
			$algn =~ s/\$//g;   
			$algn =~ s/\*//g;
		
			my $uc_ref = uc $ref;
			my $lc_ref = lc $ref;
			$algn =~ s/\./$uc_ref/g;
			$algn =~ s/\,/$lc_ref/g;
		
			my @bases = split(//, $algn);
			@bases = grep /\S/, @bases;
		
			my $A = 0; my $a = 0; my $xA = 0;
			my $C = 0; my $c = 0; my $xC = 0;	
			my $G = 0; my $g = 0; my $xG = 0;
			my $T = 0; my $t = 0; my $xT = 0;
		
			for (my $x=0; $x<scalar(@bases);$x++) {
				if ($bases[$x] =~ /A/){
					$A++;
				}
				if ($bases[$x] =~ /a/){
					$a++;
				}
				if ($bases[$x] =~ /C/){
					$C++;
				}
				if ($bases[$x] =~ /c/){
					$c++;
				}
				if ($bases[$x] =~ /G/){
					$G++;
				}
				if ($bases[$x] =~ /g/){
					$g++;
				}
				if ($bases[$x] =~ /T/){
					$T++;
				}
				if ($bases[$x] =~ /t/){
					$t++;
				}
			}
		
			if ($A >= $a) {
				$xA = $A;
			} else {
				$xA = $a
			}
		
			if ($C >= $c) {
				$xC = $C;
			} else {
				$xC = $c
			}
		
			if ($G >= $g) {
				$xG = $G;
			} else {
				$xG = $g
			}
		
			if ($T >= $t) {
				$xT = $T;
			} else {
				$xT = $t
			}
		
			if (scalar @IndelCount == 0) { 	
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","_","_","_","_"),"\n";
			
				if ( (($xA + $xC) * ($xG + $xT)) > 0 or (($xA + $xG) * ($xC + $xT)) > 0 ) {
					print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
				}
				next;
		
			} elsif (scalar @IndelCount == 1) {
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","$IndelCount[0]","$IndelType[0]","_","_"),"\n";
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
					
			} elsif (scalar @IndelCount > 1) {
				print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT","$IndelCount[0]","$IndelType[0]","$IndelCount[1]","$IndelType[1]"),"\n";
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
			
			} else {
				next;
			}
		}
		close $PILEUP;
		close $OUT1;
		close $REF;
		$pm->finish;
	}
	$pm->wait_all_children;
	
	print "DONE.\n";

	# Creating a master list of all putative variant positions in the population
	system ( "cat *.ref.txt | uniq > VerticalRefPos.txt" );

	my $posFile = "VerticalRefPos.txt";
	my $countList = "CountFileList.txt";
	open my $POS, "<", "$posFile" || die "Can't load file $!";
	open my $LIST, "<","$countList" || die "Can't load file $!";
	open my $DEST, ">", "$output_file" || die "Can't load file $!";

	my %posHash;
	print "\nCreating a comprehensive master matrix, with genotype-specific alignment summaries,for all putative variants positions ... ";

	while (<$POS>){
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
	close $POS;

	while ( my $fileName = <$LIST> ) {
		chomp $fileName;
		open my $GENO_BASE_COUNT_FILE, "<", "$fileName" or die "Can't load file $!";

		my %genoHash;

		while ( my $line = <$GENO_BASE_COUNT_FILE>){
			my @input5 = split ("\t", $line);
			chomp @input5;
			my $chr_pos_ref = join ("\t","$input5[0]", "$input5[1]","$input5[2]");
			$genoHash{$chr_pos_ref} = $input5[3];
			@input5 = ();
		}

		foreach my $chr_pos_ref ( keys %posHash ){
			if ( $genoHash{$chr_pos_ref} ) {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "$genoHash{$chr_pos_ref}");
			} else {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_,_,_,_,_,");
			}
		}

		%genoHash = ();
		close $GENO_BASE_COUNT_FILE;
	}
	close $LIST;

	foreach my $key ( sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %posHash ) {
		print $DEST "$posHash{$key}\n";
	}
	close $DEST;
	print "\nDONE.\n";
	
#########################
# Identifying SNPs only
#########################

} else {
	print "The GBS-SNP-CROP will identify and call only SNPs from the alignment dataset.\n\n";
	print "Counting nucleotides and filtering monomorphic sites for all genotypes. Even using multi-threads this process can take a while ... \n";


	foreach my $file (@files) {
    		my $mpileup_input = join (".", "$file","mpileup");
		my $count_out = join (".", "$file","count","txt");
		my $ref_file = join (".", "$file","ref","txt");
		print $CountList "$count_out\n";

		my $pid = $pm->start and next;
	
		open my $PILEUP, "<", "$mpileup_input" or die "Can't load $mpileup_input file";
		open my $OUT1, ">", "$count_out" or die "Can't initialized $count_out ouput file";
		open my $REF, ">", "$ref_file" or die "Unable to identify $ref_file\n";
	
		while (<$PILEUP>) {	
			my @input1 = split("\t", $_);
			my $ref = $input1[2];
			my $algn = $input1[4];		
			$algn =~ s/\^.\././g;
			
			my @positions;
			my @sizes;
			my @Indels;
			while ( $algn =~ /[\.|,]{1}[\+|-]{1}(\d+)/g  ) { 
				push @positions, pos($algn);
				push @sizes, $1;
			}

			my $indices = scalar @positions - 1;

			for (my $k = $indices; $k >= 0; $k--) {
				my $indel = substr ( $algn, $positions[$k] - length($sizes[$k]) - 1, 1 + length($sizes[$k]) + $sizes[$k]);
				$indel = uc($indel);
				push @Indels, $indel;
				my $start = substr ( $algn, 0, $positions[$k] - length($sizes[$k]) - 2 );
				my $end = substr ( $algn, $positions[$k] + $sizes[$k] );
				$algn = join ("", "$start", "$end" );
			}

			$algn =~ s/\$//g;   
			$algn =~ s/\*//g;
		
			my $uc_ref = uc $ref;
			my $lc_ref = lc $ref;
			$algn =~ s/\./$uc_ref/g;
			$algn =~ s/\,/$lc_ref/g;
		
			my @bases = split(//, $algn);
			@bases = grep /\S/, @bases;
		
			my $A = 0; my $a = 0; my $xA = 0;
			my $C = 0; my $c = 0; my $xC = 0;	
			my $G = 0; my $g = 0; my $xG = 0;
			my $T = 0; my $t = 0; my $xT = 0;
		
			for(my $x=0; $x<scalar(@bases);$x++) {
				if ($bases[$x] =~ /A/){
					$A++;
				}
				if ($bases[$x] =~ /a/){
					$a++;
				}
				if ($bases[$x] =~ /C/){
					$C++;
				}
				if ($bases[$x] =~ /c/){
					$c++;
				}
				if ($bases[$x] =~ /G/){
					$G++;
				}
				if ($bases[$x] =~ /g/){
					$g++;
				}
				if ($bases[$x] =~ /T/){
					$T++;
				}
				if ($bases[$x] =~ /t/){
					$t++;
				}
			}
		
			if ($A >= $a) {
				$xA = $A;
			} else {
				$xA = $a
			}
		
			if ($C >= $c) {
				$xC = $C;
			} else {
				$xC = $c
			}
		
			if ($G >= $g) {
				$xG = $G;
			} else {
				$xG = $g
			}
		
			if ($T >= $t) {
				$xT = $T;
			} else {
				$xT = $t
			}
			
			print $OUT1 join ("\t",$input1[0],$input1[1],$input1[2]),"\t",join(",","$xA","$xC","$xG","$xT"),"\n";
		
			if ( (($xA + $xC) * ($xG + $xT)) > 0 or (($xA + $xG) * ($xC + $xT)) > 0 ) {
				print $REF join ("\t", $input1[0], $input1[1], $input1[2]),"\n";
			}
		}
		close $PILEUP;
		close $OUT1;
		close $REF;
		$pm->finish;
	}
	$pm->wait_all_children;

	print "DONE.\n";
	# Creating a master list of all putative variant positions in the population
	system ( "cat *.ref.txt | uniq > VerticalRefPos.txt" );
	
	my $posFile = "VerticalRefPos.txt";
	my $countList = "CountFileList.txt";
	open my $POS, "<", "$posFile" || die "Can't load file $!";
	open my $LIST, "<","$countList" || die "Can't load file $!";
	open my $DEST, ">", "$output_file" || die "Can't load file $!";

	my %posHash;
	print "\nCreating a comprehensive master matrix, with genotype-specific alignment summaries,for all putative variants positions ... \n";

	while (<$POS>){
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
	close $POS;

	while ( my $fileName = <$LIST> ) {
		chomp $fileName;
		open my $GENO_BASE_COUNT_FILE, "<", "$fileName" or die "Can't load file $!";

		my %genoHash;

		while ( my $line = <$GENO_BASE_COUNT_FILE>){
			my @input5 = split ("\t", $line);
			chomp @input5;
			my $chr_pos_ref = join ("\t","$input5[0]", "$input5[1]","$input5[2]");
			$genoHash{$chr_pos_ref} = $input5[3];
			@input5 = ();
		}

		foreach my $chr_pos_ref ( keys %posHash ){
			if ( $genoHash{$chr_pos_ref} ) {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "$genoHash{$chr_pos_ref}");
			} else {
				$posHash{$chr_pos_ref} = join ("\t", "$posHash{$chr_pos_ref}", "_,_,_,_");
			}
		}

		%genoHash = ();
		close $GENO_BASE_COUNT_FILE;
	}
	close $LIST;

	foreach my $key ( sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %posHash ) {
		print $DEST "$posHash{$key}\n";
	}
	close $DEST;
	print "DONE.\n";
}

sub main {
	my $dir = "mpileup"; 
	unless(-e $dir, or mkdir $dir) {die "Directory $dir does not exist and cannot be created.\n";}
	}
main();

system ( "mv *.count.txt *.mpileup VerticalRefPos.txt CountFileList.txt ./mpileup" );
system ( "rm *.ref.txt" );

print "\nThe master matrix was successfully created. Please, proceed with the filtering and genotyping step.\n";
print "Elapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
