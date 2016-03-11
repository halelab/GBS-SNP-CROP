#!/usr/bin/perl

##########################################################################################################################
# GBS-SNP-CROP, Step 7. For description, please see Melo et al. (2016) BMC Bioinformatics DOI 10.1186/s12859-016-0879-y.
##########################################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum);

my $Usage = "Usage: perl GBS-SNP-CROP-7.pl -in <master matrix filename> -out <genotyping matrix filename> -mnHoDepth0 <minimum depth required for calling a homozygote when the alternative allele depth = 0>\n"
."-mnHoDepth1 <minimum depth required for calling a homozygote when the alternative allele depth = 1> -mnHetDepth <minimum depth required for each allele when calling a heterozygote>\n"
."-altStrength <Across the population for a given putative bi-allelic SNP, this alternate allele strength parameter is the minimum proportion of non-primary allele reads that are the secondary allele (numeric)>\n" 
."-mnAlleleRatio <minimum requierd ratio of less frequent to more frequent allele depht to more frequent allele depth> -mnCall <minimum acceptable proportion of genotyped individuals to retain a SNP>\n"
."-mnAvgDepth <minimum average depth of an acceptable SNP> -mxAvgDepth <maximum average depth of an acceptable SNP>\n";
my $Manual = "Please see UserManual on GBS-SNP-CROP GitHub page (https://github.com/halelab/GBS-SNP-CROP.git) or the original manuscript: Melo et al. (2016) BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n"; 

my ($SummaryFile,$output,$minHomoDepth,$minHomoDepthOneAlt,$minHeteroDepth,$AltAlleleStrength,$AlleleFreqProp,$CallRate,$minAvgDepth,$maxAvgDepth);

GetOptions(
'in=s' => \$SummaryFile,                   # file
'out=s' => \$output,                       # file
'mnHoDepth0=s' => \$minHomoDepth,          # numeric
'mnHoDepth1=s' => \$minHomoDepthOneAlt,    # numeric
'mnHetDepth=s' => \$minHeteroDepth,        # numeric
'altStrength=s' => \$AltAlleleStrength,    # numeric
'mnAlleleRatio=s' => \$AlleleFreqProp,     # numeric
'mnCall=s' => \$CallRate,                    # numeric
'mnAvgDepth=s' => \$minAvgDepth,             # numeric
'mxAvgDepth=s' => \$maxAvgDepth,           # numeric
) or die "$Usage\n$Manual\n";

print "\n#################################\n# GBS-SNP-CROP, Step 7, v.1.1\n#################################\n";
print "\nGBS-SNP-CROP is filtering SNPs and calling genotypes for your population...\n";

open my $VAR, "<", "$SummaryFile" or die "Can't load file $SummaryFile";
open my $DEST, ">", "$output" or die "Can't initialized $output output file";

my $lc = 0;

while (<$VAR>) {
	chomp;
	my @input = split("\t", $_);
	my $header = $input[0];
	my $position = $input[1];
	my $ref = $input[2];

	my %Var_depth = ( 'A' => 0,
					'C' => 0,
					'G' => 0,
					'T' => 0 );

	my %Var_cnt = ( 'A' => 0,
					 'C' => 0,
					 'G' => 0,
					 'T' => 0 );

	for ( my $x=3; $x < scalar(@input); $x++ ){
		my $geno=$input[$x];
		my @genotype = split(",", $geno);
		
		if ( $genotype[0] ne "_" ) {
			$Var_depth{'A'} = $Var_depth{'A'} + $genotype[0];
			if ($genotype[0] > 0) {
				$Var_cnt{'A'} = $Var_cnt{'A'} + 1;
			}
		}
		if ( $genotype[1] ne "_" ) {
			$Var_depth{'C'} = $Var_depth{'C'} + $genotype[1];
			if ($genotype[1] > 0) {
				$Var_cnt{'C'} = $Var_cnt{'C'} + 1;
			}
		}
		if ( $genotype[2] ne "_" ) {
			$Var_depth{'G'} = $Var_depth{'G'} + $genotype[2];
			if ($genotype[2] > 0) {
				$Var_cnt{'G'} = $Var_cnt{'G'} + 1;
			}
		}
		if ( $genotype[3] ne "_" ) {
			$Var_depth{'T'} = $Var_depth{'T'} + $genotype[3];
			if ($genotype[3] > 0) {
				$Var_cnt{'T'} = $Var_cnt{'T'} + 1;
			}
		}
	}

	my $total;
	while ( my ($key, $val) = each %Var_depth ) {
		$total += sum $val;
	}
	if ( $total == 0) {
		next;
	}

	my @rank = ();

	foreach my $variant (sort { $Var_depth{$b} <=> $Var_depth{$a} } keys %Var_depth) {
		push @rank, $variant;
	}

	my $pop_one_var = $rank[0];
	my $pop_two_var = $rank[1];
	my $pop_one_count = $Var_depth{$pop_one_var};
	my $pop_two_count = $Var_depth{$pop_two_var};

	# Initial filters based on population-level allele frequencies
	my $alt_allele_ratio = ( $pop_two_count / ($pop_one_count + $pop_two_count) );
	if ( $alt_allele_ratio < 0.05) {
		next;
	}
	
	if ( $pop_two_count == 0 ) {
		next;
	}
	
	my $alt_allele_strength = ( $pop_two_count / ($total - $pop_one_count) );
	if ( $alt_allele_strength < $AltAlleleStrength ) {
		next;
	}
	
	my $N1 = $Var_cnt{$pop_one_var};
	my $N2 = $Var_cnt{$pop_two_var};
	if ( (scalar(@input) - 3) > 20 and ($N2 < 3) ) {
		next;
	} elsif( (scalar(@input) - 3) < 20 and ($N2 < 2) ) {
		next;
	}
	
	# Secondary filters based on genotype-level allele frequencies
	my $row = "";
	my $homo_pri = 0;
	my $hets = 0;
	my $homo_alt = 0;
	my $cumulative_depth = 0;

	for ( my $x=3; $x < scalar(@input); $x++ ) {
	
		my $geno=$input[$x];
		my @genotype = split(",",$geno);
	
		my %geno_hash = ( 'A' => 0,
							'C' => 0,
							'G' => 0,
							'T' => 0 );

		if ( $genotype[0] ne "_" ) {
			$geno_hash{'A'} = $genotype[0];
		}
		if ( $genotype[1] ne "_" ) {
			$geno_hash{'C'} = $genotype[1];
		}
		if ( $genotype[2] ne "_" ) {
			$geno_hash{'G'} = $genotype[2];
		}
		if ( $genotype[3] ne "_" ) {
			$geno_hash{'T'} = $genotype[3];
		}

		my $primary_count = $geno_hash{$pop_one_var};
		my $alt_count = $geno_hash{$pop_two_var};

		if (!$primary_count) {
			$primary_count = 0;
		}
		if ($primary_count eq "_" ) {
			$primary_count = 0;
		} 
		if (!$alt_count) {
			$alt_count = 0;
		}
		if ($alt_count eq "_" ) {
			$alt_count = 0;
		} 

		my $depth = 0;

		# Calling heterozygotes
		if ( ($primary_count > $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($primary_count > $alt_count) and $primary_count >= 20 and $alt_count >= ($AlleleFreqProp * $primary_count)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($alt_count > $primary_count) and ($alt_count < 20) and ($alt_count >= $minHeteroDepth && $primary_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($alt_count > $primary_count) and $alt_count >= 20 and $primary_count >= ($AlleleFreqProp * $alt_count)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;

		# Calling homozygotes
		} elsif ( $primary_count >= $minHomoDepth && $alt_count == 0 ) {
			$homo_pri++;
			$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/0";
			$depth = $primary_count;
		} elsif ( $primary_count >= $minHomoDepthOneAlt && $alt_count == 1 ) {
			$homo_pri++;
			$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/1";
			$depth = $primary_count;
		} elsif ( $alt_count >= $minHomoDepth && $primary_count == 0 ) {
			$homo_alt++;
			$row = "$row\t$pop_two_var/$pop_two_var|0/$alt_count";
			$depth = $alt_count;

		# Declaring missing genotypes
		} elsif ( $primary_count == $minHeteroDepth || $alt_count == $minHeteroDepth ) {
			$row = "$row\t-|$primary_count/$alt_count";  
		} elsif ( $primary_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
			$row = "$row\t-|$primary_count/$alt_count";        
		} elsif ( $alt_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
			$row = "$row\t-|$primary_count/$alt_count";

		} else {
			$row = "$row\t-|-";
		}
			$cumulative_depth = $cumulative_depth + $depth;
	}

	# Calculating population-level genotype parameters for additional filtering
	my $scored_genotypes = $homo_pri + $hets + $homo_alt;
	my $percentage_scored_genotypes = sprintf("%.2f",($scored_genotypes) / ((scalar(@input) - 3))*100);

	if ( $scored_genotypes == 0 ) {
		next;
	}

	my $avgDepth = sprintf("%.2f", ($cumulative_depth / $scored_genotypes));

	# Depth filter
	if( ($avgDepth < $minAvgDepth) or ($avgDepth > $maxAvgDepth) ) {
		next;
	}
	
	# Genotype representation filter
	if ( (scalar(@input) - 3) > 20 and (($homo_pri + $hets) < 3 or ($homo_alt + $hets) < 3 ) ) {
		next;
	} elsif( (scalar(@input) - 3) < 20 and (($homo_pri + $hets) < 2 or ($homo_alt + $hets) < 2 ) ) {
		next;
	}

	# Genotype proportion filter
	if($scored_genotypes <= $CallRate * (scalar @input - 3)){
		next;
	}

	$row = join ("\t","$header","$position","$ref","$avgDepth","$pop_one_var","$pop_two_var","$percentage_scored_genotypes","$homo_pri","$hets","$homo_alt$row" );
	print $DEST "$row\n";
	$lc++;
}

close $DEST;

print "\nYour $output genotyping matrix was successfully created into SNPsCalled directory. \n";

print "GBS-SNP-CROP called $lc SNPs in your population using the following parameters:\n"
."Minimum read depth required to call homozygotes when the secondary allele depth = 0: $minHomoDepth\n"
."Minimum read depth required to call homozygotes when the secondary allele depth = 1: $minHomoDepthOneAlt\n"
."Minimum read depth for each allele equired to call heterozygotes: $minHeteroDepth\n"
."Minimum proportion of non-primary allele reads that are the secondary allele: $AltAlleleStrength\n"
."Minimum proportion of less frequent allele depth to more frequent allele depth: $AlleleFreqProp\n"
."Minimum proportion of genotyped individuals to accept a SNP: $CallRate\n"
."Minimum acceptable average read depth: $minAvgDepth\n"
."Maximum acceptable average read depth: $maxAvgDepth\n";

sub main {
	my $dir = "SNPsCalled"; 
	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
}
main();

system ( "mv $SummaryFile $output ./SNPsCalled" );

print "\n\nPlease cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

exit;
