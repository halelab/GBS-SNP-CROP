#!/usr/bin/perl

###########################################################################################################
# GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline
#
# Authors: Arthur Melo, Radhika Bartaula, Iago Hale
# Department of Agriculture, Nutrition, and Food Systems, University of New Hampshire, Durham, NH, 03824
#
# A detailed description can be found at https://github.com/halelab/GBS-SNP-CROP
# 
# For help: perl GBS-SNP-CROP-7.pl
###########################################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum);

######################
# The help function
######################
my $help = $ARGV[0];
my ($SummaryFile,$output,$type,$minHomoDepth,$minHomoDepthOneAlt,$minHeteroDepth,$AltAlleleStrength,$AlleleFreqProp,$CallRate,$minAvgDepth,$maxAvgDepth);

my $H = "\n###########################################################\n"
	."GBS-SNP-CROP: GBS SNP Calling Reference Optional Pipeline\n"
	."###########################################################\n"
	."Version: 4.0\n"
	."Step 7: Filter variants and call genotypes\n"
	."\nA detailed description can be found at https://github.com/halelab/GBS-SNP-CROP\n\n"
	."Usage: perl GBS-SNP-CROP-7.pl [options]\n\n"
	."Options:\n"
	."-in: Discovery master matrix input file. The output from step 6. File. Default: GSC.Summary.txt\n"
	."-out: Genotyping matrix output file. File. Default: GSC.GenoMatrix.txt\n"
	."-p: snp or indel: polymorphism type. SNPs only (snp) or SNPs + indels (indel). String. Required. Should be the same used on step 6.\n"
	."-mnHoDepth0: Minimum depth required for calling a homozygote when the alternative allele depth = 0. Numeric. Default: 5 (Diploid)\n"
	."-mnHoDepth1: Minimum depth required for calling a homozygote when the alternative allele depth = 1. Numeric. Default: 20 (Diploid)\n"
	."-mnHetDepth: Minimum depth required for each allele when calling a heterozygote. Numeric. Default: 3 (Diploid)\n"
	."-altStrength: Across the population for a given putative bi-allelic variant, this alternate allele strength parameter is the minimum proportion of non-primary allele reads that are the secondary allele. Numeric. Default: 0.8\n"
	."-mnAlleleRatio: Minimum required ratio of less frequent allele depth to more frequent allele depth. Numeric. Default: 0.25\n"
	."-mnCall: Minimum acceptable proportion of genotyped individuals to retain a variant. Numeric. Default: 0.75\n"
	."-mnAvgDepth: Minimum average depth of an acceptable variant. Numeric. Default: 3\n"
	."-mxAvgDepth: Maximum average depth of an acceptable variant. Numeric. Default: 200\n\n";

if (! defined $help or $help =~ "h" or $help =~ "H")  {
	print "$H";
	goto FINAL;
}

#################################
# Setting the parameters values
#################################
$SummaryFile = 'GSC.SummaryMatrix.txt';	$output = 'GSC.GenoMatrix.txt';
$minHomoDepth = 5;			$minHomoDepthOneAlt = 20;
$minHeteroDepth = 3;			$AltAlleleStrength = 0.8;
$AlleleFreqProp = 0.25;			$CallRate = 0.75;
$minAvgDepth = 3;			$maxAvgDepth = 200;

GetOptions(
'in=s' => \$SummaryFile,                   # file
'out=s' => \$output,                       # file
'p=s' => \$type,			   # string
'mnHoDepth0=s' => \$minHomoDepth,          # numeric
'mnHoDepth1=s' => \$minHomoDepthOneAlt,    # numeric
'mnHetDepth=s' => \$minHeteroDepth,        # numeric
'altStrength=s' => \$AltAlleleStrength,    # numeric
'mnAlleleRatio=s' => \$AlleleFreqProp,     # numeric
'mnCall=s' => \$CallRate,                  # numeric
'mnAvgDepth=s' => \$minAvgDepth,           # numeric
'mxAvgDepth=s' => \$maxAvgDepth,           # numeric
) or die "$H\n";

#########################
# Starting GBS-SNP-CROP
#########################
print "\n###############################\n# GBS-SNP-CROP, Step 7, v.4.0\n###############################\n";
my $sttime = time;

# Creating a directory
my $dir = "variants"; 
unless(-e $dir, or mkdir $dir) {die "Directory $dir cannot be created.\n";}

###################################
# Genotyping both SNPs and Indels  
###################################
if ($type eq "indel") {
	print "\nGBS-SNP-CROP is filtering both SNPs and indels and calling genotypes for your population...\n";

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
			if ( $genotype[4] ne "_" ) {
				if ( exists $Var_depth{$genotype[5]} ) {
					$Var_depth{$genotype[5]} = $Var_depth{$genotype[5]} + $genotype[4];
				} else {
					$Var_depth{$genotype[5]} = $genotype[4];
				}
				if ( $genotype[4] > 0 ) {
					if ( exists $Var_cnt{$genotype[5]} ) {
						$Var_cnt{$genotype[5]} = $Var_cnt{$genotype[5]} + 1;
					} else {
						$Var_cnt{$genotype[5]} = 1;
					}
				}
			}
						
			if ( $genotype[6] ne "_" ) {
				if ( exists $Var_depth{$genotype[7]} ) {
					$Var_depth{$genotype[7]} = $Var_depth{$genotype[7]} + $genotype[6];
				} else {
					$Var_depth{$genotype[7]} = $genotype[6];
				}
				if ( $genotype[6] > 0 ) {
					if ( exists $Var_cnt{$genotype[7]} ) {
						$Var_cnt{$genotype[7]} = $Var_cnt{$genotype[7]} + 1;
					} else {
						$Var_cnt{$genotype[7]} = 1;
					}
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
		
		# Delete variants which the reference doesn't match to either primary or secondary
		if (($pop_one_var ne $ref) and ($pop_two_var ne $ref)) {
			next;
		}
		
		# Initial filters based on population allele frequency.  
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
		my @maf = ();
		my @HetDepth = ();
		my @PriHetDepth = ();

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
			if ( $genotype[4] ne "_" ) {
				$geno_hash{$genotype[5]} = $genotype[4];
			} 
			if ( $genotype[6] ne "_" ) {
				$geno_hash{$genotype[7]} = $genotype[6];
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
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count > $alt_count) and $primary_count >= 20 and $alt_count >= ($AlleleFreqProp * $primary_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($alt_count > $primary_count) and ($alt_count < 20) and ($alt_count >= $minHeteroDepth && $primary_count >= $minHeteroDepth)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @maf,  1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($alt_count > $primary_count) and $alt_count >= 20 and $primary_count >= ($AlleleFreqProp * $alt_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count == $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count == $alt_count) and $primary_count >= 20 and $alt_count >= ($AlleleFreqProp * $primary_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + $alt_count)));
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;

			# Calling homozygotes
			} elsif ( $primary_count >= $minHomoDepth && $alt_count == 0 ) {
				$homo_pri++;
				$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/0";
				$depth = $primary_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + 0)));
			} elsif ( $primary_count >= $minHomoDepthOneAlt && $alt_count == 1 ) {
				$homo_pri++;
				$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/1";
				$depth = $primary_count;
				push @maf, 1 - ((2 * $primary_count) / (2 * ($primary_count + 0)));
			} elsif ( $alt_count >= $minHomoDepthOneAlt && $primary_count == 1 ) {
				$homo_alt++;
				$row = "$row\t$pop_two_var/$pop_two_var|1/$alt_count";
				$depth = $alt_count;
				push @maf, 1 - ((2 * $alt_count) / (2 * ($alt_count + 0)));
			} elsif ( $alt_count >= $minHomoDepth && $primary_count == 0 ) {
				$homo_alt++;
				$row = "$row\t$pop_two_var/$pop_two_var|0/$alt_count";
				$depth = $alt_count;
				push @maf, 1 - ((2 * $alt_count) / (2 * ($alt_count + 0)));


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
		
		# Minor Allele Frequency
		my $sum;
		$sum += $_ for @maf;
		my $final_maf = sprintf("%.2f", ($sum / $scored_genotypes));
		
		# Identifying paralogous variants by heterozygous deviation of allele-specific 1:1 ratio (McKinney et al. 2016)
		my $N = 0;
		my $NA = 0;
		my $Z;
		$N += $_ for @HetDepth;
		$NA += $_ for @PriHetDepth;
		if ( $N == 0 and $NA == 0 ) {
			$Z = "-";
		} else {
			my $SD = sqrt ($N * 0.5 * 0.5);
			$Z = sprintf("%.2f", (($N / 2) - $NA) / $SD);
		}
		
		# Depth filter
		my $avgDepth = sprintf("%.2f", ($cumulative_depth / $scored_genotypes));
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
	
		# Variant call type
		my $VarType = "";
		if ( $pop_one_var =~ /^(\+|-)/ or $pop_two_var =~ /^(\+|-)/ ) {
			$VarType = "Indel";
		} else {
			$VarType = "SNP";
		}
		
		# Alternative allele
		my $Alt_allele;
		if ($ref eq $pop_one_var) {
			$Alt_allele = $pop_two_var;
		} elsif ($ref eq $pop_two_var) {
			$Alt_allele = $pop_one_var;
		}

		$row = join ("\t","$header","$position","$VarType","$ref","$Alt_allele","$avgDepth","$percentage_scored_genotypes","$homo_pri","$hets","$homo_alt","$final_maf","$Z$row" );
		print $DEST "$row\n";
		$lc++;
	}
	close $VAR;
	close $DEST;

	print "\nYour $output genotyping matrix was successfully created into variants directory. \n";

	print "GBS-SNP-CROP called $lc SNPs + indels in your population using the following parameters:\n"
	."Minimum read depth required to call homozygotes when the secondary allele depth = 0: $minHomoDepth\n"
	."Minimum read depth required to call homozygotes when the secondary allele depth = 1: $minHomoDepthOneAlt\n"
	."Minimum read depth for each allele equired to call heterozygotes: $minHeteroDepth\n"
	."Minimum proportion of non-primary allele reads that are the secondary allele: $AltAlleleStrength\n"
	."Minimum proportion of less frequent allele depth to more frequent allele depth: $AlleleFreqProp\n"
	."Minimum proportion of genotyped individuals to accept a SNP: $CallRate\n"
	."Minimum acceptable average read depth: $minAvgDepth\n"
	."Maximum acceptable average read depth: $maxAvgDepth\n";
		
########################
# Genotyping only SNPs  
########################
} elsif ($type eq "snp") {
	print "\nGBS-SNP-CROP is filtering only SNPs and calling genotypes for your population...\n";
	
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

		for ( my $x=3; $x < scalar(@input); $x++ ) {
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
		# Delete monomorphic locus
		if ( ($pop_one_count == 0) or ($pop_two_count == 0) ) {
			next;
		}
		# Delete locus which the reference doesn't match to either primary or secondary
		if (($pop_one_var ne $ref) and ($pop_two_var ne $ref)) {
			next;
		}
		# UC 1. Alternative allele strength: minimum proportion of non-primary allele reads that are the secondary allele
		my $alt_allele_strength = ( $pop_two_count / ($total - $pop_one_count) );
		if ( $alt_allele_strength < $AltAlleleStrength ) {
			next;
		}
		# UCC 1. Minor allele count filter - MAF
		my $alt_allele_ratio = ( $pop_two_count / ($pop_one_count + $pop_two_count) );
		if ( $alt_allele_ratio < 0.05) {
			next;
		}
		# UCC 2. 
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
		my @HetDepth = ();
		my @PriHetDepth = ();

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
			# UC 2:
			# Calling heterozygotes
			if ( ($primary_count > $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count > $alt_count) and $primary_count >= 20 and $alt_count >= ($AlleleFreqProp * $primary_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($alt_count > $primary_count) and ($alt_count < 20) and ($alt_count >= $minHeteroDepth && $primary_count >= $minHeteroDepth)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($alt_count > $primary_count) and $alt_count >= 20 and $primary_count >= ($AlleleFreqProp * $alt_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count == $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;
			} elsif ( ($primary_count == $alt_count) and $primary_count >= 20 and $alt_count >= ($AlleleFreqProp * $primary_count)) {
				$hets++;
				$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
				$depth = $primary_count + $alt_count;
				push @HetDepth, $primary_count, $alt_count;
				push @PriHetDepth, $primary_count;

			# Calling homozygotes
			} elsif ( $primary_count >= $minHomoDepth && $alt_count == 0 ) {
				$homo_pri++;
				$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/0";
				$depth = $primary_count;
			} elsif ( $primary_count >= $minHomoDepthOneAlt && $alt_count == 1 ) {
				$homo_pri++;
				$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/1";
				$depth = $primary_count;
			} elsif ( $alt_count >= $minHomoDepthOneAlt && $primary_count == 1 ) {
				$homo_alt++;
				$row = "$row\t$pop_two_var/$pop_two_var|1/$alt_count";
				$depth = $alt_count;
			} elsif ( $alt_count >= $minHomoDepth && $primary_count == 0 ) {
				$homo_alt++;
				$row = "$row\t$pop_two_var/$pop_two_var|0/$alt_count";
				$depth = $alt_count;

			# Declaring missing genotypes
			#} elsif ( $primary_count == $minHeteroDepth || $alt_count == $minHeteroDepth ) {
			#	$row = "$row\t-|$primary_count/$alt_count";  
			#} elsif ( $primary_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
			#	$row = "$row\t-|$primary_count/$alt_count";        
			#} elsif ( $alt_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
			#	$row = "$row\t-|$primary_count/$alt_count";

			} else {
				$row = "$row\t-|$primary_count/$alt_count";
			}
				$cumulative_depth = $cumulative_depth + $depth;
		}

		# Calculating population-level genotype parameters for additional filtering
		my $scored_genotypes = $homo_pri + $hets + $homo_alt;
		my $percentage_scored_genotypes = sprintf("%.2f",($scored_genotypes) / ((scalar(@input) - 3))*100);
		if ( $scored_genotypes == 0 ) {
			next;
		}
		# Minor Allele Frequency
		my $MAF = sprintf("%.3f", ((2 * $homo_alt) + $hets) / (2 * ($homo_pri + $hets + $homo_alt)) );
		
		# Identifying paralogous variants by heterozygous deviation of 
		# allele-specific 1:1 ratio (McKinney et al. 2016)
		my $N = 0;
		my $NA = 0;
		my $Z;
		$N += $_ for @HetDepth;
		$NA += $_ for @PriHetDepth;
		if ( $N == 0 and $NA == 0 ) {
			$Z = "-";
		} else {
			my $SD = sqrt ($N * 0.5 * 0.5);
			$Z = sprintf("%.3f", (($N / 2) - $NA) / $SD);
		}
		
		# Average depth filter
		my $avgDepth = sprintf("%.2f", ($cumulative_depth / $scored_genotypes));
		if( ($avgDepth < $minAvgDepth) or ($avgDepth > $maxAvgDepth) ) {
			next;
		}
	
		# UCC 3. Minimum representation of genotypes. How informative is the locus in population.
		if ( (scalar(@input) - 3) > 20 and (($homo_pri + $hets) < 3 or ($homo_alt + $hets) < 3 ) ) {
			next;
		} elsif( (scalar(@input) - 3) < 20 and (($homo_pri + $hets) < 2 or ($homo_alt + $hets) < 2 ) ) {
			next;
		}

		# UC 3. Genotype proportion filter
		if($scored_genotypes <= $CallRate * (scalar @input - 3)){
			next;
		}
		
		# Variant call type
		my $VarType = "";
		if ( $pop_one_var =~ /^(\+|-)/ or $pop_two_var =~ /^(\+|-)/ ) {
			$VarType = "Indel";
		} else {
			$VarType = "SNP";
		}
		
		# Reference/Alternative alleles
		my $Alt_allele;
		if ($ref eq $pop_one_var) {
			$Alt_allele = $pop_two_var;
		} elsif ($ref eq $pop_two_var) {
			$Alt_allele = $pop_one_var;
		}

		$row = join ("\t","$header","$position","$VarType","$ref","$Alt_allele","$avgDepth","$percentage_scored_genotypes","$homo_pri","$hets","$homo_alt","$MAF","$Z$row" );
		print $DEST "$row\n";
		$lc++;
	}
	close $VAR;
	close $DEST;

	print "\nYour $output genotyping matrix was successfully created into variants directory. \n";

	print "GBS-SNP-CROP called $lc SNPs in your population using the following parameters:\n"
	."Minimum read depth required to call homozygotes when the secondary allele depth = 0: $minHomoDepth\n"
	."Minimum read depth required to call homozygotes when the secondary allele depth = 1: $minHomoDepthOneAlt\n"
	."Minimum read depth for each allele equired to call heterozygotes: $minHeteroDepth\n"
	."Minimum proportion of non-primary allele reads that are the secondary allele: $AltAlleleStrength\n"
	."Minimum proportion of less frequent allele depth to more frequent allele depth: $AlleleFreqProp\n"
	."Minimum proportion of genotyped individuals to accept a SNP: $CallRate\n"
	."Minimum acceptable average read depth: $minAvgDepth\n"
	."Maximum acceptable average read depth: $maxAvgDepth\n";
}

system ( "mv $SummaryFile $output ./variants" );
print "\nElapsed time: ", sprintf("%.2f",((time - $sttime)/60)), " min", "\n";
print "Please cite: Melo et al. (2016) GBS-SNP-CROP: A reference-optional pipeline for\n"
."SNP discovery and plant germplasm characterization using variable length, paired-end\n"
."genotyping-by-sequencing data. BMC Bioinformatics. DOI 10.1186/s12859-016-0879-y.\n\n";

FINAL:
exit;
