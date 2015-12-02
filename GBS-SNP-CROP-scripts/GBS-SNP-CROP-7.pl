#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 7. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum);

my $Usage = "Usage: perl GBS-SNP-CROP-7.pl -in <count matrix> -out <genotyping matrix> -mnHoDepth 
<minimum depth value required for homozygotes call> -mnHeDepth <minimum depth value required for heterozygotes call>
-mnAlleleProp <minimum treshold of less frequent to more frequent allele> -call <SNP call rate across the population threshold value>
-mnAvgDepth <minimum acceptable of SNP average depth> -mxAvgDepth <maximum acceptable of SNP average depth>\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($SummaryFile,$output,$minHomoDepth,$minHeteroDepth,$AlleleFreqProp,$CallRate,$minAvgDepth,$maxAvgDepth);

GetOptions(
'in=s' => \$SummaryFile,    			# file
'out=s' => \$output,				# file
'mnHoDepth=s' => \$minHomoDepth, 		# numeric
'mnHeDepth=s' => \$minHeteroDepth,	 	# numeric
'mnAlleleProp=s' => \$AlleleFreqProp,  		# numeric
'call=s' => \$CallRate,				# numeric
'mnAvgDepth=s' => \$minAvgDepth, 		# numeric
'mxAvgDepth=s' => \$maxAvgDepth, 		# numeric
) or die "$Usage\n$Manual\n";

print "\nGBS-SNP-CROP is filtering variants and calling genotypes for your population...\n";


open (VAR, "$SummaryFile") || die "cant load file";
open (DEST, ">$output") || die "cant load file";

my $lc = 0;

while (<VAR>) {
	chomp;
	my @input = split("\t", $_);
	my $header = $input[0];
	my $position = $input[1];
	my $ref = $input[2];
	
# Parse mpileup count summary file in two hashes
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
#		if ( $genotype[4] ne "_" ) {
#			if ( exists $Var_depth{$genotype[5]} ) {
#				$Var_depth{$genotype[5]} = $Var_depth{$genotype[5]} + $genotype[4];
#			} else {
#				$Var_depth{$genotype[5]} = $genotype[4];
#			}
#			if ( $genotype[4] > 0 ) {
#				if ( exists $Var_cnt{$genotype[5]} ) {
#					$Var_cnt{$genotype[5]} = $Var_cnt{$genotype[5]} + 1;
#				} else {
#					$Var_cnt{$genotype[5]} = 1;
#				}
#			}
#		}
#						
#		if ( $genotype[6] ne "_" ) {
#			if ( exists $Var_depth{$genotype[7]} ) {
#				$Var_depth{$genotype[7]} = $Var_depth{$genotype[7]} + $genotype[6];
#			} else {
#				$Var_depth{$genotype[7]} = $genotype[6];
#			}
#			if ( $genotype[6] > 0 ) {
#				if ( exists $Var_cnt{$genotype[7]} ) {
#					$Var_cnt{$genotype[7]} = $Var_cnt{$genotype[7]} + 1;
#				} else {
#					$Var_cnt{$genotype[7]} = 1;
#				}
#			}
#		}
	}
	
	# Estimating the total population depth and skip NA lines	
	my $total;
	while ( my ($key, $val) = each %Var_depth ) {
    	$total += sum $val;
	}
	if ( $total == 0) {
		next;
	}
	
	# Create the @rank array and sort Variant depth hash		
	my @rank = ();
	
	foreach my $variant (sort { $Var_depth{$b} <=> $Var_depth{$a} } keys %Var_depth) {
    	push @rank, $variant;
    }

	my $pop_one_var = $rank[0];				my $pop_two_var = $rank[1];
	my $pop_one_count = $Var_depth{$pop_one_var};		my $pop_two_count = $Var_depth{$pop_two_var};

	# Initial filters based on population allele frequency.  
	# The secondary allele ratio, strength and representativeness are ask
	my $alt_allele_ratio = ( $pop_two_count / ($pop_one_count + $pop_two_count) );
	if ( $alt_allele_ratio < 0.05) {
		next;
	}

	if ( $pop_two_count == 0 ) {
		next;
	}
	
	my $alt_allele_strength = ( $pop_two_count / ($total - $pop_one_count) );
	if ( $alt_allele_strength < 0.9 ) {
		next;
	}
	
	my $N1 = $Var_cnt{$pop_one_var};
	my $N2 = $Var_cnt{$pop_two_var};
	if ( (scalar(@input) - 3) > 20 and ($N2 < 3) ) {
		next;
	} elsif( (scalar(@input) - 3) < 20 and ($N2 < 2) ) {
		next;
	}
	
	# Define some variables to use after for loop (applying the genotyping criteria)	
	my $row = "";
	my $homo_pri = 0;
	my $hets = 0;
	my $homo_alt = 0;
	my $cumulative_depth = 0;	

	# Initialize for loop across the entire population for count primary and secondary alleles for 
	# each genotype separately. A genotype specific hash will be use to count primary and 
	# secondary alleles. The genotyping criteria will be applied ... 
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
#		if ( $genotype[4] ne "_" ) {
#			$geno_hash{$genotype[5]} = $genotype[4];
#		} 
#		if ( $genotype[6] ne "_" ) {
#			$geno_hash{$genotype[7]} = $genotype[6];
#		} 
		
		# Estimating the genotype specific primary and secondary alleles 
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

		# Genotyping criteria. Call Heterozygotes		   
		if ( ($primary_count > $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($primary_count > $alt_count) and $primary_count >= 20 && $alt_count >= ($AlleleFreqProp * $primary_count)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";						
			$depth = $primary_count + $alt_count;
		} elsif ( ($alt_count > $primary_count) and ($alt_count < 20) and ($alt_count >= $minHeteroDepth && $primary_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($alt_count > $primary_count) and $alt_count >= 20 && $primary_count >= ($AlleleFreqProp * $alt_count)) {
			$hets++;
			$row = "$row\t$pop_one_var/$pop_two_var|$primary_count/$alt_count";						
			$depth = $primary_count + $alt_count;	
		
		# Genotyping criteria. Call Homozygotes        
        	} elsif ( $primary_count >= $minHomoDepth && $alt_count == 0 ) {
			$homo_pri++;
			$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/0";
			$depth = $primary_count;
		} elsif ( $primary_count >= 50 && $alt_count == 1 ) {
			$homo_pri++;
			$row = "$row\t$pop_one_var/$pop_one_var|$primary_count/1";
			$depth = $primary_count;
		} elsif ( $alt_count >= $minHomoDepth && $primary_count == 0 ) {
			$homo_alt++;
			$row = "$row\t$pop_two_var/$pop_two_var|0/$alt_count";
			$depth = $alt_count;
	
		# Genotyping criteria. Missing data
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

	# Estimate some parameters like % of scored genotypes and applying some population level 
	# filter like variant loci depth, genotyping call and control the sequencing error (loci 
	# representation across the population)
	my $scored_genotypes = $homo_pri + $hets + $homo_alt;
	my $percentage_scored_genotypes = ($scored_genotypes) / ((scalar(@input) - 3))*100;
	
	if ( $scored_genotypes == 0 ) {
		next;
	}
		
	my $avgDepth = $cumulative_depth / $scored_genotypes;
			
	# Depth filter
	if( ($avgDepth < $minAvgDepth) or ($avgDepth > $maxAvgDepth) ) {
		next;
	}	
	
	# Control sequencing error - Representation
	if ( (scalar(@input) - 3) > 20 and (($homo_pri + $hets) < 3 or ($homo_alt + $hets) < 3 ) ) {
		next;
	} elsif( (scalar(@input) - 3) < 20 and (($homo_pri + $hets) < 2 or ($homo_alt + $hets) < 2 ) ) {
		next;
	}

	# Genotype call control
	if($scored_genotypes <= $CallRate * (scalar @input - 3)){	
		next;
	}

	# Create a new column about variant type.
	# if some Indel pattern (+|-) appear on pop_one or pop_two base, this loci is an Indel call
	my $VarType = "";
	if ( $pop_one_var =~ /^(\+|-)/ or $pop_two_var =~ /^(\+|-)/ ) {
		$VarType = "Indel";
	} else {
		$VarType = "SNP";
	}
	
	$row = join ("\t","$header","$position","$VarType","$ref","$avgDepth","$pop_one_var","$pop_two_var","$percentage_scored_genotypes","$homo_pri","$hets","$homo_alt$row" );
	print DEST "$row\n";
	$lc++;
	
}

close DEST;

print "Your $output genotyping SNP matrix was successfully created. \n\n";

print "GBS-SNP-CROP called $lc SNPs in your population using the follow parameters:\n
Minimum read depth required to call homozygotes alleles: $minHomoDepth
Minimum read depth required to call heterozygotes alleles: $minHeteroDepth
Minimum proportion of less frequent to more frequent alleles: $AlleleFreqProp
Minimum threshold to call some variant rate across the population: $CallRate
Minimum acceptable average read depth: $minAvgDepth
Maximum acceptable average read depth: $maxAvgDepth\n";

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;
