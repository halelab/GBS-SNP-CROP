#!/usr/bin/perl

##########################################################################################
# GBS-SNP-CROP, Step 6. For description, please see Melo et al. (2015) DOI XXX
##########################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $Usage = "Usage: perl GBS-SNP-CROP-7.pl -in <input file> -out <output file> -mnHoDepth 
<minimum depth value required for homozygotes call> -mnHeDepth <minimum depth value required for heterozygotes call>
-mnAlleleProp <minimum treshold of less frequent to more frequent allele> -call <SNP call rate across the population threshold value>
-mnAvgDepth <minimum acceptable of SNP average depth> -mxAvgDepth <maximum acceptable of SNP average depth>\n";
my $Manual = "Please see Additional File 2 (User Manual) from Melo et al. (2015) BMC Bioinformatics. DOI XXX\n"; 

my ($SummaryFile,$output,$minHomoDepth,$minHeteroDepth,$AlleleFreqProp,$CallRate,$minAvgDepth,$maxAvgDepth);

GetOptions(
'in=s' => \$SummaryFile,    			# file
'out=s' => \$output,					# file
'mnHoDepth=s' => \$minHomoDepth, 		# numeric
'mnHeDepth=s' => \$minHeteroDepth,	 	# numeric
'mnAlleleProp=s' => \$AlleleFreqProp,  	# numeric
'call=s' => \$CallRate,					# numeric
'mnAvgDepth=s' => \$minAvgDepth, 		# numeric
'mxAvgDepth=s' => \$maxAvgDepth, 		# numeric
) or die "$Usage\n$Manual\n";

print "\nGBS-SNP-CROP is filtering SNPs and calling genotypes for your population...\n";


open (SNP, "$SummaryFile") || die "cant load file";
open (DEST, ">$output") || die "cant load file";

my $lc = 0;

while (<SNP>)
{
	chomp;
	my @input = split("\t", $_);
	my $header = $input[0];
	my $position = $input[1];
	my $ref = $input[2];

	my $APop=0;
	my $CPop=0;
	my $GPop=0;
	my $TPop=0;

	my $NA=0;
	my $NC=0;
	my $NG=0;
	my $NT=0;

	for ( my $x=3; $x < scalar(@input); $x++ ){
	
		my $geno=$input[$x];
		my @genotype = split(",", $geno);

		if ( $genotype[0] ne "_" ) {
			$APop = $APop + $genotype[0];
			if ( $genotype[0] > 0 ) {
				$NA++;
			}
		}	

		if ( $genotype[1] ne "_" ) {
			$CPop = $CPop + $genotype[1];
			if ( $genotype[1] > 0 ) {
				$NC++;
			}
		}
		
		if ( $genotype[2] ne "_" ) {
			$GPop = $GPop + $genotype[2];
			if ( $genotype[2] > 0 ) {
				$NG++;
			}
		}
		
		if ( $genotype[3] ne "_" ) {
			$TPop = $TPop + $genotype[3];
			if ( $genotype[3] > 0 ) {
				$NT++;
			}
		}		
		
	}
	
	if ( $APop == 0 && $CPop == 0 && $GPop == 0 && $TPop == 0 ) {
		next;
	}
	
	my %rank = (	'A' => $APop,
					'C' => $CPop,
					'G' => $GPop,
					'T' => $TPop    );
	
	my @rank = ();
	
	foreach my $base (sort { $rank{$b} <=> $rank{$a} } keys %rank) {
    	push @rank, $base;
    }

	my $pop_one_base = $rank[0];
	my $pop_two_base = $rank[1];
	my $pop_three_base = $rank[2];
	my $pop_four_base = $rank[3];
	
	my $pop_one_count = $rank{$pop_one_base};
	my $pop_two_count = $rank{$pop_two_base};
	my $pop_three_count = $rank{$pop_three_base};
	my $pop_four_count = $rank{$pop_four_base};

	my $alt_allele_ratio = $pop_two_count / ($pop_one_count + $pop_two_count );
	
	if ( $alt_allele_ratio < 0.05) {
		next;
	}
	
	if ( $pop_two_count == 0 ) {
		next;
	}
	
	my $alt_allele_strength = ( $pop_two_count / ( $pop_two_count + $pop_three_count + $pop_four_count ) );
	
	if ( $alt_allele_strength < 0.9 ) {
		next;
	}
	
	my %geno_count = (	'A' => $NA,
						'C' => $NC,
						'G' => $NG,
						'T' => $NT     );
	
	my $N1 = $geno_count{$pop_one_base};
	my $N2 = $geno_count{$pop_two_base};
	
	if ( $N2 < 2 ) {
		next;
	}
		
	my $row = "";	

	my %base_index = (	'A' => 0,
						'C' => 1,
						'G' => 2,
						'T' => 3    );
	
	my $homo_pri = 0;
	my $hets = 0;
	my $homo_alt = 0;
	my $cumulative_depth = 0;
	
	for ( my $x=3; $x < scalar(@input); $x++ ){
	
		my $geno=$input[$x];
		my @genotype = split(",",$geno);
		
		my $primary_count = $genotype[$base_index{$pop_one_base}];
		my $alt_count = $genotype[$base_index{$pop_two_base}];
		
		if ( $primary_count eq "_" ) {
			$primary_count = 0;
		}

		if ( $alt_count eq "_" ) {
			$alt_count = 0;
		}
		
		my $depth = 0;
		   
		if ( ($primary_count > $alt_count) and ($primary_count < 20) and ($primary_count >= $minHeteroDepth && $alt_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_base/$pop_two_base|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($primary_count > $alt_count) and $primary_count >= 20 && $alt_count >= ($AlleleFreqProp * $primary_count)) {
			$hets++;
			$row = "$row\t$pop_one_base/$pop_two_base|$primary_count/$alt_count";						
			$depth = $primary_count + $alt_count;
			
		} elsif ( ($alt_count > $primary_count) and ($alt_count < 20) and ($alt_count >= $minHeteroDepth && $primary_count >= $minHeteroDepth)) {
			$hets++;
			$row = "$row\t$pop_one_base/$pop_two_base|$primary_count/$alt_count";
			$depth = $primary_count + $alt_count;
		} elsif ( ($alt_count > $primary_count) and $alt_count >= 20 && $primary_count >= ($AlleleFreqProp * $alt_count)) {
			$hets++;
			$row = "$row\t$pop_one_base/$pop_two_base|$primary_count/$alt_count";						
			$depth = $primary_count + $alt_count;			
    	} elsif ( $primary_count == $minHeteroDepth || $alt_count == $minHeteroDepth ) {
           $row = "$row\t-|$primary_count/$alt_count";
        
        } elsif ( $primary_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
           $row = "$row\t-|$primary_count/$alt_count";
        
        } elsif ( $alt_count == 0 && ($primary_count + $alt_count) < $minHomoDepth ) {
            $row = "$row\t-|$primary_count/$alt_count";
        
        } elsif ( $primary_count >= $minHomoDepth && $alt_count == 0 ) {
			$homo_pri++;
			$row = "$row\t$pop_one_base/$pop_one_base|$primary_count/0";
			$depth = $primary_count + $alt_count;
			
		} elsif ( $alt_count >= $minHomoDepth && $primary_count == 0 ) {
			$homo_alt++;
			$row = "$row\t$pop_two_base/$pop_two_base|0/$alt_count";
			$depth = $primary_count + $alt_count;
		
		} else {
			$row = "$row\t-|-";
		}
		
		$cumulative_depth = $cumulative_depth + $depth;
	}

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
	if( ($homo_pri + $hets) < 3 or ($homo_alt + $hets) < 3 ) {
		next;
	}

	# Genotype call control
	if($scored_genotypes <= $CallRate * (scalar @input - 3)){	
		next;
	}
	
	$row = "$header\t$position\t$ref\t$avgDepth\t$pop_one_base\t$pop_two_base\t$percentage_scored_genotypes\t$homo_pri\t$hets\t$homo_alt$row";
	print DEST "$row\n";
	$lc++;

}

close DEST;

sub main {
   	my $dir = "SNPsCalled"; 
   	unless(-e $dir, or mkdir $dir) {die "Directory $dir just exist.\n";}
}
main();

system ( "mv *.txt ./SNPsCalled" );

print "Your $output genotyping matrix was successfully created. Please, see 'SNPsCalled' directory.\n\n";

print "GBS-SNP-CROP called $lc SNPs in your population using the follow parameters:\n
Minimum read depth required for calling homozygotes: $minHomoDepth
Minimum read depth for each allele required for calling heterozygotes: $minHeteroDepth
Minimum proportion of less frequent to more frequent allele: $AlleleFreqProp
Minimum threshold for SNP call rate across the population: $CallRate
Minimum acceptable SNP average read depth: $minAvgDepth
Maximum acceptable SNP average read depth: $maxAvgDepth\n";

print "\n\nPlease cite: Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for
SNP discovery and plant germplasm characterization using variable length, paired-end
genotyping-by-sequencing data. BMC Bioinformatics. DOI XXX.\n\n";

exit;