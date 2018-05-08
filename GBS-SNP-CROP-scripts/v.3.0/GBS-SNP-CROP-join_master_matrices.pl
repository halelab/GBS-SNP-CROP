#!/usr/bin/perl

#this script joins several master matrices into a big one.
#To avoid having to read all the data into memory, each matrix is read
#one line at a time. We take advandage of the internal sorting of the
#matrices (alphabetical for chromosome, then ascending for position).
#In the spirit of other GBS-SNP-CROP scripts, no function/OOP was used
#(even there was ample room for both).

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Zlib;
use Array::Utils qw(array_diff);

my $Usage = "Usage: perl GBS-SNP-CROP-join_master_matrices.pl\n" 
	."-l <file containing the list of master matrix files to be joined, one file per line>\n" 
	."-o <output file name>\n";

my ($infile_list,$outfile);

GetOptions(
	'l=s' => \$infile_list,	 # file with the list of files to be joined
	'o=s' => \$outfile,			 # output file
) or die "$Usage\n";

#let's just check the outfile, so that to avoid parsing the matrices
#only to discover that we cannot write
open my $OUTFILE, ">", "$outfile" or die "[ERROR] Cannot write file $outfile";

#read list of input files
open my $INFILE_LIST, '<', $infile_list or die "[ERROR] Cannot open list file $infile_list";
chomp(my @infiles = <$INFILE_LIST>);
close $INFILE_LIST;

#from here on the general approach is to use the file name (complete with
#path) as a hash key for each input matrix

#file pointers, opened
my %FPs;
foreach my $infile (@infiles){
	$FPs{$infile} = new IO::Zlib;
	if (!$FPs{$infile}->open($infile, 'rb')){
		die "[ERROR] Cannot open $infile\n";
	}
}

#filling up the buffer with first lines from each file
my %lines;
foreach my $infile (@infiles){
	$lines{$infile} = $FPs{$infile}->getline();
	chomp($lines{$infile});
}

#this is the stop flag, it goes to zero when all files are processed
my $continue = 1;

#these hashes are all indexed with matrix full filename
my %samples_num; #counter for sample numbers for each matrix
my %SNP_file;    #current SNP for each specific matrix
my %genos;       #current genotypes for each specific matrix

#this hash is indexed by available SNPs signatures, so that
#it is easy to get the minimum set with no redundancy
my %SNPs;

while ($continue){
	#extracting the SNP signature (chrom, pos, ref) and genotypes from each line
	foreach my $infile (@infiles){
		#should we process this line?
		#$lines{$infile} is undef if we already processed all $infile
		#$SNP_file{$infile} is defined if we already processed the first line from $infile
		if (!defined $lines{$infile} or defined $SNP_file{$infile}){
			next;
		}
		
		#some chomping is always welcome
		chomp($lines{$infile});
		
		#$pieces[0] = chromosome name
		#$pieces[1] = position in chromosome
		#$pieces[2] = reference allele
		#$pieces[3+] = sample readings depths (genotypes)
		my @pieces = split(/\t/, $lines{$infile});
		
		#number of samples for this file
		$samples_num{$infile} = (scalar @pieces) - 3;
		
		#signature for this SNP (first three fields)
		my $sign = join("\t", shift(@pieces), shift(@pieces), shift(@pieces));
		
		#taking notes that we found this SNP in this file
		$SNPs{$sign} = '';
		$SNP_file{$infile} = $sign;

		#genotypes for this SNP/file combo
		$genos{$infile} = join("\t", @pieces);
	}

	#finding what is the first SNP to be processed
	my @SNPs_sorted = sort {(split /\t/, $a)[0] cmp (split /\t/, $b)[0] || (split /\t/, $a)[1] <=> (split /\t/, $b)[1]} keys %SNPs;
	my $first = $SNPs_sorted[0];
	
	#now for each file, if it refers to the best SNP we use its genotype, 
	#otherwise we put the missing string
	my $outline = $first;
	foreach my $infile (@infiles){
		if (!defined $SNP_file{$infile} or $SNP_file{$infile} ne $first){
			#missing
			$outline .= "\t_,_,_,_" x $samples_num{$infile};
		}else{
			#putting in the genotype
			$outline .= "\t".$genos{$infile};
			
			#no need for this SNP anymore 
			$SNP_file{$infile} = undef;
			delete $SNPs{$first};
			
			#reading the new line for this file
			$lines{$infile} = $FPs{$infile}->getline();
		}
	}
	
	#ready to write the outline
	print $OUTFILE "$outline\n";
	
	#are we done with the files?
	$continue = 0;
	foreach my $infile (@infiles){
		if (defined $lines{$infile}){
			$continue = 1;
		}
	}
}

print "All is good under the sun\n";
select()->flush();

#file pointers, closed
foreach my $infile (@infiles){
	$FPs{$infile}->close();
}
close ($OUTFILE);


