#!/bin/perl
############################################################
#perl script that queries editing level of known sites in a BAM file

use warnings;
use strict;
require "parse_pileup_query.pl"; #NEED PARSE PILEUP LIBRARY

if (@ARGV != 3) {
	die "need to provide 3 input:Edit Site list, INDEXED BAM alignment file and output file name\n";
}
my ($inputfile, $bamfile, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2]);

#GLOBAL VARIABLES - PLEASE MODIFY THESE

my $minbasequal = 20; # MINIMUM BASE QUALITY SCORE
my $minmapqual = 255; # MINIMUM READ MAPPING QUALITY SCORE. 255 FOR UNIQUE MAPPING WITH STAR
my $sampath = "samtools"; #PATH TO THE SAMTOOLS EXECUTABLE
my $genomepath = "Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta"; #PATH TO REFERENCE GENOME

my $offset = 33; #BASE QUALITY SCORE OFFSET - 33 FOR SANGER SCALE, 64 FOR ILLUMINA SCALE

##END GLOBAL VARIABLES

my $bedtemp = join '', $outputfile, '.bed';
system("awk \'\$1\!\=\"chromosome\"\{print \$1\"\t\"\$2-1\"\t\"\$2\}\' $inputfile \> $bedtemp");
my $piletemp = join '', $outputfile, '.pileup';
system("$sampath mpileup -A -B -d 1000000 -q $minmapqual -Q $minbasequal -f $genomepath -l $bedtemp $bamfile \> $piletemp");

my %sitehash;
open (my $PILEUP, "<", $piletemp);
while(<$PILEUP>) {
	chomp;
	my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
	my $location = join '_', $chr, $position;
	my ($refnuccount, $acount, $tcount, $ccount, $gcount) = &parse_pileup($_, $minbasequal, $offset);# parse each line of pileup
	my $counts = join ',', $refnuccount, $ccount, $gcount;
	$sitehash{$location} = $counts;
}
system("rm $bedtemp");
system("rm $piletemp");

open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);
print $OUTPUT "#chrom\tposition\tgene\tstrand\tannot1\tannot2\tcoverage\teditedreads\teditlevel\n";

while (<$INPUT>) { #READ IN LIST OF KNOWN EDITED SITES AND QUERY EDITING STATUS
	chomp;
	my @fields = split;
	next if ($fields[0] eq 'chromosome');
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $location = join '_', $chr, $position;
	my ($gene, $strand, $annot1, $annot2) = ($fields[2], $fields[3],$fields[4], $fields[5]);

	if ($sitehash{$location}) { #PRINT OUT RESULT
		my ($refcount, $ccount, $gcount) = split(/\,/,$sitehash{$location});
		my ($newcov, $newmismatch) = (0,0);
		if ($strand eq '+') {
			$newmismatch = $gcount;
		} else {
			$newmismatch = $ccount;
		}
		$newcov = $refcount + $newmismatch;
		if ($newcov) {		
			my $varfreq = 0;
			$varfreq = sprintf("%.3f", $newmismatch/$newcov);
			print $OUTPUT "$fields[0]\t$fields[1]\t$newcov\t$newmismatch\t$varfreq\n";
		} else {
			print $OUTPUT "$fields[0]\t$fields[1]\t0\t0\tN/A\n";
		}
	} else {
		print $OUTPUT "$fields[0]\t$fields[1]\t0\t0\tN/A\n";
	}
}
close $INPUT;	
close $OUTPUT;
