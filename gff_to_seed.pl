#!/usr/bin/env perl
#########################################################################################
# Date: 2 February 2016
# Purpose: convert nucleotide FASTA and GFF inputs to SEED-format
# Author: Danielle Campbell
#########################################################################################
use strict;
use warnings;
#use File::Copy;

unless(@ARGV){
	print "Usage: gff_to_seed.pl <input_gff> <input_nt_fasta> <output_directory>\n";
	exit;
}

my $out = "_SEED_" . $ARGV[2];
mkdir $out or die "Cannot make directory $out!\n";
my $gff = $ARGV[0];

my $featdir = $out . "/Features";
mkdir $featdir;
	
my $afunct = $out . "/assigned_functions";
open(AFUNCT,">$afunct") or die "Cannot open file $afunct\n";
	
my $pegdir = $featdir . "/peg";
mkdir $pegdir;
	
my $pegtbl = $pegdir . "/tbl";
open(PEG,">$pegtbl") or die "Cannot open file $pegtbl\n";
	
my $rnadir = $featdir . "/rna";
mkdir $rnadir;
	
my $rnatbl = $rnadir . "/tbl";
open(RNA,">$rnatbl") or die "Cannot open file $rnatbl\n";
		
open(GFF,"<$gff") or die "Cannot open file $gff\n";
	
my $prodigal = 0;	
while(my $line=<GFF>){
	chomp $line;
	if($line =~ /^\#/){
		if($line =~ /Prodigal/){
			$prodigal = 1;
		}
	}else{
		my @gfftokens = split("\t",$line);
		my $contig = $gfftokens[0];
		my $gene;
		my $type = $gfftokens[2];
		my $product;
		if($prodigal == 1){
			$type = "PEG";
			$product = "hypothetical protein";
		}
		my $start = $gfftokens[3];
		my $stop = $gfftokens[4];
		my $strand = $gfftokens[6];
		my @infotokens = split(";",$gfftokens[8]);
		if($prodigal == 1){
			$infotokens[0] =~ s/^ID=//;
			$gene = $infotokens[0];
		}else{
			if($type eq "gene"){
				foreach my $token(@infotokens){			
					if($token =~ /locus_tag/){
						$token =~ s/locus_tag//;
						$gene = $token;
					}elsif($token =~ /^product/){
						$token =~ s/^product//;
						$product = $token;
						if($token =~ /RNA/){
							$type = "RNA";
						}else{
							$type = "PEG";
						}
					}	
				}
			}
		}
		print AFUNCT "$gene\t$product\n";
		if($type eq "PEG"){
			print PEG "$gene\t$contig" . "_" . "$start" . "_" . "$stop\n"
		}elsif($type eq "RNA"){
			print RNA "$gene\t$contig" . "_" . "$start" . "_" . "$stop\n"
		}
	}
}

close(AFUNCT);
close(PEG);
close(RNA);
close(GFF);

my $fna = $ARGV[1];
my $contigs = $out . "/contigs";
open(FNA,"<$fna") or die "Cannot open file $fna!\n";
open(CONTIGS,">$contigs");

while(my $l=<FNA>){
	chomp $l;
	if($l =~ /^\>/){
		print CONTIGS "$l\n";
	}else{
		$l=~s/(.{60})/$1\n/g;
		print CONTIGS "$l";
		if(length($l)<60){
			print CONTIGS "\n";
		}		
	}
}

close(FNA);
close(CONTIGS);
	
		#copy($fasta,$destination) or die "Cannot copy $fasta!\n";
		
	
