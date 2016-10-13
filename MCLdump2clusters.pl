#!/usr/bin/perl
use strict;
use warnings;

unless(@ARGV){
	print "This program reformats the MCL dump file.\n";
	print "Usage: MCLdump2clusters.pl <in> <out>\n";
	exit;
}

my $in = $ARGV[0];
my $out = $ARGV[1];

open(IN,"<$in") or die "Cannot open $in\n";
open(OUT,">$out");

my $cluster = 0;
while(my $line = <IN>){
	chomp $line;
	my @elements = split("\t",$line);
	foreach(@elements){
		print OUT "$cluster\t$_\n";
	}
	$cluster++;
}
close(IN);
close(OUT);
