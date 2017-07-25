#!/usr/bin/perl
use strict;
use warnings;

unless(@ARGV){
	print "This program reformats the MCL dump file.\n";
	print "Usage: MCLdump2clusters.pl <in> <out> [<prefix>]\n";
	exit;
}

my $in = $ARGV[0];
my $out = $ARGV[1];
my $prefix = "";
if(exists $ARGV[2]){
	$prefix = $ARGV[2];
}

open(IN,"<$in") or die "Cannot open $in\n";
open(OUT,">$out");

my $cluster = 0;
while(my $line = <IN>){
	chomp $line;
	my @elements = split("\t",$line);
	foreach(@elements){
		print OUT "$prefix$cluster\t$_\n";
	}
	$cluster++;
}
close(IN);
close(OUT);
