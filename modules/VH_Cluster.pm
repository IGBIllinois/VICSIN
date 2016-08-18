#!/usr/bin/perl

package VH_Cluster;

use File::Path qw(make_path);
use VH_helpers;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant CLUSTER_DIR => "Cluster_Output";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;
	my $predictions = shift;

	print VH_helpers->current_time()."Starting clustering...\n";
	make_path($params->{"output_path"}."/".CLUSTER_DIR);
	print VH_helpers->current_time()."\tGenerating query fasta file... ";
	# Generate fasta file with predictions
	my $query_fasta_file_name = $params->{"output_path"}."/".CLUSTER_DIR."/query.fna" or die "Could not truncate query file.";
	open(my $queryfh, '>', $query_fasta_file_name);
	foreach my $prefix (@{$prefixes}){
		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$prefix.fna";
		open(my $fastafh, '<', $fasta_file_name) or die "Could not open fasta file.";
		my $seqname = "";
		my $sequence = "";
		
		while (my $row = <$fastafh>) {
			chomp $row;
			if(substr($row,0,1) eq '>'){
				if($seqname ne ""){
					# Grab predictions from sequence
					my $predCount = 1;
					for (my $bin = 0; $bin < 4; $bin++) {
						foreach my $prediction (@{$predictions->{$prefix}[$bin]}){
							# TODO check if virsorter/phispy?
							if ($seqname eq $prediction->{'sequence'}){
								print $queryfh '>'.$prefix.'-'.$seqname.'-'.$predCount."\n";
								print $queryfh substr($sequence,$prediction->{'start'}-1,$prediction->{'end'}-$prediction->{'start'})."\n";
								$predCount++;
							}
						}
					}
				}
				$seqname = substr($row,1);
				$sequence = "";
			} else {
				if($seqname ne ""){
					$sequence .= $row;
				}
			}
		}
	}
	print "Done.\n";
	# Run blastn against all virsorter/phispy predictions and curated database phage
	print VH_helpers->current_time()."\tBlasting predictions against curated database... ";
	my $blast_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.aln";
	# TODO this shouldnt be known viral types
	`$params->{blastn} -db $params->{known_viral_types} -query $query_fasta_file_name -num_threads 1 -outfmt 6 -max_target_seqs 1000000 -evalue 0.0001 -out $blast_file_name`;
	print "Done.\n";
	# Run mcl with blast.aln file (filter blast hits at a certain percent identity)
	print VH_helpers->current_time()."\tRunning mcl... ";
	my $abc_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.abc";
	my $mci_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.mci";
	my $tab_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.tab";
	my $cluster_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/out.blast.mci";
	my $dump_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/dump.blast.mci";
	`$params->{cut} -f 1,2,12 $blast_file_name > $abc_file_name`;
	`$params->{mcxload} -abc $abc_file_name --stream-mirror -o $mci_file_name -write-tab $tab_file_name > /dev/null 2>&1`;
	`$params->{mcl} $mci_file_name -I $params->{mcl_inflation} -o $cluster_file_name -q x -V all`;
	`$params->{mcxdump} -icl $cluster_file_name -tabr $tab_file_name -o $dump_file_name > /dev/null 2>&1`;
	print "Done.\n";

	print VH_helpers->current_time()."\tParsing mcl results... ";
	open(my $dump_fh, '<', $dump_file_name);
	while($dump_line = <$dump_fh>){
		chomp $dump_line;
		my @segments = split "\t", $dump_line;

		
	}
	print "Done.\n";
}

1;