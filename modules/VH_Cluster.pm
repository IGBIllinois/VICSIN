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

	print VH_helpers->current_time()."Clustering...\n";
	make_path($params->{"output_path"}."/".CLUSTER_DIR);
	print VH_helpers->current_time()."\tGenerating query fasta file... ";
	# 1. Generate fasta file with predictions
	my $query_fasta_file_name = $params->{"output_path"}."/".CLUSTER_DIR."/all.fasta" or die "Could not truncate query file.";
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
	# Add known_viral_types
	if($params->{'known_viral_types'} ne "" and -f $params->{'known_viral_types'}){
		open(my $knowntypesfh, '<', $params->{'known_viral_types'});
		while(my $row= <$knowntypesfh>){
			print $queryfh $row;
		}
	}
	print "Done.\n";

	# 2. Generate blast database
	print VH_helpers->current_time()."\tGenerating blast database... ";
	my $dbname = $params->{"output_path"}."/".CLUSTER_DIR."/all_db";
	`$params->{makeblastdb} -in $query_fasta_file_name -parse_seqids -dbtype nucl -out $dbname`;
	print "Done.\n";

	# 3. Run blastn against all virsorter/phispy predictions and curated database phage
	print VH_helpers->current_time()."\tBlasting predictions against curated database... ";
	my $blast_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.aln";
	`$params->{blastn} -db $dbname -query $query_fasta_file_name -num_threads 1 -outfmt '6 std qlen' -max_target_seqs 1000000 -evalue 0.0001 -out $blast_file_name`;
	print "Done.\n";

	# 4. Convert blast.aln to out.tbl
	print VH_helpers->current_time()."\tConverting to MCL input... ";
	my $out_tbl_name = $params->{'output_path'}."/".CLUSTER_DIR."/out.tbl";
	`$params->{blast_to_mcl} $blast_file_name > $out_tbl_name`;
	print "Done.\n";

	# 5. Threshold & use the correct clustering parameter
	print VH_helpers->current_time()."\tApplying clustering parameter... ";
	my $mcl_in_name = $params->{'output_path'}."/".CLUSTER_DIR."/mcl_in.abc";
	if($params->{'clustering_parameter'} eq "total_length_aligned"){
		`awk '{if (\$4 >= $params->{cluster_min_length}) print \$1 "\t" \$2 "\t" \$4}' $out_tbl_name > $mcl_in_name`;
	} elsif ($params->{'clustering_parameter'} eq "total_bit_score") {
		`awk '{if (\$5 >= $params->{cluster_min_bit_score}) print \$1 "\t" \$2 "\t" \$5}' $out_tbl_name > $mcl_in_name`;
	} else { # Percent_length_allowed is default
		`awk '{if (\$3 >= $params->{cluster_min_perc_length}) print \$1 "\t" \$2 "\t" \$3}' $out_tbl_name > $mcl_in_name`;
	}
	print "Done.\n";

	# 6. Run mcl with blast.aln file (filter blast hits at a certain percent identity)
	print VH_helpers->current_time()."\tRunning mcl... ";
	my $mci_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.mci";
	my $tab_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.tab";
	my $cluster_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/out.blast.mci";
	my $dump_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/dump.blast.mci";
	my $reformat_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/dump.reformat.blast.mci";
	`$params->{mcxload} -abc $mcl_in_name --stream-mirror -o $mci_file_name -write-tab $tab_file_name > /dev/null 2>&1`;
	`$params->{mcl} $mci_file_name -I $params->{mcl_inflation} -o $cluster_file_name -q x -V all`;
	`$params->{mcxdump} -icl $cluster_file_name -tabr $tab_file_name -o $dump_file_name > /dev/null 2>&1`;
	# 7. Reformat MCL dump file
	`$params->{mcldump2clusters} $dump_file_name $reformat_file_name`;
	print "Done.\n";
}

1;