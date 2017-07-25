#!/usr/bin/perl

# Prediction Clustering module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_Cluster;

use File::Path qw(make_path);
use VH_helpers;
use Tie::File;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant CLUSTER_DIR => "Cluster_Output";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;
	my $predictions = shift;

	my %clusternames;

	### 6. CLUSTER PREDICTIONS
	VH_helpers->log($params,"Clustering...");
	make_path($params->{"output_path"}."/".CLUSTER_DIR);
	VH_helpers->log($params,"\t\tGenerating query fasta file... ",2);
	# 6.1. Generate fasta file with predictions
	my $query_fasta_file_name = $params->{"output_path"}."/".CLUSTER_DIR."/all.fasta";
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
						for (my $pindex=0;$pindex<scalar(@{$predictions->{$prefix}[$bin]});$pindex++){
							# TODO check if virsorter/phispy?
							if ($seqname eq $predictions->{$prefix}[$bin][$pindex]{'sequence'}){
								print $queryfh '>'.$prefix.'-'.$seqname.'-'.$predCount."\n";
								$clusternames{$prefix.'-'.$seqname.'-'.$predCount} = {'prefix'=>$prefix,'sequence'=>$seqname,'bin'=>$bin,'index'=>$pindex, 'length'=>abs($predictions->{$prefix}[$bin][$pindex]{'start'}-$predictions->{$prefix}[$bin][$pindex]{'end'})};
								print $queryfh substr($sequence,$predictions->{$prefix}[$bin][$pindex]{'start'}-1,$predictions->{$prefix}[$bin][$pindex]{'end'}-$predictions->{$prefix}[$bin][$pindex]{'start'})."\n";
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
		close($fastafh)
	}
	#print Dumper(%clusternames);
	# Add known_viral_types
	if($params->{'known_viral_types'} ne "" and -f $params->{'known_viral_types'}){
		open(my $knowntypesfh, '<', $params->{'known_viral_types'});
		while(my $row= <$knowntypesfh>){
			print $queryfh $row;
		}
		close($knowntypesfh);
	}
	close($queryfh);

	# 6.2. Generate blast database
	VH_helpers->log($params,"\t\tGenerating blast database... ",2);
	my $dbname = $params->{"output_path"}."/".CLUSTER_DIR."/all_db";
	`$params->{makeblastdb} -in $query_fasta_file_name -parse_seqids -dbtype nucl -out $dbname`;

	# 6.3. Run blastn against all virsorter/phispy predictions and curated database phage
	VH_helpers->log($params,"\t\tBlasting predictions against curated database... ",2);
	my $blast_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.aln";
	`$params->{blastn} -db $dbname -query $query_fasta_file_name -num_threads 1 -outfmt '6 std qlen' -max_target_seqs 1000000 -evalue 0.0001 -out $blast_file_name`;

	# 6.4. Convert blast.aln to out.tbl
	VH_helpers->log($params,"\t\tConverting to MCL input... ",2);
	my $out_tbl_name = $params->{'output_path'}."/".CLUSTER_DIR."/out.tbl";
	`$params->{blast_to_mcl} $blast_file_name > $out_tbl_name`;

	# 6.5. Threshold & use the correct clustering parameter
	VH_helpers->log($params,"\t\tApplying clustering parameter... ",2);
	my $mcl_in_name = $params->{'output_path'}."/".CLUSTER_DIR."/mcl_in.abc";
	if($params->{'clustering_parameter'} eq "total_length_aligned"){
		`awk '{if (\$4 >= $params->{cluster_min_length}) print \$1 "\t" \$2 "\t" \$4}' $out_tbl_name > $mcl_in_name`;
	} elsif ($params->{'clustering_parameter'} eq "total_bit_score") {
		`awk '{if (\$5 >= $params->{cluster_min_bit_score}) print \$1 "\t" \$2 "\t" \$5}' $out_tbl_name > $mcl_in_name`;
	} else { # Percent_length_allowed is default
		`awk '{if (\$3 >= $params->{cluster_min_perc_length}) print \$1 "\t" \$2 "\t" \$3}' $out_tbl_name > $mcl_in_name`;
	}

	# 6.5.a. Filter out small predictions from mcl_in.abc
	# TODO put in known types
	VH_helpers->log($params,"\t\tFiltering out small predictions... ",2);
	my $mcl_in_large_name = $params->{'output_path'}."/".CLUSTER_DIR."/mcl_in_large.abc";
	open(my $mcl_in_fh, '<', $mcl_in_name);
	open(my $mcl_in_large_fh, '>', $mcl_in_large_name);
	while (my $row = <$mcl_in_fh>){
		chomp $row;
		my @mclrow = split "\t", $row;
		if(exists $clusternames{$mclrow[0]} and $clusternames{$mclrow[0]}{'length'}>$params->{'cluster_size_threshold'}
			and exists $clusternames{$mclrow[1]} and $clusternames{$mclrow[1]}{'length'}>$params->{'cluster_size_threshold'}){
			print $mcl_in_large_fh $row."\n";
		}
	}

	# 6.6. Run mcl with blast.aln file (filter blast hits at a certain percent identity)
	VH_helpers->log($params,"\t\tRunning mcl... ",2);
	my $mci_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.mci";
	my $tab_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.tab";
	my $cluster_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/out.blast.mci";
	my $dump_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/dump.blast.mci";
	my $reformat_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/dump.reformat.blast.mci";
	`$params->{mcxload} -abc $mcl_in_large_name --stream-mirror -o $mci_file_name -write-tab $tab_file_name > /dev/null 2>&1`;
	`$params->{mcl} $mci_file_name -I $params->{mcl_inflation} -o $cluster_file_name -q x -V all`;
	`$params->{mcxdump} -icl $cluster_file_name -tabr $tab_file_name -o $dump_file_name > /dev/null 2>&1`;
	# 6.7. Reformat MCL dump file
	`$params->{mcldump2clusters} $dump_file_name $reformat_file_name`;
	print "\n";

	#6.8 Add small predictions to clusters
	VH_helpers->log($params,"\t\tAdding small predictions to existing clusters... ",2);
	open(my $reformat_fh, '<', $reformat_file_name);
	while(my $row = <$reformat_fh>){
		chomp $row;
		my @reformat_row = split "\t", $row;
		if(exists $clusternames{$reformat_row[1]}){
			$clusternames{$reformat_row[1]}{'cluster'} = $reformat_row[0];
		}
	}
	close($reformat_fh);

	tie my @dump_file, 'Tie::File', $dump_file_name;
	tie my @reformat_file, 'Tie::File', $reformat_file_name;
	foreach my $pred (keys %clusternames){
		if($clusternames{$pred}{'length'} <= $params->{'cluster_size_threshold'}){
			seek($mcl_in_fh, 0, SEEK_SET);
			my $cluster = -1;
			while(my $row = <$mcl_in_fh>){
				chomp $row;
				my @mclrow = split "\t", $row;
				if($mclrow[0] eq $pred and exists $clusternames{$mclrow[1]} and $clusternames{$mclrow[1]}{'length'}>$params->{'cluster_size_threshold'}){
					if($cluster<0 or $cluster==$clusternames{$mclrow[1]}{'cluster'}){
						$cluster = $clusternames{$mclrow[1]}{'cluster'};
					} else {
						$cluster = -1;
						break;
					}
				}
				if($mclrow[1] eq $pred and exists $clusternames{$mclrow[0]} and $clusternames{$mclrow[0]}{'length'}>$params->{'cluster_size_threshold'}){
					if($cluster<0 or $cluster==$clusternames{$mclrow[0]}{'cluster'}){
						$cluster = $clusternames{$mclrow[0]}{'cluster'};
					} else {
						$cluster = -1;
						break;
					}
				}
			}
			if($cluster>=0){
				# Add to dump file
				$dump_file[$cluster] .= "\t".$pred;
				$clusternames{$pred}{'cluster'} = $cluster;
			}
		}
	}
	# Reformat dump file again, since it's changed
	`$params->{mcldump2clusters} $dump_file_name $reformat_file_name`;

	#6.9 Cluster remaining small predictions
	VH_helpers->log($params, "\t\tClustering remaining small predictions... ",2);
	my $mcl_in_small_name = $params->{'output_path'}."/".CLUSTER_DIR."/mcl_in_small.abc";
	seek($mcl_in_fh, 0, SEEK_SET);
	open(my $mcl_in_small_fh, '>', $mcl_in_small_name);
	while (my $row = <$mcl_in_fh>){
		chomp $row;
		my @mclrow = split "\t", $row;
		if(exists $clusternames{$mclrow[0]} and not exists $clusternames{$mclrow[0]}{'cluster'}
			and exists $clusternames{$mclrow[1]} and not exists $clusternames{$mclrow[1]}{'cluster'}){
			print $mcl_in_small_fh $row."\n";
		}
	}

	my $mci_small_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast_small.mci";
	my $tab_small_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast_small.tab";
	my $cluster_small_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/out_small.blast.mci";
	my $dump_small_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/dump_small.blast.mci";
	my $reformat_small_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/dump_small.reformat.blast.mci";
	`$params->{mcxload} -abc $mcl_in_small_name --stream-mirror -o $mci_small_file_name -write-tab $tab_small_file_name > /dev/null 2>&1`;
	`$params->{mcl} $mci_small_file_name -I $params->{mcl_inflation} -o $cluster_small_file_name -q x -V all`;
	`$params->{mcxdump} -icl $cluster_small_file_name -tabr $tab_small_file_name -o $dump_small_file_name > /dev/null 2>&1`;
	# 6.10. Reformat small MCL dump file
	`$params->{mcldump2clusters} $dump_small_file_name $reformat_small_file_name S`;

	# 6.11. Add small clusters to previous dump files
	tie my @dump_small_file, 'Tie::File', $dump_small_file_name;
	tie my @reformat_small_file, 'Tie::File', $reformat_small_file_name;
	push @dump_file, @dump_small_file;
	push @reformat_file, @reformat_small_file;

	untie @dump_file;
	untie @reformat_file;
	untie @dump_small_file;
	untie @reformat_small_file;


	### 7. DEFINE CORE GENOME OF EACH CLUSTER
	VH_helpers->log($params,"Defining Core Genomes...");
	VH_helpers->log($params,"\t\tParsing MCL Dump... ",2);
	my @clusters = ();
	open(my $dump_fh, '<', $dump_file_name);
	while(my $row = <$dump_fh>){
		chomp $row;
		my @cluster = split("\t",$row);
		push @clusters, \@cluster;
	}
	close($dump_fh);

	my @cluster_core;
	for(my $i=0; $i<scalar(@clusters); $i++){
		VH_helpers->log($params,"\t\tCreating cluster blast report $i/".(scalar(@clusters)-1)."...",2);
		# Create blast report containing only the results from the cluster
		my $cluster_blast_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast_cluster_$i.aln";
		open(my $blast_fh, '<', $blast_file_name);
		open(my $cluster_blast_fh, '>', $cluster_blast_name);
		while (my $row = <$blast_fh>){
			chomp $row;
			my @blastrow = split "\t", $row;
			for(my $j=0; $j<scalar(@{$clusters[$i]}); $j++){
				if($blastrow[1] eq $clusters[$i][$j]){
					print $cluster_blast_fh $row."\n";
				}
			}
		}
		close $blast_fh;
		close $cluster_blast_fh;

		my $cluster_core_name = $params->{"output_path"}."/".CLUSTER_DIR."/core_$i.txt";

		`$params->{core_genome} $cluster_blast_name $params->{percent_id_min_core} $params->{cluster_core_max_distance} > $cluster_core_name`;
		unlink $cluster_blast_name;
	}
	print "\n";
	return \@cluster_core;
}

1;