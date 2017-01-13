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

	my %clusternames;

	### 6. CLUSTER PREDICTIONS
	VH_helpers->log($params,"Clustering...");
	make_path($params->{"output_path"}."/".CLUSTER_DIR);
	VH_helpers->log($params,"\tGenerating query fasta file... ",1);
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
								$clusternames{$prefix.'-'.$seqname.'-'.$predCount} = {'prefix'=>$prefix,'sequence'=>$seqname,'bin'=>$bin,'index'=>$pindex};
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
	VH_helpers->log($params,"\tGenerating blast database... ",1);
	my $dbname = $params->{"output_path"}."/".CLUSTER_DIR."/all_db";
	`$params->{makeblastdb} -in $query_fasta_file_name -parse_seqids -dbtype nucl -out $dbname`;

	# 6.3. Run blastn against all virsorter/phispy predictions and curated database phage
	VH_helpers->log($params,"\tBlasting predictions against curated database... ",1);
	my $blast_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.aln";
	`$params->{blastn} -db $dbname -query $query_fasta_file_name -num_threads 1 -outfmt '6 std qlen' -max_target_seqs 1000000 -evalue 0.0001 -out $blast_file_name`;

	# 6.4. Convert blast.aln to out.tbl
	VH_helpers->log($params,"\tConverting to MCL input... ",1);
	my $out_tbl_name = $params->{'output_path'}."/".CLUSTER_DIR."/out.tbl";
	`$params->{blast_to_mcl} $blast_file_name > $out_tbl_name`;

	# 6.5. Threshold & use the correct clustering parameter
	VH_helpers->log($params,"\tApplying clustering parameter... ",1);
	my $mcl_in_name = $params->{'output_path'}."/".CLUSTER_DIR."/mcl_in.abc";
	if($params->{'clustering_parameter'} eq "total_length_aligned"){
		`awk '{if (\$4 >= $params->{cluster_min_length}) print \$1 "\t" \$2 "\t" \$4}' $out_tbl_name > $mcl_in_name`;
	} elsif ($params->{'clustering_parameter'} eq "total_bit_score") {
		`awk '{if (\$5 >= $params->{cluster_min_bit_score}) print \$1 "\t" \$2 "\t" \$5}' $out_tbl_name > $mcl_in_name`;
	} else { # Percent_length_allowed is default
		`awk '{if (\$3 >= $params->{cluster_min_perc_length}) print \$1 "\t" \$2 "\t" \$3}' $out_tbl_name > $mcl_in_name`;
	}

	# 6.6. Run mcl with blast.aln file (filter blast hits at a certain percent identity)
	VH_helpers->log($params,"\tRunning mcl... ",1);
	my $mci_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.mci";
	my $tab_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/blast.tab";
	my $cluster_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/out.blast.mci";
	my $dump_file_name = 		$params->{'output_path'}."/".CLUSTER_DIR."/dump.blast.mci";
	my $reformat_file_name = 	$params->{'output_path'}."/".CLUSTER_DIR."/dump.reformat.blast.mci";
	`$params->{mcxload} -abc $mcl_in_name --stream-mirror -o $mci_file_name -write-tab $tab_file_name > /dev/null 2>&1`;
	`$params->{mcl} $mci_file_name -I $params->{mcl_inflation} -o $cluster_file_name -q x -V all`;
	`$params->{mcxdump} -icl $cluster_file_name -tabr $tab_file_name -o $dump_file_name > /dev/null 2>&1`;
	# 6.7. Reformat MCL dump file
	`$params->{mcldump2clusters} $dump_file_name $reformat_file_name`;
	print "\n";

	### 7. DEFINE CORE GENOME OF EACH CLUSTER
	### Adapted from a script by Patrick Degnan <pdegnan@illinois.edu>
	VH_helpers->log($params,"Defining Core Genomes...");
	VH_helpers->log($params,"\tParsing MCL Dump... ",1);
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
		VH_helpers->log($params,"\tCluster $i/".(scalar(@clusters)-1)."...",1);
		# Create blast report containing only the results from the cluster
		VH_helpers->log($params,"\t\tCreating cluster blast report... ",2);
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

		VH_helpers->log($params,"\t\tBuilding histogram... ",2);
		my %SEQUENCE;
		for(my $j=0; $j<scalar(@{$clusters[$i]}); $j++){
			for(my $k=$j+1; $k<scalar(@{$clusters[$i]}); $k++){
				open($cluster_blast_fh, '<', $cluster_blast_name);
				my @matches;
				while(my $cluster_blast_line = <$cluster_blast_fh>){
					chomp $cluster_blast_line;
					if(index($cluster_blast_line, $clusters[$i][$j]) >= 0 and 
						index($cluster_blast_line, $clusters[$i][$k]) >= 0){
						push @matches, $cluster_blast_line;
					}
				}
				# $allmatches = `awk '\$1 ~/$clusters[$i][$j]/' $cluster_blast_name | awk '\$2 ~/$clusters[$i][$k]/' `;
				# @matches = split(/\n/,$allmatches);
				foreach my $m (@matches){
					@cols = split(/\t/,$m);
					if($cols[2] >= $params->{'percent_id_min_core'}){
						# query [0] [6] [7]
						foreach my $x ($cols[6]..$cols[7]){
							$SEQUENCE{$cols[0]}[$x]++;
						}
						# subject [1] [8] [9]
						if($cols[8] < $cols[9]){
							for my $x ($cols[8]..$cols[9]){
								$SEQUENCE{$cols[1]}[$x]++;
							}
						} else {
							for my $x ($cols[9]..$cols[8]){
								$SEQUENCE{$cols[1]}[$x]++;
							}
						}
					}
				}
			}
		}

		my $cluster_core_name = $params->{"output_path"}."/".CLUSTER_DIR."/core_$i.txt";
		open(my $core_fh, '>', $cluster_core_name);

		my $total = scalar(@{$clusters[$i]});

		VH_helpers->log($params,"\t\tWriting core file... ",2);
		my $currentseq = 0;
		foreach my $seqname (@{$clusters[$i]}){
			print $core_fh "[$seqname]\n";
			$cluster_core[$i][$currentseq] = $clusternames{$seqname};
			$cluster_core[$i][$currentseq]{'name'} = $seqname;
			$cluster_core[$i][$currentseq]{'core'} = [];
			my $interval_start = -1;
			my $interval_stop = -1;
			my $interspace=0;
			my $currentcore = 0;
			for (my $j = 1; $j<scalar(@{$SEQUENCE{$seqname}}); $j++){
				# Print stop position once distance threshold is exceeded
				if($interspace == $params->{"cluster_core_max_distance"}+1 and $interval_start>=0){
					$interval_start = -1;
					print $core_fh "$interval_stop\n";
					$cluster_core[$i][$currentseq]{'core'}[$currentcore]{'stop'} = $interval_stop;
					$currentcore++;
				}

				$ratio = ($SEQUENCE{$seqname}[$j]+1)/$total;
				if($ratio >= $params->{'cluster_core_congruence'}){
					if($interval_start < 0){
						$interval_start = $j;
						print $core_fh "$interval_start\t";
						$cluster_core[$i][$currentseq]{'core'}[$currentcore]{'start'} = $interval_start;
					}
					$interspace = 0;
				} else {
					if($interspace==0){
						$interval_stop = $j-1;
					}
					$interspace++;
				}
			}
			# Print last stop position if not already done
			if($interspace == 0 and $interval_start >= 0){
				print $core_fh scalar(@{$SEQUENCE{$seqname}})."\n";
				$cluster_core[$i][$currentseq]{'core'}[$currentcore]{'stop'} = scalar(@{$SEQUENCE{$seqname}});
			} elsif ($interspace <= $params->{"cluster_core_max_distance"}+1 and $interval_start >= 0){
				print $core_fh "$interval_stop\n";
				$cluster_core[$i][$currentseq]{'core'}[$currentcore]{'stop'} = $interval_stop;
			}
			$currentseq++;
		}
		close($core_fh);
		# unlink $cluster_blast_name;
	}
	print "\n";
	return \@cluster_core;
}

1;