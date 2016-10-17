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

	### 6. CLUSTER PREDICTIONS
	print VH_helpers->current_time()."Clustering...\n";
	make_path($params->{"output_path"}."/".CLUSTER_DIR);
	print VH_helpers->current_time()."\tGenerating query fasta file... ";
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
	print "Done.\n";

	# 6.2. Generate blast database
	print VH_helpers->current_time()."\tGenerating blast database... ";
	my $dbname = $params->{"output_path"}."/".CLUSTER_DIR."/all_db";
	`$params->{makeblastdb} -in $query_fasta_file_name -parse_seqids -dbtype nucl -out $dbname`;
	print "Done.\n";

	# 6.3. Run blastn against all virsorter/phispy predictions and curated database phage
	print VH_helpers->current_time()."\tBlasting predictions against curated database... ";
	my $blast_file_name = $params->{'output_path'}."/".CLUSTER_DIR."/blast.aln";
	`$params->{blastn} -db $dbname -query $query_fasta_file_name -num_threads 1 -outfmt '6 std qlen' -max_target_seqs 1000000 -evalue 0.0001 -out $blast_file_name`;
	print "Done.\n";

	# 6.4. Convert blast.aln to out.tbl
	print VH_helpers->current_time()."\tConverting to MCL input... ";
	my $out_tbl_name = $params->{'output_path'}."/".CLUSTER_DIR."/out.tbl";
	`$params->{blast_to_mcl} $blast_file_name > $out_tbl_name`;
	print "Done.\n";

	# 6.5. Threshold & use the correct clustering parameter
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

	# 6.6. Run mcl with blast.aln file (filter blast hits at a certain percent identity)
	print VH_helpers->current_time()."\tRunning mcl... ";
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
	print "Done.\n";
	print "\n";

	### 7. DEFINE CORE GENOME OF EACH CLUSTER
	### Adapted from a script by Patrick Degnan <pdegnan@illinois.edu>
	print VH_helpers->current_time()."Defining Core Genomes...\n";
	print VH_helpers->current_time()."\tParsing MCL Dump... ";
	my @clusters = ();
	open(my $dump_fh, '<', $dump_file_name);
	while(my $row = <$dump_fh>){
		chomp $row;
		my @cluster = split("\t",$row);
		push @clusters, \@cluster;
	}
	close($dump_fh);
	print "Done.\n";

	for(my $i=0; $i<scalar(@clusters); $i++){
		print VH_helpers->current_time()."\tCluster $i... ";
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

		my %SEQUENCE;
		for(my $j=0; $j<scalar(@{$clusters[$i]}); $j++){
			for(my $k=$j+1; $k<scalar(@{$clusters[$i]}); $k++){
				$allmatches = `awk '\$1 ~/$clusters[$i][$j]/' $cluster_blast_name | awk '\$2 ~/$clusters[$i][$k]/' `;
				# print "[$j][$k]\n$allmatches\n";
				@matches = split(/\n/,$allmatches);
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
								$SEQUENCE{$cols[1]}[$x];
							}
						} else {
							for my $x ($cols[9]..$cols[8]){
								$SEQUENCE{$cols[1]}[$x];
							}
						}
					}
				}
			}
		}

		my $cluster_core_name = $params->{"output_path"}."/".CLUSTER_DIR."/core_$i.txt";
		open(my $core_fh, '>', $cluster_core_name);

		my $total = scalar(@{$clusters[$i]});

		foreach my $seqname (@{$clusters[$i]}){
			print $core_fh "[$seqname]\n";
			my $interval_start = -1;
			my $interval_stop = -1;
			my $interspace=0;
			for (my $j = 1; $j<scalar(@{$SEQUENCE{$seqname}}); $j++){
				# Print stop position once distance threshold is exceeded
				if($interspace == $params->{"cluster_core_max_distance"}+1 and $interval_start>=0){
					$interval_start = -1;
					print $core_fh "$interval_stop\n";
				}

				$ratio = ($SEQUENCE{$seqname}[$j]+1)/$total;
				if($ratio >= $params->{'cluster_core_congruence'}){
					if($interval_start < 0){
						$interval_start = $j;
						print $core_fh "$interval_start\t";
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
			} elsif ($interspace <= $params->{"cluster_core_max_distance"} and $interval_start >= 0){
				print $core_fh "$interval_stop\n";
			}
		}
		close($core_fh);
		unlink $cluster_blast_name;
		print "Done.\n";
	}
	print "\n";
}

1;