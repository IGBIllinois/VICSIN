#!/usr/bin/perl

package VH_ReBlast;

use File::Path qw(make_path);
use VH_helpers;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant REBLAST_DIR => "ReBlast_Runs";

sub run {
	# TODO redo this function to do one prefix at a time, not all together
	shift;
	my $params = shift;
	my $prefixes = shift;
	my $predictions = shift;
	my $masks = shift;
	my $contigs = shift;

	my %return_predictions;

	print VH_helpers->current_time()."Starting blastn against predictions...\n";
	make_path($params->{"output_path"}."/".REBLAST_DIR);
	foreach my $curprefix (@$prefixes){
		print VH_helpers->current_time()."\t$curprefix...\n";
		# Generate fasta file with predictions to query against
		print VH_helpers->current_time()."\t\tGenerating fasta file for re-blast... ";
		my $query_fasta_file_name = $params->{"output_path"}."/".REBLAST_DIR."/$curprefix-query.fna" or die "Could not truncate query file.";
		open(my $queryfh, '>', $query_fasta_file_name);
		my %query_lengths;
		foreach my $prefix (@{$prefixes}){
			if($curprefix ne $prefix){
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
									if ($seqname eq $prediction->{'sequence'}){
										print $queryfh '>'.$prefix.'-'.$seqname.'-'.$predCount."\n";
										print $queryfh substr($sequence,$prediction->{'start'}-1,$prediction->{'end'}-$prediction->{'start'})."\n";
										$query_lengths{$prefix.'-'.$seqname.'-'.$predCount} = $prediction->{'end'}-$prediction->{'start'}+1;
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

				# One more time to parse the final sequence
				if($seqname ne ""){
					# Grab predictions from sequence
					my $predCount = 1;
					for (my $bin = 0; $bin < 2; $bin++) {
						foreach my $prediction (@{$predictions->{$prefix}[$bin]}){
							if ($seqname eq $prediction->{'sequence'}){
								print $queryfh '>'.$prefix.'-'.$seqname.'-'.$predCount."\n";
								print $queryfh substr($sequence,$prediction->{'start'}-1,$prediction->{'end'}-$prediction->{'start'})."\n";
								$query_lengths{$prefix.'-'.$seqname.'-'.$predCount} = $prediction->{'end'}-$prediction->{'start'}+1;
								$predCount++;
							}
						}
					}
				}
			}
		}
		print "Done.\n";

		# Re-blast predictions against all prefixes
		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$curprefix.fna";
		my $query_file_name = $params->{"output_path"}."/".REBLAST_DIR."/$curprefix-query.fna";
		my $br_file_name = $params->{"output_path"}."/".REBLAST_DIR."/$curprefix.br";
		# Run blastn
		print VH_helpers->current_time()."\t\tRunning blastn... ";
		`$params->{blastn} -query $query_file_name -subject $fasta_file_name -outfmt 6 -out $br_file_name`;
		print "Done.\n";

		# Delete fasta file (it's huge)
		unlink($query_file_name);

		# Parse Results
		print VH_helpers->current_time()."\t\tParsing blast output... ";
		my @blast_predictions;
		open my $br_fh, '<', $br_file_name;
		while(my $br_line = <$br_fh>){
			chomp $br_line;
			if ($br_line =~ m/^(.+?)\s(.+?)\s(.+?)\s.+?\s.+?\s.+?\s.+?\s.+?\s(.+?)\s(.+?)\s/m){
				my $query_seq = $1;
				my $seq_name = $2;
				my $perc_identity = $3;
				my $start;
				my $end;
				if($4<$5){
					$start = $4;
					$end = $5;	
				} else {
					$start = $5;
					$end = $4;
				}
				

				push @blast_predictions, {'sequence'=>$seq_name,'methods'=>'R','start'=>$start,'end'=>$end,'query_seq'=>$query_seq,'perc_identity'=>$perc_identity};
			}
		}
		print scalar(@blast_predictions)." ";
		print "Done.\n";

		# Mask results w/ masking file
		if(exists $masks->{$curprefix}){
			print VH_helpers->current_time()."\t\tApplying mask file... ";
			foreach my $sequence ( keys %{$masks->{$curprefix}} ){
				foreach my $mask (@{$masks->{$curprefix}{$sequence}} ){
					foreach my $prediction (@blast_predictions){
						if($prediction->{'sequence'} eq $sequence){
							if($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'end'}){
								# Prediction entirely contained within mask; ignore prediction
								$prediction->{'masked'} = 1;
							} elsif ($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'start'}){
								# Prediction's left side overlaps mask; trim left side
								$prediction->{'start'} = $mask->{'end'}+1;
							} elsif ($mask->{'start'}<=$prediction->{'end'} and $mask->{'end'}>=$prediction->{'end'}){
								# Prediction's right side overlaps mask; trim right side
								$prediction->{'end'} = $mask->{'start'}-1;
							} elsif ($mask->{'start'}>$prediction->{'start'} and $mask->{'end'}<$prediction->{'end'} and not exists $prediction->{'masked'}){
								# Prediction contains a masked area within it; set to ignore prediction
								$prediction->{'masked'} = 1;
							}
						}
					}
				}
			}

			print "Done.\n";
		}


		# Mask results w/ previous predictions
		print VH_helpers->current_time()."\t\tMasking previous predictions... ";
		foreach my $prefix ( @{$prefixes} ){
			if($curprefix eq $prefix){
				for (my $bin = 0; $bin < 4; $bin++) {
					foreach my $mask (@{$predictions->{$prefix}[$bin]}){
						foreach my $prediction (@blast_predictions){
							if($prediction->{'sequence'} eq $mask->{'sequence'}){
								if($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'end'}){
									# Prediction entirely contained within mask; ignore prediction
									$prediction->{'masked'} = 1;
								} elsif ($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'start'}){
									# Prediction's left side overlaps mask; trim left side
									$prediction->{'start'} = $mask->{'end'}+1;
								} elsif ($mask->{'start'}<=$prediction->{'end'} and $mask->{'end'}>=$prediction->{'end'}){
									# Prediction's right side overlaps mask; trim right side
									$prediction->{'end'} = $mask->{'start'}-1;
								} elsif ($mask->{'start'}>$prediction->{'start'} and $mask->{'end'}<$prediction->{'end'} and not exists $prediction->{'masked'}){
									# Prediction contains a masked area within it; set to ignore prediction
									$prediction->{'masked'} = 1;
								}
							}
						}
					}
				}
			}
		}

		print "Done.\n";

		# Remove masked predictions from the list so the next step doesn't take quite so long
		my @masked_predictions;
		foreach my $prediction (@blast_predictions){
			if(not exists $prediction->{'masked'}){
				push @masked_predictions, $prediction;
			}
		}


		# Combine overlapping results
		print VH_helpers->current_time()."\t\tCombining results... ";
		my $found_overlaps;
		my @combined_predictions;
		do {
			print scalar(@masked_predictions)." ";
			@combined_predictions = ();
			$found_overlaps = 0;
			for (my $i = 0; $i < scalar(@masked_predictions); $i++) {
				my $found_match = 0;
				for (my $j = 0; $j < scalar(@combined_predictions); $j++) {
					#If results overlap
					if($found_match == 0 and $masked_predictions[$i]{'sequence'}==$combined_predictions[$j]{'sequence'} 
						and $masked_predictions[$i]{'start'}<=$combined_predictions[$j]{'end'}+$params->{'reblast_distance'} 
						and $combined_predictions[$j]{'start'}<=$masked_predictions[$i]{'end'}+$params->{'reblast_distance'}){
						$found_match = 1;
						$found_overlaps = 1;
						# Expand bounds of j to contain i
						if($masked_predictions[$i]{'start'}<$combined_predictions[$j]{'start'}){
							$combined_predictions[$j]{'start'} = $masked_predictions[$i]{'start'};
						}
						if($masked_predictions[$j]{'end'}>$masked_predictions[$i]{'end'}){
							$masked_predictions[$i]{'end'} = $masked_predictions[$j]{'end'};
						}
					}
				}
				if($found_match == 0){
					# If this prediction doesn't overlap with any others, add it to the array
					push @combined_predictions, $masked_predictions[$i];
				}
			}
			@masked_predictions = ();
			push @masked_predictions,@combined_predictions;
		} while ($found_overlaps == 1);
		print scalar(@combined_predictions)." ";
		print "Done.\n";

		# Weed out results based on length
		print VH_helpers->current_time()."\t\tExcluding based on length/% identity... ";
		foreach my $prediction (@combined_predictions){
			# TODO implement new rules based on end-ness
			if($prediction->{'start'} <= $params->{'reblast_edge_distance'} xor $prediction->{'end'}>=$contigs->{$curprefix}{$prediction->{'sequence'}}{'length'}-$params->{'reblast_edge_distance'}){
				if($prediction->{'perc_identity'} < $params->{"reblast_min_perc_id"}){
					# If hit touches either contig end (not both), min 90% ident
					$prediction->{'masked'} = 1;
				}
			} elsif ($prediction->{'start'}>$params->{'reblast_edge_distance'} and $prediction->{'end'}<$contigs->{$curprefix}{$prediction->{'sequence'}}{'length'}-$params->{'reblast_edge_distance'}){
				if($prediction->{'perc_identity'}<$params->{"reblast_min_perc_id"} or $prediction->{'end'}-$prediction->{'start'}<($params->{"reblast_min_perc_length"}/100.0)*$query_lengths{$prediction->{'query_seq'}}){
					# If hit touches neither contig end, min 90% ident, min 50% of query length
					$prediction->{'masked'} = 1;
				}
			}
		}
		print "Done.\n";

		$return_predictions{$curprefix} = \@combined_predictions;
	}
	print "\n";
	return \%return_predictions;
}

1;