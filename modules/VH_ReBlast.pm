#!/usr/bin/perl

# Re-BLAST module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

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

	my %raw_predictions;
	my %return_predictions;

	VH_helpers->log($params,"Starting blastn against predictions...");
	make_path($params->{"output_path"}."/".REBLAST_DIR);
	foreach my $curprefix (@$prefixes){
		VH_helpers->log($params,"\tRunning re-blast for $curprefix...",1);
		# Generate fasta file with predictions to query against
		VH_helpers->log($params,"\t\tGenerating fasta file for re-blast... ",2);
		my $query_fasta_file_name = $params->{"output_path"}."/".REBLAST_DIR."/$curprefix-query.fna";
		open(my $queryfh, '>', $query_fasta_file_name) or die "Could not truncate query file.";
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
									if ($seqname eq $prediction->{'sequence'} and $prediction->{'end'}-$prediction->{'start'}>=$params->{'reblast_min_contig_length'}){
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

		# Re-blast predictions against all prefixes
		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$curprefix.fna";
		my $br_file_name = $params->{"output_path"}."/".REBLAST_DIR."/$curprefix.br";
		
		# Run blastn
		VH_helpers->log($params,"\t\tRunning blastn... ",2);
		my $blast_cmd = "$params->{blastn} -query $query_fasta_file_name -subject $fasta_file_name -outfmt 6 -out $br_file_name";
		VH_helpers->log($params, "\t\t$blast_cmd", 2);
		`$blast_cmd`;

		# Delete fasta file (it's huge)
		unlink($query_fasta_file_name);

		# Parse Results
		VH_helpers->log($params,"\t\tParsing blast output... ",2);
		my @blast_predictions;
		open my $br_fh, '<', $br_file_name;
		while(my $br_line = <$br_fh>){
			chomp $br_line;
			if ($br_line =~ m/^(.+?)\s(.+?)\s(.+?)\s.+?\s.+?\s.+?\s.+?\s.+?\s(.+?)\s(.+?)\s/m){
				my @br_array = split "\t", $br_line;
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

				push @blast_predictions, {'sequence'=>$seq_name,'methods'=>'R','start'=>$start,'end'=>$end,'query_seq'=>$query_seq,'perc_identity'=>$perc_identity, 'gap'=>$br_array[5], 'mismatch'=>$br_array[4], 'query_start'=>$br_array[6], 'query_stop'=>$br_array[7], 'bit'=>$br_array[11], 'evalue'=>$br_array[10]};
			}
		}

		# Mask results w/ masking file
		if(exists $masks->{$curprefix}){
			VH_helpers->log($params,"\t\tApplying mask file... ",2);
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
		}


		# Mask results w/ previous predictions
		VH_helpers->log($params,"\t\tMasking previous predictions... ",2);
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

		# Remove masked predictions from the list so the next step doesn't take quite so long
		my @masked_predictions;
		foreach my $prediction (@blast_predictions){
			if(not exists $prediction->{'masked'}){
				push @masked_predictions, $prediction;
			}
		}

		# Combine overlapping results
		VH_helpers->log($params,"\t\tCombining results... ",2);
		my $found_overlaps;
		my @combined_predictions;
		do {
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

		# Weed out results based on length
		VH_helpers->log($params,"\t\tExcluding based on length/% identity... ",2);
		my @final_predictions;
		my $cur = 0;
		foreach my $prediction (@combined_predictions){
			# TODO implement new rules based on end-ness
			if($prediction->{'start'} <= $params->{'reblast_edge_distance'} xor $prediction->{'end'}>=$contigs->{$curprefix}{$prediction->{'sequence'}}{'length'}-$params->{'reblast_edge_distance'}){
				if($prediction->{'perc_identity'} >= $params->{"reblast_min_perc_id"}){
					# If hit touches either contig end (not both), min 90% ident
					$prediction->{'reblast'} = [$cur];
					@final_predictions[$cur] = $prediction;
					$cur++;
				}
			} elsif ($prediction->{'start'}>$params->{'reblast_edge_distance'} and $prediction->{'end'}<$contigs->{$curprefix}{$prediction->{'sequence'}}{'length'}-$params->{'reblast_edge_distance'}){
				if($prediction->{'perc_identity'}>=$params->{"reblast_min_perc_id"} and $prediction->{'end'}-$prediction->{'start'}>=($params->{"reblast_min_perc_length"}/100.0)*$query_lengths{$prediction->{'query_seq'}}){
					# If hit touches neither contig end, min 90% ident, min 50% of query length
					$prediction->{'reblast'} = [$cur];
					@final_predictions[$cur] = $prediction;
					$cur++;
				}
			}
		}

		$return_predictions{$curprefix} = \@final_predictions;
	}
	print "\n";
	return \%return_predictions;
}

1;