#!/usr/bin/perl

# Re-BLAST module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_ReBlast;

use File::Path qw(make_path);
use List::Util qw(reduce);
use VICSIN;
use VH_helpers;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant REBLAST_DIR => "ReBlast_Runs";

sub run {
	# TODO redo this function to do one prefix at a time, not all together
	my $prefixes = shift;
	my $predictions = shift;
	my $masks = shift;
	my $contigs = shift;

	my %raw_predictions;
	my %return_predictions;

	VH_helpers::log("Starting blastn against predictions...");
	make_path(VICSIN::param("output_path")."/".REBLAST_DIR);
	foreach my $curprefix (@$prefixes){
		VH_helpers::log("\tRunning re-blast for $curprefix...",1);
		# Generate fasta file with predictions to query against
		VH_helpers::log("\t\tGenerating fasta file for re-blast... ",2);
		my $query_fasta_file_name = VICSIN::param("output_path")."/".REBLAST_DIR."/$curprefix-query.fna";
		open(my $queryfh, '>', $query_fasta_file_name) or die "Could not truncate query file.";
		my %query_lengths;
		foreach my $prefix (@{$prefixes}){
			if($curprefix ne $prefix){
				my $fasta_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$prefix.fna";
				open(my $fastafh, '<', $fasta_file_name) or die "Could not open fasta file.";
				my $seqname = "";
				my $sequence = "";
				
				while (my $row = <$fastafh>) {
					chomp $row;
					if(substr($row,0,1) eq '>'){
						if($seqname ne ""){
							# Grab predictions from sequence
							my $predCount = 1;
							for (my $bin = 0; $bin < 2; $bin++) {
								foreach my $prediction (@{$predictions->{$prefix}[$bin]}){
									if ($seqname eq $prediction->{'sequence'} and $prediction->{'end'}-$prediction->{'start'}>=VICSIN::param('reblast_min_contig_length')){
										print $queryfh '>'.$prediction->{'name'}."\n";
										print $queryfh substr($sequence,$prediction->{'start'}-1,$prediction->{'end'}-$prediction->{'start'})."\n";
										$query_lengths{$prediction->{'name'}} = $prediction->{'end'}-$prediction->{'start'}+1;
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
								print $queryfh '>'.$prediction->{'name'}."\n";
								print $queryfh substr($sequence,$prediction->{'start'}-1,$prediction->{'end'}-$prediction->{'start'})."\n";
								$query_lengths{$prediction->{'name'}} = $prediction->{'end'}-$prediction->{'start'}+1;
								$predCount++;
							}
						}
					}
				}
			}
		}

		# Re-blast predictions against all prefixes
		my $fasta_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$curprefix.fna";
		my $br_file_name = VICSIN::param("output_path")."/".REBLAST_DIR."/$curprefix.br";
		
		# Run blastn
		VH_helpers::log("\t\tRunning blastn... ",2);
		VH_helpers::run_cmd(VICSIN::param('blastn')." -query $query_fasta_file_name -subject $fasta_file_name -outfmt 6 -out $br_file_name");
		
		# Delete fasta file (it's huge)
		unlink($query_fasta_file_name);

		# Parse Results
		VH_helpers::log("\t\tParsing blast output... ",2);
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
			VH_helpers::log("\t\tApplying mask file... ",2);
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
		VH_helpers::log("\t\tMasking previous predictions... ",2);
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

		# get existing prediction counts
		my %predCount;
		for (my $bin = 0; $bin <= 3; $bin++) {
			foreach my $prediction (@{$predictions->{$curprefix}[$bin]}){
				if(not exists $predCount{$prediction->{'sequence'}}){
					$predCount{$prediction->{'sequence'}} = 0;
				}
				$predCount{$prediction->{'sequence'}}++;
			}
		}

		# For each R hit
		VH_helpers::log("\t\tChecking each contig for full coverage...",2);
		foreach my $contig (keys(%{$contigs->{$curprefix}})) {
			# If hits cover 90% of entire contig
			my @histogram = [];
			for(my $i=0; $i<$contigs->{$curprefix}{$contig}{'length'}; $i++){
				$histogram[$i] = 0;
			}
			for (my $i=0; $i<scalar(@masked_predictions); $i++){
				if($masked_predictions[$i]{'sequence'} eq $contig){
					my $start, $end;
					if($masked_predictions[$i]{'start'} < $masked_predictions[$i]{'end'}){
						$start = $masked_predictions[$i]{'start'};
						$end = $masked_predictions[$i]{'end'};
					} else {
						$start = $masked_predictions[$i]{'end'};
						$end = $masked_predictions[$i]{'start'};
					}
					for(my $j=$start; $j<=$end; $j++){
						$histogram[$j] = 1;
					}
				}
			}
			my $coverage = reduce {$a+$b} @histogram;
			VH_helpers::log("\t\t\tContig $contig coverage: ".$coverage." out of ".($contigs->{$curprefix}{$contig}{'length'}),2);
			if($coverage >= $contigs->{$curprefix}{$contig}{'length'} * 0.9 ) {
				# Add entire contig as prediction to bin 2?3?
				if(not exists $predCount{$contig}){
					$predCount{$contig} = 0;
				}
				push @{$predictions->{$curprefix}[3]}, {'name'=>$curprefix.'-'.$contig.'-'.$predCount{$contig},'sequence'=>$contig,'methods'=>'R','start'=>1,'end'=>$contigs->{$curprefix}{$contig}{'length'}};
				$predCount{$contig}++;
			}
		}

		VH_helpers::log("\t\tChecking each existing consensus prediction for extension by ReBLAST hit...",2);
		#For each prediction P:
		for (my $bin = 0; $bin <= 3; $bin++) {
			foreach my $prediction (@{$predictions->{$curprefix}[$bin]}){
				my $did_extend = 0;
				my $set_method = (index($prediction->{'methods'}, 'R') != -1);
				# Find any non end-adjacent overlapping hits
				do {
					my @overlaps = ();
					# For each hit R:
					for (my $i=0; $i<scalar(@masked_predictions); $i++){
						# %id and %length must be above thresholds
						if($masked_predictions[$i]{'sequence'} eq $prediction->{'sequence'}
							and $masked_predictions[$i]{'perc_identity'} >= VICSIN::param('reblast_min_perc_id')
							and abs($masked_predictions[$i]{'end'}-$masked_predictions[$i]{'start'})>=(VICSIN::param("reblast_min_perc_length")/100.0)*$query_lengths{$prediction->{'query_seq'}}){
							# if R overlaps P, add R to list of extenders
							my $start, $end;
							if($masked_predictions[$i]{'start'} < $masked_predictions[$i]{'end'}){
								$start = $masked_predictions[$i]{'start'};
								$end = $masked_predictions[$i]{'end'};
							} else {
								$start = $masked_predictions[$i]{'end'};
								$end = $masked_predictions[$i]{'start'};
							}
							# check for hit overlapping one or more ends, not totally contained
							if( $end>$prediction->{'start'}-VICSIN::param('reblast_distance') 
								and $start<$prediction->{'end'}+VICSIN::param('reblast_distance') 
								and not ($start>=$prediction->{'start'}-VICSIN::param('reblast_distance') and $end<=$prediction->{'end'}+VICSIN::param('reblast_distance')) ){
								push @overlaps, $masked_predictions[$i];
							}
						}
					}
					if(scalar(@overlaps)>0){
						$did_extend = 1;
						# Find largest extender
						my $maxlength = 0;
						my $largestextender;
						for (my $i=0; $i<scalar(@overlaps); $i++){
							if( abs($overlaps[$i]{'start'}-$overlaps[$i]{'end'})+1 >= $maxlength ){
								$maxlength = abs($overlaps[$i]{'start'}-$overlaps[$i]{'end'});
								$largestextender = $overlaps[$i];
							}
						}
						# Extend P with largest extender
						if($prediction->{'start'}>$largestextender->{'start'}){
							$prediction->{'start'} = $largestextender->{'start'};
						}
						if($prediction->{'end'}<$largestextender->{'end'}){
							$prediction->{'end'} = $largestextender->{'end'};
						}
						
						if(not $set_method){
							$prediction->{'methods'} = $prediction->{'methods'}.",R";
							$set_method = 1;
						}
					} else {
						$did_extend = 0;
					}
				} while ($did_extend); # While still something to extend with

				# One more pass to find edge-adjacent hits which overlap the now fully-extended prediction
				for(my $i=0; $i<scalar(@masked_predictions); $i++){
					# %id must be above threshold
					if($masked_predictions[$i]{'sequence'} eq $prediction->{'sequence'}
						and $masked_predictions[$i]{'perc_identity'} >= VICSIN::param('reblast_min_per_id')){
						my $start, $end;
						if($masked_predictions[$i]{'start'} < $masked_predictions[$i]{'end'}){
							$start = $masked_predictions[$i]{'start'};
							$end = $masked_predictions[$i]{'end'};
						} else {
							$start = $masked_predictions[$i]{'end'};
							$end = $masked_predictions[$i]{'start'};
						}
						# Check if hit spans from prediction to contig end
						if($start<$prediction->{'end'}+VICSIN::param('reblast_distance') and $end > $contigs->{$curprefix}{$masked_predictions[$i]{'sequence'}}{'length'} - VICSIN::param('reblast_edge_distance')){
							$prediction->{'end'} = $contigs->{$curprefix}{$masked_predictions[$i]{'sequence'}}{'length'};
							if(not $set_method){
								$prediction->{'methods'} = $prediction->{'methods'}.",R";
								$set_method = 1;
							}
						}
						# Check if hit spans from prediction to contig start
						if($end>$prediction->{'start'}-VICSIN::param('reblast_distance') and $start < 1 + VICSIN::param('reblast_edge_distance')){
							$prediction->{'start'} = 1;
							if(not $set_method){
								$prediction->{'methods'} = $prediction->{'methods'}.",R";
								$set_method = 1;
							}
						}
					}
				}
			}
		}
	}
	print "\n";
	return $predictions;
}

1;