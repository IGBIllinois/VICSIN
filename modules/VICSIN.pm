#! /usr/bin/env perl

# Contains main VICSIN functionality
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VICSIN;

use strict;
use Getopt::Long qw(GetOptions);
use YAML qw(LoadFile);

use VH_helpers;

use Data::Dumper;

my %params = (
	"input_path"=>".",
	"output_path"=>".",
	"genbank_to_seed"=>"genbank_to_seed.py",
	"genbank_to_fasta"=>"genbank_to_fasta.py",
	"core_genome"=>"core_genome.py",
	"gff_to_seed"=>"gff_to_seed.pl",
	"prodigal"=>"prodigal",
	"virsorter"=>"wrapper_phage_contigs_sorter_iPlant.pl",
	"virsorter_data_dir"=>"/data", # This is the virsorter's builtin default
	"phispy"=>"phispy.py",
	"blastn"=>"blastn",
	"makeblastdb"=>"makeblastdb",
	"spine"=>"spine.pl",
	"agent"=>"AGEnt.pl",
	"cut"=>"cut",
	"mcxload"=>"mcxload",
	"mcl"=>"mcl",
	"mcxdump"=>"mcxdump",
	"blast_to_mcl"=>"Blast_to_MCL.1.py",
	"mcldump2clusters"=>"MCLdump2clusters.pl",
	"pamprotopatterngrab"=>"PAMProtoPatternGrab_full.py",
	"phispy_windowsize"=>40,
	"phispy_threshold"=>20,
	"spine_percent_input"=>100,
	"spine_max_distance"=>10,
	"spine_agent_min_perc_id"=>85,
	"spine_agent_min_size_core"=>10,
	"spine_core_file"=>"",
	"spacer_fasta_file"=>"",
	"crispr_match_threshold"=>0.9,
	"known_viral_types"=>"",
	"virsorter_database"=>2,
	"masking_file"=>"",
	"mcl_inflation"=>2.0,
	"overhang_threshold"=>100,
	"merge_threshold"=>1500,
	"clustering_parameter"=>"percent_length_aligned",
	"reblast_min_perc_id"=>90,
	"reblast_min_perc_length"=>50,
	"reblast_edge_distance"=>5,
	"reblast_distance"=>5,
	"reblast_min_contig_length"=>10000,
	"cluster_core_congruence"=>0.66,
	"percent_id_min_core"=>85,
	"cluster_min_perc_length"=>0.5,
	"cluster_min_length"=>100,
	"cluster_min_bit_score"=>100.0,
	"cluster_core_max_distance"=>10,
	"cluster_size_threshold"=>2000,
	"use_database"=>'false',
	"database_host"=>'',
	"database_name"=>'',
	"database_port"=>3306,
	"database_user"=>'',
	"database_pass"=>'',
	"verbosity"=>0,
	"stop"=>'',
	"skip"=>'',
	"num_threads"=>1
);

our @methods = (
	{'key'=>'agent','abbr'=>'A','bin'=>3,'extend'=>1},
	{'key'=>'virsorter','abbr'=>'V','bin'=>2,'extend'=>1},
	{'key'=>'phispy','abbr'=>'P','bin'=>3,'extend'=>0},
	{'key'=>'blast','abbr'=>'B','bin'=>2,'extend'=>1},
	{'key'=>'crispr','abbr'=>'C','bin'=>4,'extend'=>0},
	{'key'=>'reblast','abbr'=>'R','bin'=>2,'extend'=>0}
);

sub parseParameters {
	my $config_file = shift;

	# First, read in the config file, if provided
	if($config_file ne ""){
		my $config_hash = LoadFile($config_file);
		foreach my $k (keys %params) {
			$params{$k} = $config_hash->{$k} if exists $config_hash->{$k};
		}
	}

	# Override with command line arguments, if provided
	GetOptions(	"input_path=s"=>\$params{"input_path"},
				"output_path=s"=>\$params{"output_path"},
				"genbank_to_seed=s"=>\$params{"genbank_to_seed"},
				"genbank_to_fasta=s"=>\$params{"genbank_to_fasta"},
				"core_genome=s"=>\$params{"core_genome"},
				"gff_to_seed=s"=>\$params{"gff_to_seed"},
				"prodigal=s"=>\$params{"prodigal"},
				"virsorter=s"=>\$params{"virsorter"},
				"virsorter_data_dir=s"=>\$params{"virsorter_data_dir"},
				"phispy=s"=>\$params{"phispy"},
				"blastn=s"=>\$params{"blastn"},
				"makeblastdb=s"=>\$params{"makeblastdb"},
				"spine=s"=>\$params{"spine"},
				"agent=s"=>\$params{"agent"},
				"cut=s"=>\$params{"cut"},
				"mcxload=s"=>\$params{"mcxload"},
				"mcl=s"=>\$params{"mcl"},
				"mcxdump=s"=>\$params{"mcxdump"},
				"blast_to_mcl=s"=>\$params{"blast_to_mcl"},
				"mcldump2clusters=s"=>\$params{"mcldump2clusters"},
				"pamprotopatterngrab=s"=>\$params{"pamprotopatterngrab"},
				"phispy_windowsize=i"=>\$params{"phispy_windowsize"},
				"phispy_threshold=i"=>\$params{"phispy_threshold"},
				"spine_percent_input=i"=>\$params{"spine_percent_input"},
				"spine_max_distance=i"=>\$params{"spine_max_distance"},
				"spine_agent_min_perc_id=i"=>\$params{"spine_agent_min_perc_id"},
				"spine_agent_min_size_core=i"=>\$params{"spine_agent_min_size_core"},
				"spine_core_file=s"=>\$params{"spine_core_file"},
				"spacer_fasta_file=s"=>\$params{"spacer_fasta_file"},
				"crispr_match_threshold=f"=>\$params{"crispr_match_threshold"},
				"known_viral_types=s"=>\$params{"known_viral_types"},
				"virsorter_database=i"=>\$params{"virsorter_database"},
				"masking_file=s"=>\$params{"masking_file"},
				"mcl_inflation=f"=>\$params{"mcl_inflation"},
				"overhang_threshold=i"=>\$params{"overhang_threshold"},
				"merge_threshold=i"=>\$params{"merge_threshold"},
				"clustering_parameter=s"=>\$params{"clustering_parameter"},
				"reblast_min_perc_id=i"=>\$params{"reblast_min_perc_id"},
				"reblast_min_perc_length=i"=>\$params{"reblast_min_perc_length"},
				"reblast_edge_distance=i"=>\$params{"reblast_edge_distance"},
				"reblast_distance=i"=>\$params{"reblast_distance"},
				"reblast_min_contig_length=i"=>\$params{"reblast_min_contig_length"},
				"cluster_core_congruence=f"=>\$params{"cluster_core_congruence"},
				"percent_id_min_core=i"=>\$params{"percent_id_min_core"},
				"cluster_min_perc_length=f"=>\$params{"cluster_min_perc_length"},
				"cluster_min_length=i"=>\$params{"cluster_min_length"},
				"cluster_min_bit_score=f"=>\$params{"cluster_min_bit_score"},
				"cluster_core_max_distance=i"=>\$params{"cluster_core_max_distance"},
				"cluster_size_threshold=i"=>\$params{"cluster_size_threshold"},
				"use_database=s"=>\$params{"use_database"},
				"database_host=s"=>\$params{"database_host"},
				"database_name=s"=>\$params{"database_name"},
				"database_port=i"=>\$params{"database_port"},
				"database_user=s"=>\$params{"database_user"},
				"database_pass=s"=>\$params{"database_pass"},
				"verbosity=i"=>\$params{"verbosity"},
				"stop=s"=>\$params{"stop"},
				"skip=s"=>\$params{"skip"},
				"num_threads=i"=>\$params{"num_threads"}
				);
}

sub logParameters {
	foreach my $k (sort keys %params) {
		VH_helpers::log("\t$k: $params{$k}",1);
	}
}

sub param {
	my $paramname = shift;
	if(exists $params{$paramname}) {
		return $params{$paramname};
	} else {
		return;
	}
}
sub setParam {
	my $paramname = shift;
	my $paramvalue = shift;
	if(exists $params{$paramname}){
		$params{$paramname} = $paramvalue;
	}
}

##### Consensus processing functions #####
# Finds overlaps between prediction a and list of predictions b; trims overhanging edges in b
# TODO this function has way too many side effects
sub overlap_exists {
	my ($a,$b) = @_;
	
	my $found_overlap = 0;
	my @predictions_to_add;
	if( not exists $a->{'used'} and not exists $a->{'masked'} ){
		foreach my $prediction (@$b){
			my $this_prediction_hit = 0;
			if(not (exists $prediction->{'used'} and $prediction->{'used'} == 1) and not exists $prediction->{'masked'} ){

				if($a->{'start'}<=$prediction->{'start'} and $a->{'end'}>=$prediction->{'end'}){
					# Prediction entirely contained within a
					$found_overlap = 1;
					$this_prediction_hit = 1;
					$prediction->{'used'} = 1;
				} elsif ($a->{'start'}<=$prediction->{'start'} and $a->{'end'}>=$prediction->{'start'}){
					#Prediction overhangs right side of a
					$found_overlap = 1;
					$this_prediction_hit = 1;
					if($prediction->{'end'}-$a->{'end'}>=VICSIN::param('overhang_threshold')){
						push @predictions_to_add, {'start'=>$a->{'end'}+1,'end'=>$prediction->{'end'}};
						$prediction->{'end'} = $a->{'end'};
					}
					$prediction->{'used'} = 1;
				} elsif ($a->{'start'}<=$prediction->{'end'} and $a->{'end'}>=$prediction->{'end'}){
					# Prediction overhangs left side of a
					$found_overlap = 1;
					$this_prediction_hit = 1;
					if($a->{'start'}-$prediction->{'start'}>=VICSIN::param('overhang_threshold')){
						push @predictions_to_add, {'start'=>$prediction->{'start'},'end'=>$a->{'start'}-1};
						$prediction->{'start'} = $a->{'start'};
					}
					$prediction->{'used'} = 1;
				} elsif ($a->{'start'}>$prediction->{'start'} and $a->{'end'}<$prediction->{'end'}){
					# Prediction overhangs a on both sides
					$found_overlap = 1;
					$this_prediction_hit = 1;
					if($a->{'start'}-$prediction->{'start'}>=VICSIN::param('overhang_threshold')){
						push @predictions_to_add, {'start'=>$prediction->{'start'},'end'=>$a->{'start'}-1};
						$prediction->{'start'} = $a->{'start'};
					}
					if($prediction->{'end'}-$a->{'end'}>=VICSIN::param('overhang_threshold')){
						push @predictions_to_add, {'start'=>$a->{'end'}+1,'end'=>$prediction->{'end'}};
						$prediction->{'end'} = $a->{'end'};
					}
					$prediction->{'used'} = 1;
				}
			}
			if($this_prediction_hit == 1){
				# Keep track of which hits are incorporated here
				foreach my $method (@VICSIN::methods){
					if(exists $prediction->{$method->{'key'}}){
						if(not exists $a->{$method->{'key'}}){
							$a->{$method->{'key'}} = [];
						}
						push @{$a->{$method->{'key'}}}, @{$prediction->{$method->{'key'}}};
					}
				}
			}
		}
		push @{$b},@predictions_to_add;
	}
	return $found_overlap;
}

# Merges the predictions in the given array as much as possible
# TODO ignore masked-out predictions
sub merge_predictions {
	my $mergeable_predictions = shift;

	# Merge until there ain't no more to merge
	my $found_merges;
	my %merged_predictions;
	foreach my $sequence (keys %{$mergeable_predictions}){
		do {
			$merged_predictions{$sequence} = [];
			$found_merges = 0;
			for (my $i = 0; $i < scalar(@{$mergeable_predictions->{$sequence}}); $i++){
				if ( not (exists $mergeable_predictions->{$sequence}[$i]{'masked'} and $mergeable_predictions->{$sequence}[$i]{'masked'} == 1) ){ # Ignore masked out predictions
					my $found_match = 0; 
					for (my $j=0; $j < scalar(@{$merged_predictions{$sequence}}); $j++){
						# If predictions overlap or are within merge_threshold of each other
						if($found_match == 0 
							and not (exists $merged_predictions{$sequence}[$j]{'masked'} and $merged_predictions{$sequence}[$j]{'masked'} == 1) # Ignore masked out predictions
							and $mergeable_predictions->{$sequence}[$i]{'start'}<=$merged_predictions{$sequence}[$j]{'end'}+VICSIN::param('merge_threshold')
							and $merged_predictions{$sequence}[$j]{'start'}<=$mergeable_predictions->{$sequence}[$i]{'end'}+VICSIN::param('merge_threshold')){

							$found_match = 1;
							$found_merges = 1;
							# Expand bounds of j to contain i
							if($mergeable_predictions->{$sequence}[$i]{'start'}<$merged_predictions{$sequence}[$j]{'start'}){
								$merged_predictions{$sequence}[$j]{'start'} = $mergeable_predictions->{$sequence}[$i]{'start'};
							}
							if($mergeable_predictions->{$sequence}[$i]{'end'}>$merged_predictions{$sequence}[$j]{'end'}){
								$merged_predictions{$sequence}[$j]{'end'} = $mergeable_predictions->{$sequence}[$i]{'end'};
							}
							# merge methods string
							foreach my $method (@{$mergeable_predictions->{$sequence}[$i]{'methods'}}){
								if(not grep(/^$method$/,@{$merged_predictions{$sequence}[$j]{'methods'}})){
									push @{$merged_predictions{$sequence}[$j]{'methods'}}, $method;
								}
							}

							# Keep track of which hits are incorporated here
							foreach my $method (@VICSIN::methods){
								if(exists $mergeable_predictions->{$sequence}[$i]{$method->{'key'}}){
									if(not exists $merged_predictions{$sequence}[$j]{$method->{'key'}}){
										$merged_predictions{$sequence}[$j]{$method->{'key'}} = [];
									}
									push @{$merged_predictions{$sequence}[$j]{$method->{'key'}}}, @{$mergeable_predictions->{$sequence}[$i]{$method->{'key'}}};
								}
							}
						}
					}
					if($found_match == 0 ){
						# If this prediction didn't overlap with any others, add it to the array
						push @{$merged_predictions{$sequence}}, $mergeable_predictions->{$sequence}[$i];
					}
				}
			}
			$mergeable_predictions->{$sequence} = [];
			push @{$mergeable_predictions->{$sequence}},@{$merged_predictions{$sequence}};
		} while ($found_merges == 1);
	}

	return %merged_predictions;
}

sub bin_predictions {
	my $merged_predictions = shift;
	my $predictions = shift;
	my $prefix = shift;

	my @binned_predictions = ([],[],[],[],[]);
	# First pass: compare non-extending predictions with extending predictions
	my %predCount;
	foreach my $sequence (keys %{$merged_predictions}){
		$predCount{$sequence} = 0;
		for (my $j = 0; $j < scalar(@{$merged_predictions->{$sequence}}); $j++){
			for(my $i=0; $i<scalar(@VICSIN::methods); $i++){
				if($VICSIN::methods[$i]{'extend'} == 0){
					if(exists $predictions->{$prefix}{$VICSIN::methods[$i]{'key'}}{$sequence}){
						if(VICSIN::overlap_exists($merged_predictions->{$sequence}[$j],$predictions->{$prefix}{$VICSIN::methods[$i]{'key'}}{$sequence}) == 1){
							push @{$merged_predictions->{$sequence}[$j]{'methods'}},$VICSIN::methods[$i]{'abbr'};
						}
					}
				}
			}
			my $bin;
			if(scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) > 2){
				$bin = 0;
			} elsif(scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) == 2) {
				$bin = 1;
			} elsif (scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) == 1){
				for(my $k=0; $k<scalar(@VICSIN::methods); $k++){
					if($VICSIN::methods[$k]{'abbr'} eq $merged_predictions->{$sequence}[$j]{'methods'}[0]){
						$bin = $VICSIN::methods[$k]{'bin'};
					}
				}
			}

			my %final_prediction = ('name'=>$prefix.'-'.$sequence.'-'.$predCount{$sequence},'sequence'=>$sequence,'methods'=>join(',',@{$merged_predictions->{$sequence}[$j]{'methods'}}),'start'=>$merged_predictions->{$sequence}[$j]{'start'},'end'=>$merged_predictions->{$sequence}[$j]{'end'});
			$predCount{$sequence}++;
			foreach my $method (@VICSIN::methods){
				if(exists $merged_predictions->{$sequence}[$j]{$method->{'key'}}){
					$final_prediction{$method->{'key'}} = $merged_predictions->{$sequence}[$j]{$method->{'key'}};
				}
			}
			push @{$binned_predictions[$bin]}, \%final_prediction;
		}
	}

	# Second pass: compare non-extending predictions amongst themselves
	for(my $i=0;$i<scalar @VICSIN::methods; $i++) {
		if($VICSIN::methods[$i]{'extend'}==0){ # Only compare methods we haven't already merged
			foreach my $sequence (keys %{ $predictions->{$prefix}{$VICSIN::methods[$i]{'key'}}}){
				foreach my $prediction ( @{ $predictions->{$prefix}{$VICSIN::methods[$i]{'key'}}{$sequence} }){
					if(not exists $prediction->{'used'} and $prediction->{'start'}>0 and $prediction->{'end'}>0){ 
						my $used_methods = $VICSIN::methods[$i]{'abbr'};
						my $matches = 1;
						for(my $j=$i+1; $j<scalar @VICSIN::methods; $j++){
							if($VICSIN::methods[$j]{'extend'}==0){
								if(exists $predictions->{$prefix}{$VICSIN::methods[$j]{'key'}}{$sequence}){
									if(VICSIN::overlap_exists($prediction,$predictions->{$prefix}{$VICSIN::methods[$j]{'key'}}{$sequence}) == 1){
										$matches++;
										$used_methods .= ",".$VICSIN::methods[$j]{'abbr'};
									}
								}
							}
						}
						my $bin;
						if($matches > 2){
							$bin = 0;
						} elsif($matches == 2){
							$bin = 1;
						} elsif($matches == 1){
							$bin = $VICSIN::methods[$i]{'bin'};
						}

						my %final_prediction = ('name'=>$prefix.'-'.$sequence.'-'.$predCount{$sequence},'sequence'=>$sequence,'methods'=>$used_methods,'start'=>$prediction->{'start'},'end'=>$prediction->{'end'});
						$predCount{$sequence}++;
						foreach my $method (@VICSIN::methods){
							if(exists $prediction->{$method->{'key'}}){
								$final_prediction{$method->{'key'}} = $prediction->{$method->{'key'}};
							}
						}
						push @{$binned_predictions[$bin]}, \%final_prediction;
					}
				}
			}
		}
	}

	return \@binned_predictions;
}

sub apply_mask {
	my $predictions = shift;
	my $prefix = shift;
	my $masks = shift;

	if(exists $masks->{$prefix}){
		VH_helpers::log("\t\tApplying mask... ",1);
		foreach my $sequence ( keys %{$masks->{$prefix}} ){
			foreach my $mask ( @{$masks->{$prefix}{$sequence}} ){
				foreach my $method ( keys %{$predictions->{$prefix}} ){
					if(exists $predictions->{$prefix}{$method}{$sequence}){
						my @predictions_to_add;
						foreach my $prediction (@{$predictions->{$prefix}{$method}{$sequence}}){
							if($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'end'}){
								# Prediction entirely contained within mask; set to ignore prediction
								$prediction->{'start'} = -1; # Ignoring is _much_ easier than removing from the array
								$prediction->{'end'} = -1;
								$prediction->{'masked'} = 1;
							} elsif ($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'start'}) {
								# Prediction's left side overlaps mask; trim left side
								$prediction->{'start'} = $mask->{'end'}+1;
							} elsif ($mask->{'start'}<=$prediction->{'end'} and $mask->{'end'}>=$prediction->{'end'}) {
								# Prediction's right side overlaps mask; trim right side
								$prediction->{'end'} = $mask->{'start'}-1;
							} elsif ($mask->{'start'}>$prediction->{'start'} and $mask->{'end'}<$prediction->{'end'} and not exists $prediction->{'used'}) {
								# Prediction contains a masked area within it; set to ignore prediction
								$prediction->{'start'} = -1; # Ignoring is _much_ easier than removing from the array
								$prediction->{'end'} = -1;
								$prediction->{'masked'} = 1;
								# Prediction contains a masked area within it; split into two predictions. Don't split ignored predictions.
								# push @predictions_to_add, {'start'=>$mask->{'end'}+1,'end'=>$prediction->{'end'}};
								# $prediction->{'end'} = $mask->{'start'}-1;
							}
						}
						push @{$predictions->{$prefix}{$method}{$sequence}}, @predictions_to_add;
					}
				}
			}
		}
	}
}

1;