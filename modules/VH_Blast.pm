#!/usr/bin/perl

# Homology BLAST module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_Blast;

use File::Path qw(make_path);
use VICSIN;
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant KNOWN_TYPES_DIR => "Known_Type_Runs";

sub run {
	my $prefixes = shift;

	VH_helpers::log("Starting blastn against known viral types...");
	make_path(VICSIN::param("output_path")."/".KNOWN_TYPES_DIR);
	foreach(@$prefixes){
		my $fasta_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$_.fna";
		my $output_file_name = VICSIN::param("output_path")."/".KNOWN_TYPES_DIR."/$_.br";
		my $log_file_name = VICSIN::param("output_path")."/".KNOWN_TYPES_DIR."/$_-log.txt";
		my $lock_file_name = VICSIN::param("output_path")."/".KNOWN_TYPES_DIR."/${_}_lock";
		if( -f $output_file_name and not -f $lock_file_name){
			VH_helpers::log("\t$_ homology blast already completed. Skipping.");
		} else {
			VH_helpers::log("\tRunning homology blast for $_... ",1);
			# Create a lockfile to signify that the blastn run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;
			
			# Run blastn
			VH_helpers::run_cmd(VICSIN::param('blastn')." -query ".VICSIN::param('known_viral_types')." -subject $fasta_file_name -outfmt 6 -out $output_file_name 2>&1 >$log_file_name");
			
			unlink $lock_file_name;
		}
	}
	print "\n";
}

sub get_predictions {
	my $prefix = shift;

	my $br_file_name = VICSIN::param("output_path")."/".KNOWN_TYPES_DIR."/${prefix}.br";

	my %predictions;
	open my $br_fh, '<', $br_file_name;
	while(my $br_line = <$br_fh>){
		chomp $br_line;
		if ($br_line =~ m/^.+?\s(.+?)\s.+?\s.+?\s.+?\s.+?\s.+?\s.+?\s(.+?)\s(.+?)\s/m){
			my @br_array = split "\t", $br_line;
			my $seq_name = $1;
			my $start;
			my $end;
			if($2<$3){
				$start = $2;
				$end = $3;	
			} else {
				$start = $3;
				$end = $2;
			}

			if(not exists $predictions{$seq_name}){
				$predictions{$seq_name}=[];
			}

			my $current_prediction = scalar(@{$predictions{$seq_name}});
			$predictions{$seq_name}[$current_prediction]{'start'} = $start;
			$predictions{$seq_name}[$current_prediction]{'end'} = $end;
			$predictions{$seq_name}[$current_prediction]{'query'} = $br_array[0];
			$predictions{$seq_name}[$current_prediction]{'perc_id'} = $br_array[2];
			$predictions{$seq_name}[$current_prediction]{'gap'} = $br_array[5];
			$predictions{$seq_name}[$current_prediction]{'mismatch'} = $br_array[4];
			$predictions{$seq_name}[$current_prediction]{'query_start'} = $br_array[6];
			$predictions{$seq_name}[$current_prediction]{'query_stop'} = $br_array[7];
			$predictions{$seq_name}[$current_prediction]{'bit'} = $br_array[11];
			$predictions{$seq_name}[$current_prediction]{'evalue'} = $br_array[10];
			$predictions{$seq_name}[$current_prediction]{'blast'} = [$current_prediction];
		}
	}
	return \%predictions;
}

1;