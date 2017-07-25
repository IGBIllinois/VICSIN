#!/usr/bin/perl

# CRISPR Match BLAST module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_CRISPR;

use File::Path qw(make_path);
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant CRISPR_DIR => "CRISPR_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	VH_helpers->log($params,"Starting CRISPR runs...");

	foreach(@$prefixes) {
		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna";
		my $crispr_file_name = $params->{"output_path"}."/".CRISPR_DIR."/${_}_CRISPR.br";
		my $lock_file_name = $params->{'output_path'}.'/'.CRISPR_DIR."/${_}_CRISPR_lock";
		make_path($params->{"output_path"}."/".CRISPR_DIR);

		# If the crispr file exists but not the lock file, the CRISPR run was already complete
		if (-f $crispr_file_name and not -f $lock_file_name){
			VH_helpers->log($params,"$_ CRISPR already completed. Skipping.",1);
		} else {
			VH_helpers->log($params,"\tRunning CRISPR blast for $_...",1);
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			my $crispr_cmd = "$params->{blastn} -task \"blastn-short\" -subject $fasta_file_name -query $params->{spacer_fasta_file} -outfmt 6 -out $crispr_file_name";
			VH_helpers->log($params, "\t\t$crispr_cmd", 2);
			`$crispr_cmd`;
			
			unlink $lock_file_name;
		}
	}
	print "\n";
}

sub get_predictions {
	shift;
	my $params = shift;
	my $prefix = shift(@_);

	my $br_file_name = $params->{"output_path"}."/".CRISPR_DIR."/${prefix}_CRISPR.br";

	my %predictions;
	open my $br_fh, '<', $br_file_name;
	while(my $br_line = <$br_fh>){
		chomp $br_line;
		if ($br_line =~ m/^.+?\s(.+?)\s.+?\s.+?\s.+?\s.+?\s.+?\s.+?\s(.+?)\s(.+?)\s/){
			my @br_array = split "\t", $br_line;
			my $seq_name = $1;
			my $start = $2;
			my $end = $3;

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
			$predictions{$seq_name}[$current_prediction]{'crispr'} = [$current_prediction];
		}
	}
	return \%predictions;
}

1;