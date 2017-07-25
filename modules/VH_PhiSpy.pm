#!/usr/bin/perl

# PhiSpy module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_PhiSpy;

use File::Path qw(make_path);
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant PHISPY_DIR => "PhiSpy_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	VH_helpers->log($params,"Starting PhiSpy runs...");
	foreach(@$prefixes) {
		my $seed_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_";
		my $phispy_dir_name = $params->{"output_path"}."/".PHISPY_DIR."/$_";

		my $tbl_file_name = $params->{"output_path"}."/".PHISPY_DIR."/$_/prophage.tbl";
		my $lock_file_name = $params->{"output_path"}."/".PHISPY_DIR."/$_/${_}_PhiSpy_lock";
		make_path($phispy_dir_name);

		# Check for lockfile to determine if phispy already complete
		if (-f $tbl_file_name and not -f $lock_file_name) {
			VH_helpers->log($params,"\t$_ PhiSpy already completed. Skipping.",1);
		} else {
			VH_helpers->log($params,"\tRunning PhiSpy for $_... ",1);
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			my $phispy_cmd = "python $params->{phispy} -i $seed_file_name -n $params->{phispy_threshold} -o $phispy_dir_name -w $params->{phispy_windowsize} -qt";
			VH_helpers->log($params, "\t\t$phispy_cmd",2);
			`$phispy_cmd`;

			VH_helpers->clean_folder($phispy_dir_name,[$tbl_file_name]);
			if ($? == 0){
			} else {
				VH_helpers->log($params,"PhiSpy returned an error on $_.");
			}
		}
	}
	print "\n";
}

sub get_predictions {
	shift;
	my $params = shift;
	my $prefix = shift(@_);

	my $tbl_file_name = $params->{"output_path"}."/".PHISPY_DIR."/$prefix/prophage.tbl";

	my %predictions;
	open my $tbl_fh, '<', $tbl_file_name;
	while(my $tbl_line = <$tbl_fh>){
		chomp $tbl_line;
		if ($tbl_line =~ m/^.+?\s(.+)_(.+?)_(.+?)$/m) {
			my $seq_name = $1;
			my $start = $2;
			my $end = $3;

			if(not exists $predictions{$seq_name}){
				$predictions{$seq_name} = [];
			}

			my $current_prediction = scalar(@{$predictions{$seq_name}});
			$predictions{$seq_name}[$current_prediction]{'start'} = $start;
			$predictions{$seq_name}[$current_prediction]{'end'} = $end;
			$predictions{$seq_name}[$current_prediction]{'phispy'} = [$current_prediction];
		}
	}
	return \%predictions;
}

1;