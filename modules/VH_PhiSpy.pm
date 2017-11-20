#!/usr/bin/perl

# PhiSpy module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_PhiSpy;

use File::Path qw(make_path);
use VICSIN;
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant PHISPY_DIR => "PhiSpy_Runs";

sub run {
	my $prefixes = shift;

	VH_helpers::log("Starting PhiSpy runs...");
	foreach(@$prefixes) {
		my $seed_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/_SEED_$_";
		my $phispy_dir_name = VICSIN::param("output_path")."/".PHISPY_DIR."/$_";

		my $tbl_file_name = VICSIN::param("output_path")."/".PHISPY_DIR."/$_/prophage.tbl";
		my $lock_file_name = VICSIN::param("output_path")."/".PHISPY_DIR."/$_/${_}_PhiSpy_lock";
		my $log_file_name = VICSIN::param("output_path")."/".PHISPY_DIR."/$_/log.txt";
		make_path($phispy_dir_name);

		# Check for lockfile to determine if phispy already complete
		if (-f $tbl_file_name and not -f $lock_file_name) {
			VH_helpers::log("\t$_ PhiSpy already completed. Skipping.",1);
		} else {
			VH_helpers::log("\tRunning PhiSpy for $_... ",1);
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			VH_helpers::run_cmd("python ".VICSIN::param('phispy')." -i $seed_file_name -n ".VICSIN::param('phispy_threshold')." -o $phispy_dir_name -w ".VICSIN::param('phispy_windowsize')." 2>&1 >$log_file_name");

			VH_helpers::clean_folder($phispy_dir_name,[$tbl_file_name,$log_file_name]);
			if ($? == 0){
			} else {
				VH_helpers::log("PhiSpy returned an error on $_.");
			}
		}
	}
	print "\n";
}

sub get_predictions {
	my $prefix = shift(@_);

	my $tbl_file_name = VICSIN::param("output_path")."/".PHISPY_DIR."/$prefix/prophage.tbl";

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