#!/usr/bin/perl

package VH_PhiSpy;

use File::Path qw(make_path);
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant PHISPY_DIR => "PhiSpy_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	print VH_helpers->current_time()."Starting PhiSpy runs...\n";
	foreach(@$prefixes) {
		print VH_helpers->current_time()."\t$_... ";

		my $seed_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_";
		my $phispy_dir_name = $params->{"output_path"}."/".PHISPY_DIR."/$_";

		my $tbl_file_name = $params->{"output_path"}."/".PHISPY_DIR."/$_/prophage.tbl";
		my $lock_file_name = $params->{"output_path"}."/".PHISPY_DIR."/$_/${_}_PhiSpy_lock";
		make_path($phispy_dir_name);

		# Check for lockfile to determine if phispy already complete
		if (-f $tbl_file_name and not -f $lock_file_name) {
			print "$_ PhiSpy already completed. Skipping.\n";
		} else {
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			`python $params->{phispy} -i $seed_file_name -n $params->{phispy_threshold} -o $phispy_dir_name -w $params->{phispy_windowsize} -qt`;

			VH_helpers::clean_folder($phispy_dir_name,[$tbl_file_name]);
			if ($? == 0){
				print "Done.\n";
			} else {
				print "PhiSpy returned an error.\n";
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
		}
	}
	return \%predictions;
}

1;