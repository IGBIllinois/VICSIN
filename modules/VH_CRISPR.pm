#!/usr/bin/perl

package VH_CRISPR;

use File::Path qw(make_path);
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant CRISPR_DIR => "CRISPR_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	print VH_helpers->current_time()."Starting CRISPR runs...\n";

	foreach(@$prefixes) {
		print VH_helpers->current_time()."\t$_...";

		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna";
		my $crispr_file_name = $params->{"output_path"}."/".CRISPR_DIR."/${_}_CRISPR.br";
		my $lock_file_name = $params->{'output_path'}.'/'.CRISPR_DIR."/${_}_CRISPR_lock";
		make_path($params->{"output_path"}."/".CRISPR_DIR);

		# If the crispr file exists but not the lock file, the CRISPR run was already complete
		if (-f $crispr_file_name and not -f $lock_file_name){
			print "$_ CRISPR already completed. Skipping.\n";
		} else {
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			`$params->{blastn} -task "blastn-short" -subject $fasta_file_name -query $params->{spacer_fasta_file} -outfmt 6 -out $crispr_file_name`;

			unlink $lock_file_name;
			print "Done.\n";
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
			my $seq_name = $1;
			my $start = $2;
			my $end = $3;

			if(not exists $predictions{$seq_name}){
				$predictions{$seq_name}=[];
			}

			my $current_prediction = scalar(@{$predictions{$seq_name}});
			$predictions{$seq_name}[$current_prediction]{'start'} = $start;
			$predictions{$seq_name}[$current_prediction]{'end'} = $end;
		}
	}
	return \%predictions;
}

1;