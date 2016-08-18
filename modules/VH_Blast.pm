#!/usr/bin/perl

package VH_Blast;

use File::Path qw(make_path);
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant KNOWN_TYPES_DIR => "Known_Type_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	print VH_helpers->current_time()."Starting blastn against known viral types...\n";
	make_path($params->{"output_path"}."/".KNOWN_TYPES_DIR);
	foreach(@$prefixes){
		print VH_helpers->current_time()."\t$_... ";
		my $fasta_file_name = $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna";
		my $output_file_name = $params->{"output_path"}."/".KNOWN_TYPES_DIR."/$_.br";
		my $lock_file_name = $params->{"output_path"}."/".KNOWN_TYPES_DIR."/${_}_lock";
		if( -f $output_file_name and not -f $lock_file_name){
			print "$_ blastn already completed. Skipping. \n";
		} else {
			# Create a lockfile to signify that the blastn run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;
			# Run blastn
			`$params->{blastn} -query $params->{known_viral_types} -subject $fasta_file_name -outfmt 6 -out $output_file_name`;
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

	my $br_file_name = $params->{"output_path"}."/".KNOWN_TYPES_DIR."/${prefix}.br";

	my %predictions;
	open my $br_fh, '<', $br_file_name;
	while(my $br_line = <$br_fh>){
		chomp $br_line;
		if ($br_line =~ m/^.+?\s(.+?)\s.+?\s.+?\s.+?\s.+?\s.+?\s.+?\s(.+?)\s(.+?)\s/m){
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
		}
	}
	return \%predictions;
}

1;