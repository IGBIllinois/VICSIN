#!/usr/bin/perl

package VH_VirSorter;

use File::Path qw(make_path);
use File::Copy;
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant VIRSORTER_DIR => "Virsorter_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	print VH_helpers->current_time()."Starting VirSorter runs...\n";
	foreach(@$prefixes) {
		print VH_helpers->current_time()."\t$_... ";

		my $fasta_file_name =  File::Spec->rel2abs( $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna" );
		my $wdir = $params->{"output_path"}."/".VIRSORTER_DIR."/$_";
		my $data_dir = File::Spec->rel2abs( $params->{"virsorter_data_dir"} );

		my $lock_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/${_}/${_}_VirSorter_lock";
		my $csv_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/${_}/${_}_global-phage-signal.csv";
		my $mga_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/${_}/Fasta_files/${_}_mga_final.predict";
		my $mga_dest_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/${_}/${_}_mga_final.predict";

		# If the virsorter files exist but not the lock file, virsorter previously completed
		if ( -f $csv_file_name and -f $mga_dest_file_name and not -f $lock_file_name ){ 
			print "$_ VirSorter already completed. Skipping.\n";
		} else {
			make_path($wdir);
			# Create a lockfile to signify that the VirSorter run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			`cd $wdir; $params->{virsorter} -d $_ --fna $fasta_file_name --db $params->{virsorter_database} --data-dir $data_dir 2>&1; cd -`;
			
			# Move mga file to its final destination
			move($mga_file_name,$mga_dest_file_name);
			VH_helpers::clean_folder($wdir,[$csv_file_name,$mga_dest_file_name,$wdir."/log_out",$wdir."/log_err"]);

			print "Done.\n";
		}
	}
	print "\n";
}

sub get_predictions {
	shift;
	my $params = shift;
	my $prefix = shift;

	my $csv_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/$prefix/${prefix}_global-phage-signal.csv";
	my $mga_file_name = $params->{"output_path"}."/".VIRSORTER_DIR."/$prefix/${prefix}_mga_final.predict";

	# Parse csv
	my %genes;
	my %genefrags;
	open my $csv_fh, '<', $csv_file_name;
	while (my $csv_line = <$csv_fh>){
		if($csv_line =~ m/^([^#>]+?),.*?,(.+?),(.+?),/){
			my $sequence = $1;
			my $fragment = $2;
			my $fragsize = $3;
			if(not exists $genefrags{$sequence}){
				$genefrags{$sequence} = 0;
			} else {
				$genefrags{$sequence} = $genefrags{$sequence}+1;
			}
			if($fragment =~ m/.+?-gene_(.+?)-gene_(.+?)$/m){
				$genes{$sequence}[$genefrags{$sequence}]{'start'} = $1;
				$genes{$sequence}[$genefrags{$sequence}]{'end'} = $2;
			} else {
				$genes{$sequence}[$genefrags{$sequence}]{'start'} = 1;
				$genes{$sequence}[$genefrags{$sequence}]{'end'} = $fragsize;
			}
			
		}
	}

	# Parse predict
	my %predictions;
	my $current_prediction = 0;
	my $current_sequence = "";
	my $current_seqname = "";
	my $in_comment_block = 0;
	my $in_prediction = 0;
	open my $mga_fh, '<', $mga_file_name;
	while (my $mga_line = <$mga_fh>){
		chomp $mga_line;
		if($mga_line =~ m/^>(.+?)\s/m){
			if($in_comment_block == 0){
				# First line of comment, sequence name
				# Check if we care about this sequence

				if(exists $genes{$1}){
					$current_sequence = $1;
					$current_seqname = substr($current_sequence, length($prefix)+1);
					if (not exists $predictions{$current_seqname}){
						$predictions{$current_seqname} = [];
					}
				} else {
					$current_sequence = "";
				}
			}
			$in_comment_block = 1;
		} else {
			# Now out of comment block
			$in_comment_block = 0;
			# If we're in a sequence we need, check and see if we're at a start or end of a fragment we need
			if($current_sequence ne ""){
				if ($mga_line =~ m/^gene_(\d+?)\s(\d+?)\s(\d+?)\s/m){
					my $gene = $1;
					my $start = $2;
					my $end = $3;
					foreach my $frag (@{$genes{$current_sequence}}){
						if ($gene eq $frag->{'start'}){
							#We're at the start of a fragment
							$current_prediction = scalar(@{$predictions{$current_seqname}});
							$predictions{$current_seqname}[$current_prediction]{'start'} = $start;
						}
						if ($gene eq $frag->{'end'}){
							#We're at the end of a fragment
							$predictions{$current_seqname}[$current_prediction]{'end'} = $end;
							$current_prediction = $current_prediction + 1;
						}
					}
				}
			}
		}
	}
	return \%predictions;
}

1;