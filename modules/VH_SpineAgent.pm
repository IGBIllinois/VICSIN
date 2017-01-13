#!/usr/bin/perl

package VH_SpineAgent;

use File::Path qw(make_path);

use VH_helpers;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant SPINE_DIR => "Spine_Runs";
use constant AGENT_DIR => "Agent_Runs";

sub run {
	shift;
	my $params = shift;
	my $prefixes = shift;

	my $core_file_given = 0;
	# If no spine core file is given (or found), run spine to generate one
	# TODO If spine core file needs to be generated, we should regenerate it in case something has changed.
	if ($params->{"spine_core_file"} ne ""){
		VH_helpers->log($params,"Spine core file given. Skipping spine.");
		$core_file_given = 1;
	# } elsif ($params->{"spine_core_file"} eq "" and not -f $params->{"output_path"}."/".SPINE_DIR."/spine_lock" and -f $params->{"output_path"}."/".SPINE_DIR."/output.backbone.fasta"){
	# 	$params->{"spine_core_file"} = $params->{"output_path"}."/".SPINE_DIR."/output.backbone.fasta";
	# 	print "\nSpine core file found. Skipping spine.\n";
	} else {
		VH_helpers->log($params,"Starting Spine run... ");
		my $spine_input_file = "spine_input.txt";
		my $lock_file_name = $params->{"output_path"}."/".SPINE_DIR."/spine_lock";
		make_path($params->{"output_path"}."/".SPINE_DIR);
		# Generate spine input file from prefixes
		open(my $spinefh, '>', $params->{'output_path'}."/".SPINE_DIR."/".$spine_input_file);
		foreach(@$prefixes){
			my $fasta_file_name = File::Spec->rel2abs( $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna");
			say $spinefh "$fasta_file_name\t$_\tfasta";
		}
		close $spinefh;
		# Create a lockfile to signify that the Spine run is in progress
		open(my $lockfh, '>', $lock_file_name);
		say $lockfh "$$";
		close $lockfh;
		# Run spine
		my $wdir = $params->{'output_path'}.'/'.SPINE_DIR;
		`cd $wdir; perl $params->{spine} -f $spine_input_file -p $params->{spine_agent_min_perc_id} -s $params->{spine_agent_min_size_core} -a $params->{spine_percent_input} -g $params->{spine_max_distance} 2>&1; cd -`;
		unlink $lock_file_name;
		$params->{"spine_core_file"} = $params->{"output_path"}."/".SPINE_DIR."/output.backbone.fasta";

		print "Done.\n";
	}

	# For each file:
	VH_helpers->log($params,"Starting AGEnt runs...");
	foreach(@$prefixes){
		my $fasta_file_name = File::Spec->rel2abs( $params->{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna" );
		my $core_file_name = File::Spec->rel2abs( $params->{"spine_core_file"} );
		my $lock_file_name = $params->{"output_path"}."/".AGENT_DIR."/$_/${_}_agent_lock";
		my $wdir = $params->{"output_path"}."/".AGENT_DIR."/$_";
		make_path($wdir);

		if( $core_file_given and -f $params->{"output_path"}."/".AGENT_DIR."/$_/AGENT_${_}_non-core.fasta" and not -f $lock_file_name){
			VH_helpers->log($params,"$_ AGEnt already completed. Skipping.",1);
		} else {
			VH_helpers->log($params,"\tRunning AGEnt for $_... ",1);
			# Create a lockfile to signify that the AGEnt run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;
			# Run AGEnt
			`cd $wdir; perl $params->{agent} -Q F -q $fasta_file_name -R F -r $core_file_name -o AGENT_$_ -m $params->{spine_agent_min_perc_id} -s $params->{spine_agent_min_size_core} 2>&1; cd -`;
			
			unlink $lock_file_name;

			print "Done.\n";
		}
	}
	print "\n";
}

sub get_predictions {
	shift;
	my $params = shift;
	my $prefix = shift;

	my $agent_file_name = $params->{"output_path"}."/".AGENT_DIR."/$prefix/AGENT_${prefix}_non-core.key.txt";

	my %predictions;
	open my $agent_fh, '<', $agent_file_name;
	while(my $agent_line = <$agent_fh>){
		chomp $agent_line;
		if ($agent_line =~ m/^.+?\s.+?\s.+?\s(.+?)\s(\d+?)\s([\d\?]+?)\s([\d\?]+?)$/){
			my @agent_array = split "\t", $agent_line;
			my $seq_name = $1;
			my $seq_length = $2;
			my $start = $3;
			my $end = $4;
			if($start eq "?"){
				$start = 1;
			}
			if($end eq "?"){
				$end = $seq_length;
			}

			if(not exists $predictions{$seq_name}){
				$predictions{$seq_name}=[];
			}
			my $current_prediction = scalar(@{$predictions{$seq_name}});
			$predictions{$seq_name}[$current_prediction]{'start'} = $start;
			$predictions{$seq_name}[$current_prediction]{'end'} = $end;
			$predictions{$seq_name}[$current_prediction]{'gc'} = $agent_array[2];
			$predictions{$seq_name}[$current_prediction]{'agent'} = [$current_prediction];
		}
	}
	return \%predictions;
}

1;