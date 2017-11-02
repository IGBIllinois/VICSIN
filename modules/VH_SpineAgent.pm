#!/usr/bin/perl

# Spine/AGEnt module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_SpineAgent;

use File::Path qw(make_path);

use VICSIN;
use VH_helpers;

use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant SPINE_DIR => "Spine_Runs";
use constant AGENT_DIR => "Agent_Runs";

sub run {
	my $prefixes = shift;

	my $core_file_given = 0;
	# If no spine core file is given (or found), run spine to generate one
	# TODO If spine core file needs to be generated, we should regenerate it in case something has changed.
	if (VICSIN::param("spine_core_file") ne ""){
		VH_helpers::log("Spine core file given. Skipping spine.");
		$core_file_given = 1;
	# } elsif (VICSIN::param("spine_core_file") eq "" and not -f VICSIN::param("output_path")."/".SPINE_DIR."/spine_lock" and -f VICSIN::param("output_path")."/".SPINE_DIR."/output.backbone.fasta"){
	# 	VICSIN::param("spine_core_file") = VICSIN::param("output_path")."/".SPINE_DIR."/output.backbone.fasta";
	# 	print "\nSpine core file found. Skipping spine.\n";
	} else {
		VH_helpers::log("Starting Spine run... ");
		my $spine_input_file = "spine_input.txt";
		my $lock_file_name = VICSIN::param("output_path")."/".SPINE_DIR."/spine_lock";
		make_path(VICSIN::param("output_path")."/".SPINE_DIR);
		# Generate spine input file from prefixes
		open(my $spinefh, '>', VICSIN::param('output_path')."/".SPINE_DIR."/".$spine_input_file);
		foreach(@$prefixes){
			my $fasta_file_name = File::Spec->rel2abs( VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$_.fna");
			say $spinefh "$fasta_file_name\t$_\tfasta";
		}
		close $spinefh;
		# Create a lockfile to signify that the Spine run is in progress
		open(my $lockfh, '>', $lock_file_name);
		say $lockfh "$$";
		close $lockfh;
		
		# Run spine
		my $wdir = VICSIN::param('output_path').'/'.SPINE_DIR;
		VH_helpers::run_cmd("cd $wdir; perl ".VICSIN::param('spine')." -f $spine_input_file -p ".VICSIN::param('spine_agent_min_perc_id')." -s ".VICSIN::param('spine_agent_min_size_core')." -a ".VICSIN::param('spine_percent_input')." -g ".VICSIN::param('spine_max_distance')." -t ".VICSIN::param('num_threads')." 2>&1; cd -;");
		
		unlink $lock_file_name;
		VICSIN::setParam("spine_core_file",VICSIN::param("output_path")."/".SPINE_DIR."/output.backbone.fasta");

		print "Done.\n";
	}

	# For each file:
	VH_helpers::log("Starting AGEnt runs...");
	foreach(@$prefixes){
		my $fasta_file_name = File::Spec->rel2abs( VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$_.fna" );
		my $core_file_name = File::Spec->rel2abs( VICSIN::param("spine_core_file") );
		my $lock_file_name = VICSIN::param("output_path")."/".AGENT_DIR."/$_/${_}_agent_lock";
		my $wdir = VICSIN::param("output_path")."/".AGENT_DIR."/$_";
		make_path($wdir);

		if( $core_file_given and -f VICSIN::param("output_path")."/".AGENT_DIR."/$_/AGENT_${_}.AGENT_${_}.accessory.fasta" and not -f $lock_file_name){
			VH_helpers::log("$_ AGEnt already completed. Skipping.",1);
		} else {
			VH_helpers::log("\tRunning AGEnt for $_... ",1);
			# Create a lockfile to signify that the AGEnt run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;
			
			# Run AGEnt
			VH_helpers::run_cmd("cd $wdir; perl ".VICSIN::param('agent')." -Q F -q $fasta_file_name -R F -r $core_file_name -o AGENT_$_ -m ".VICSIN::param('spine_agent_min_perc_id')." -s ".VICSIN::param('spine_agent_min_size_core')." 2>&1; cd -;");

			unlink $lock_file_name;
		}
	}
	print "\n";
}

sub get_predictions {
	my $prefix = shift;

	my $agent_file_name = VICSIN::param("output_path")."/".AGENT_DIR."/$prefix/AGENT_$prefix.AGENT_$prefix.accessory_coords.txt";

	my %predictions;
	open my $agent_fh, '<', $agent_file_name;
	<$agent_fh>; # First line is just a header
	while(my $agent_line = <$agent_fh>){
		chomp $agent_line;
		my @agent_array = split "\t", $agent_line;
		my $seq_name = $agent_array[0];
		my $seq_length = $agent_array[1];
		my $start = $agent_array[2];
		my $end = $agent_array[3];
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
	return \%predictions;
}

1;