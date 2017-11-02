#!/usr/bin/perl

# CRISPR Match BLAST module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_CRISPR;

use File::Path qw(make_path);
use File::Basename qw(basename);
use VICSIN;
use VH_helpers;

no define CONVERTED_INPUT_DIR =>;
use constant CRISPR_DIR => "CRISPR_Runs";

sub run {
	my $prefixes = shift;

	VH_helpers::log("Starting CRISPR runs...");

	foreach(@$prefixes) {
		my $fasta_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$_.fna";
		my $wdir = VICSIN::param('output_path').'/'.CRISPR_DIR."/${_}";
		my $crispr_file_name = $wdir."/${_}_CRISPR.aln";
		my $lock_file_name = $wdir."/${_}_CRISPR_lock";
		my $db_file_name = $wdir."/${_}_db";
		my $db_name = "${_}_db";
		my $spacer_fasta_file = File::Spec->rel2abs(VICSIN::param('spacer_fasta_file'));
		my $pamproto_name = $db_name."_vs_".basename($spacer_fasta_file);
		my $pamproto_out = $wdir."/$pamproto_name.dir/$pamproto_name.extra.aln";
		
		make_path($wdir);

		# If the crispr file exists but not the lock file, the CRISPR run was already complete
		if (-f $crispr_file_name and not -f $lock_file_name){
			VH_helpers::log("$_ CRISPR already completed. Skipping.",1);
		} else {
			VH_helpers::log("\tRunning CRISPR blast for $_...",1);
			# Create a lockfile to signify that the CRISPR run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			# Create a blast database from a single genome
			VH_helpers::run_cmd(VICSIN::param('makeblastdb')." -in $fasta_file_name -dbtype nucl -parse_seqids -out $db_file_name");

			# Run PAMProtoPatternGrab_full
			File::Path::rmtree(glob($wdir."/$pamproto_name.dir"));
			VH_helpers::run_cmd("cd $wdir; ".VICSIN::param('pamprotopatterngrab')." $spacer_fasta_file $db_name; cd -");

			# Filter ...extra.aln with awk
			VH_helpers::run_cmd("awk '{if (\$18>=".VICSIN::param('crispr_match_threshold').") print}' $pamproto_out > $crispr_file_name");
			
			unlink $lock_file_name;
		}
	}
	print "\n";
}

sub get_predictions {
	my $prefix = shift(@_);

	my $aln_file_name = VICSIN::param("output_path")."/".CRISPR_DIR."/${prefix}/${prefix}_CRISPR.aln";

	my %predictions;
	open my $aln_fh, '<', $aln_file_name;
	while(my $aln_line = <$aln_fh>){
		chomp $aln_line;
		my @aln_array = split "\t", $aln_line;
		my $seq_name = $aln_array[1];

		if(not exists $predictions{$seq_name}){
			$predictions{$seq_name}=[];
		}

		my $current_prediction = scalar(@{$predictions{$seq_name}});
		$predictions{$seq_name}[$current_prediction]{'start'} = $aln_array[8];
		$predictions{$seq_name}[$current_prediction]{'end'} = $aln_array[9];
		$predictions{$seq_name}[$current_prediction]{'query'} = $aln_array[0];
		$predictions{$seq_name}[$current_prediction]{'perc_id'} = $aln_array[2];
		$predictions{$seq_name}[$current_prediction]{'gap'} = $aln_array[5];
		$predictions{$seq_name}[$current_prediction]{'mismatch'} = $aln_array[4];
		$predictions{$seq_name}[$current_prediction]{'query_start'} = $aln_array[6];
		$predictions{$seq_name}[$current_prediction]{'query_stop'} = $aln_array[7];
		$predictions{$seq_name}[$current_prediction]{'bit'} = $aln_array[11];
		$predictions{$seq_name}[$current_prediction]{'evalue'} = $aln_array[10];
		$predictions{$seq_name}[$current_prediction]{'crispr'} = [$current_prediction];
	}
	return \%predictions;
}

1;