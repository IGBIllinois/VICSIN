#!/usr/bin/perl

# VirSorter module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_VirSorter;

use strict;
use File::Path qw(make_path);
use File::Copy qw(mv);
use VICSIN;
use VH_helpers;
use Data::Dumper;

no define CONVERTED_INPUT_DIR =>;
use constant VIRSORTER_DIR => "Virsorter_Runs";

sub run {
	my $prefixes = shift;

	VH_helpers::log("Starting VirSorter runs...");
	foreach(@$prefixes) {
		my $fasta_file_name =  File::Spec->rel2abs( VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$_.fna" );
		my $wdir = VICSIN::param("output_path")."/".VIRSORTER_DIR."/$_";
		my $data_dir = VICSIN::param("virsorter_data_dir");

		my $lock_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/${_}/${_}_VirSorter_lock";
		my $csv_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/${_}/${_}_global-phage-signal.csv";
		my $mga_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/${_}/fasta/${_}_mga_final.predict";
		my $mga_dest_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/${_}/${_}_mga_final.predict";

		# If the virsorter files exist but not the lock file, virsorter previously completed
		if ( -f $csv_file_name and -f $mga_dest_file_name and not -f $lock_file_name ){ 
			VH_helpers::log("\t$_ VirSorter already completed. Skipping.",1);
		} else {
			VH_helpers::log("\tRunning VirSorter for $_... ",1);
			make_path($wdir);
			# Create a lockfile to signify that the VirSorter run is in progress
			open(my $lockfh, '>', $lock_file_name);
			say $lockfh "$$";
			close $lockfh;

			VH_helpers::run_cmd("cd $wdir; ".VICSIN::param('virsorter')." -d $_ --fna $fasta_file_name --db ".VICSIN::param('virsorter_database')." --data-dir $data_dir --ncpu ".VICSIN::param('num_threads')." 2>&1; cd -;");
			
			# Move mga file to its final destination
			mv($mga_file_name,$mga_dest_file_name);
			unlink($lock_file_name);
			VH_helpers::clean_folder($wdir,[$csv_file_name,$mga_dest_file_name,$wdir."/logs"]);
		}
	}
	print "\n";
}

sub get_predictions {
	my $prefix = shift;

	# Workaround for VirSorter's renaming of all the sequences. Thanks, VirSorter.
	my %unescapedSeqnames;

	my $fasta_file_name = VICSIN::param("output_path")."/".CONVERTED_INPUT_DIR."/$prefix.fna";
	open(my $fasta_fh, '<', $fasta_file_name);
	while(my $row = <$fasta_fh>){
		chomp $row;
		if($row =~ /^>(.+)$/){
			my $seqname = $1;
			my $escapedname = $seqname =~ s/[\/\.,\|\s?!\*%]/_/rg;
			# my $escapedname = $1;
			$seqname =~ /^([^\s]*)\s/;
			$seqname = $1;
			$unescapedSeqnames{$escapedname} = $seqname;
		}
	}

	my $csv_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/$prefix/${prefix}_global-phage-signal.csv";
	my $mga_file_name = VICSIN::param("output_path")."/".VIRSORTER_DIR."/$prefix/${prefix}_mga_final.predict";

	# Parse csv
	my %genes;
	my %genefrags;
	open my $csv_fh, '<', $csv_file_name;
	while (my $csv_line = <$csv_fh>){
		chomp $csv_line;
		if($csv_line =~ m/^([^#>]+?),.*?,(.+?),(.+?),/){
			my @csv_array = split ',', $csv_line;
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

			$genes{$sequence}[$genefrags{$sequence}]{'category'} = 							$csv_array[4];
			$genes{$sequence}[$genefrags{$sequence}]{'genes_predicted'} = 					$csv_array[3];
			$genes{$sequence}[$genefrags{$sequence}]{'viral_hallmark_genes'} = 				$csv_array[5];
			$genes{$sequence}[$genefrags{$sequence}]{'viral_gene_enrich'} = 				$csv_array[6];
			$genes{$sequence}[$genefrags{$sequence}]{'noncaudovirales_gene_enrich'} = 		$csv_array[7];
			$genes{$sequence}[$genefrags{$sequence}]{'pfam_depletion'} = 					$csv_array[8];
			$genes{$sequence}[$genefrags{$sequence}]{'uncharacterized_gene_enrichment'} = 	$csv_array[9];
			$genes{$sequence}[$genefrags{$sequence}]{'strand_switch_depletion'} = 			$csv_array[10];
			$genes{$sequence}[$genefrags{$sequence}]{'short_gene_enrichment'} = 			$csv_array[11];
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
					if(exists $unescapedSeqnames{$current_seqname}){
						$current_seqname = $unescapedSeqnames{$current_seqname};
					}
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

							$predictions{$current_seqname}[$current_prediction]{'category'} = 							$frag->{'category'};
							$predictions{$current_seqname}[$current_prediction]{'genes_predicted'} = 					$frag->{'genes_predicted'};
							$predictions{$current_seqname}[$current_prediction]{'viral_hallmark_genes'} = 				$frag->{'viral_hallmark_genes'};
							$predictions{$current_seqname}[$current_prediction]{'viral_gene_enrich'} = 					$frag->{'viral_gene_enrich'};
							$predictions{$current_seqname}[$current_prediction]{'noncaudovirales_gene_enrich'} = 		$frag->{'noncaudovirales_gene_enrich'};
							$predictions{$current_seqname}[$current_prediction]{'pfam_depletion'} = 					$frag->{'pfam_depletion'};
							$predictions{$current_seqname}[$current_prediction]{'uncharacterized_gene_enrichment'} = 	$frag->{'uncharacterized_gene_enrichment'};
							$predictions{$current_seqname}[$current_prediction]{'strand_switch_depletion'} = 			$frag->{'strand_switch_depletion'};
							$predictions{$current_seqname}[$current_prediction]{'short_gene_enrichment'} = 				$frag->{'short_gene_enrichment'};
							$predictions{$current_seqname}[$current_prediction]{'virsorter'} = 							[$current_prediction];
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
