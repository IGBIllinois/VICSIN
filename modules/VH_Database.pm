#!/usr/bin/perl

# MySQL Database interface module for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_Database;

use strict;
use DBI;
use VH_helpers;
use Data::Dumper;

sub insert {
	shift;
	my $params = shift;
	my $prefixes = shift;
	my $genomes = shift;
	my $sequences = shift;
	my $predictions = shift;
	my $reblast_predictions = shift;
	my $binned_predictions = shift;
	my $clusters = shift;

	my $dbh = DBI->connect("DBI:mysql:database=".$params->{'database_name'}.";host=".$params->{'database_host'}.";port=".$params->{'database_port'}, $params->{'database_user'}, $params->{'database_pass'});

	my $login = getpwuid($<) || "unknown";

	VH_helpers->log($params,"Inserting into database... ");

	# Insert run
	my $run_stmt = $dbh->prepare('INSERT INTO runs (date, phispy_windowsize, phispy_threshold, spine_percent_input, spine_max_distance, spine_agent_min_perc_id, spine_agent_min_size_core, spine_core_file, spacer_fasta_file, known_viral_types, virsorter_database, clustering_parameter, reblast_min_perc_id, reblast_min_perc_length, reblast_distance, reblast_edge_distance, cluster_core_congruence, user) values (NOW(),?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)');
	$run_stmt->execute($params->{'phispy_windowsize'}, $params->{'phispy_threshold'}, $params->{'spine_percent_input'}, $params->{'spine_max_distance'}, $params->{'spine_agent_min_perc_id'}, $params->{'spine_agent_min_size_core'}, $params->{'spine_core_file'}, $params->{'spacer_fasta_file'}, $params->{'known_viral_types'}, $params->{'virsorter_database'}, $params->{'clustering_parameter'}, $params->{'reblast_min_perc_id'}, $params->{'reblast_min_perc_length'}, $params->{'reblast_distance'}, $params->{'reblast_edge_distance'}, $params->{'cluster_core_congruence'}, $login);
	my $run_id = $run_stmt->{'mysql_insertid'};

	my $genome_stmt = $dbh->prepare('INSERT INTO genomes (run_id,genome,length,scaffolds,genes,input_format,version,definition,accession,dblink,keywords, organism, strain) values (?,?,?,?,?,?,?,?,?,?,?,?,?)');
	my $scaffold_stmt = $dbh->prepare('INSERT INTO scaffolds (genome_id, name) values (?,?)');
	my $virsorter_stmt = $dbh->prepare('INSERT INTO virsorter (genome_id,scaffold,start,stop,category,genes_predicted,viral_hallmark_genes,viral_gene_enrich,noncaudovirales_gene_enrich,pfam_depletion,uncharacterized_gene_enrichment,strand_switch_depletion,short_gene_enrichment) values (?,?,?,?,?,?,?,?,?,?,?,?,?)');
	my $phispy_stmt = $dbh->prepare('INSERT INTO phispy (genome_id,scaffold,start,stop) values (?,?,?,?)');
	my $agent_stmt = $dbh->prepare('INSERT INTO agent (genome_id, scaffold, start, stop, gc_perc) values (?,?,?,?,?)');
	my $crispr_stmt = $dbh->prepare('INSERT INTO crispr (genome_id, query, subject, percent_id, gap, mismatch, query_start, query_stop, subject_start, subject_stop, bit, evalue) values (?,?,?,?,?,?,?,?,?,?,?,?)');
	my $blast_stmt = $dbh->prepare('INSERT INTO homology (genome_id, query, subject, percent_id, gap, mismatch, query_start, query_stop, subject_start, subject_stop, bit, evalue) values (?,?,?,?,?,?,?,?,?,?,?,?)');
	my $reblast_stmt = $dbh->prepare('INSERT INTO rescreen (genome_id, query, subject, percent_id, gap, mismatch, query_start, query_stop, subject_start, subject_stop, bit, evalue) values (?,?,?,?,?,?,?,?,?,?,?,?)');
	my $mge_stmt = $dbh->prepare('INSERT INTO mge (genome_id, scaffold, type, start, stop) values (?,?,?,?,?)');
	my $mge_virsorter_stmt = $dbh->prepare('INSERT IGNORE INTO mge_virsorter (mge_id,virsorter_id) values (?,?)'); 
	my $mge_phispy_stmt = $dbh->prepare('INSERT IGNORE INTO mge_phispy (mge_id,phispy_id) values (?,?)'); 
	my $mge_agent_stmt = $dbh->prepare('INSERT IGNORE INTO mge_agent (mge_id,agent_id) values (?,?)');
	my $mge_crispr_stmt = $dbh->prepare('INSERT IGNORE INTO mge_crispr (mge_id,crispr_id) values (?,?)'); 
	my $mge_homology_stmt = $dbh->prepare('INSERT IGNORE INTO mge_homology (mge_id,homology_id) values (?,?)'); 
	my $mge_rescreen_stmt = $dbh->prepare('INSERT IGNORE INTO mge_rescreen (mge_id,rescreen_id) values (?,?)'); 
	my $cluster_stmt = $dbh->prepare('INSERT INTO clusters (run_id, cluster_id) values (?,?)');
	my $cluster_core_stmt = $dbh->prepare('INSERT INTO cluster_core (run_id, cluster_id, mge_id, start, stop) values (?,?,?,?,?)');
	foreach my $prefix (@{$prefixes}){
		VH_helpers->log($params,"\t".$prefix."... ",1);
		# Insert genomes
		 # TODO what is `genes`?
		$genome_stmt->execute($run_id, $genomes->{$prefix}{'name'}, $genomes->{$prefix}{'length'}, $genomes->{$prefix}{'scaffolds'}, $genomes->{$prefix}{'genes'}, $genomes->{$prefix}{'format'}, $genomes->{$prefix}{'version'}, $genomes->{$prefix}{'definition'}, $genomes->{$prefix}{'accession'}, $genomes->{$prefix}{'dblink'}, $genomes->{$prefix}{'keywords'}, $genomes->{$prefix}{'organism'}, $genomes->{$prefix}{'strain'}) or die "execution failed: ".$dbh->errstr();
		my $genome_id = $genome_stmt->{'mysql_insertid'};
	
		foreach my $sequence (keys %{$sequences->{$prefix}}){
			# Insert scaffolds
			$scaffold_stmt->execute($genome_id,$sequence);

			# TODO Insert masks

			# Insert virsorter hits
			foreach my $virsorter (@{$predictions->{$prefix}{'virsorter'}{$sequence}}){
				$virsorter_stmt->execute($genome_id, $sequence, $virsorter->{'start'}, $virsorter->{'end'}, $virsorter->{'category'}, $virsorter->{'genes_predicted'}, $virsorter->{'viral_hallmark_genes'}eq''?0:$virsorter->{'viral_hallmark_genes'}, $virsorter->{'viral_gene_enrich'}, $virsorter->{'noncaudovirales_gene_enrich'}, $virsorter->{'pfam_depletion'}, $virsorter->{'uncharacterized_gene_enrichment'}, $virsorter->{'strand_switch_depletion'}, $virsorter->{'short_gene_enrichment'});
				$virsorter->{'id'} = $virsorter_stmt->{'mysql_insertid'};
			}

			# Insert phispy hits
			foreach my $phispy (@{$predictions->{$prefix}{'phispy'}{$sequence}}){
				$phispy_stmt->execute($genome_id, $sequence, $phispy->{'start'}, $phispy->{'end'});
				$phispy->{'id'} = $phispy_stmt->{'mysql_insertid'};
			}

			# Insert agent hits
			foreach my $agent (@{$predictions->{$prefix}{'agent'}{$sequence}}){
				$agent_stmt->execute($genome_id, $sequence, $agent->{'start'}, $agent->{'end'}, $agent->{'gc'});
				$agent->{'id'} = $agent_stmt->{'mysql_insertid'};
			}

			# Insert crispr hits
			foreach my $crispr (@{$predictions->{$prefix}{'crispr'}{$sequence}}){
				$crispr_stmt->execute($genome_id, $crispr->{'query'}, $sequence, $crispr->{'perc_id'}, $crispr->{'gap'}, $crispr->{'mismatch'}, $crispr->{'query_start'}, $crispr->{'query_stop'}, $crispr->{'start'}, $crispr->{'end'}, $crispr->{'bit'}, $crispr->{'evalue'});
				$crispr->{'id'} = $crispr_stmt->{'mysql_insertid'};
			}

			# Insert blast hits
			foreach my $blast (@{$predictions->{$prefix}{'blast'}{$sequence}}){
				$blast_stmt->execute($genome_id, $blast->{'query'}, $sequence, $blast->{'perc_id'}, $blast->{'gap'}, $blast->{'mismatch'}, $blast->{'query_start'}, $blast->{'query_stop'}, $blast->{'start'}, $blast->{'end'}, $blast->{'bit'}, $blast->{'evalue'}) or die "execution failed: ".$dbh->errstr();
				$blast->{'id'} = $blast_stmt->{'mysql_insertid'};
			}
		}
		# Insert reblast hits
		foreach my $reblast (@{$reblast_predictions->{$prefix}}){
			$reblast_stmt->execute($genome_id, $reblast->{'query_seq'}, $reblast->{'sequence'}, $reblast->{'perc_identity'}, $reblast->{'gap'}, $reblast->{'mismatch'}, $reblast->{'query_start'}, $reblast->{'query_stop'}, $reblast->{'start'}, $reblast->{'end'}, $reblast->{'bit'}, $reblast->{'evalue'}) or die "database insert failed: ".$dbh->errstr()."\n".Dumper($reblast);
			$reblast->{'id'} = $reblast_stmt->{'mysql_insertid'};
		}

		# Insert merged predictions
		for (my $bin=0; $bin<4; $bin++){
			foreach my $prediction (@{$binned_predictions->{$prefix}[$bin]}){
				if(not exists $prediction->{'masked'}){
					$mge_stmt->execute($genome_id,$prediction->{'sequence'},$bin+1,$prediction->{'start'},$prediction->{'end'}) or die "database insert failed: ".$dbh->errstr()."\nOffending entry: \n".Dumper($prediction);
					my $mge_id = $mge_stmt->{'mysql_insertid'};
					$prediction->{'id'} = $mge_id;
					foreach my $virsorter_index (@{$prediction->{'virsorter'}}){
						$mge_virsorter_stmt->execute($mge_id,$predictions->{$prefix}{'virsorter'}{$prediction->{'sequence'}}[$virsorter_index]{'id'});
					}
					foreach my $phispy_index (@{$prediction->{'phispy'}}){
						$mge_phispy_stmt->execute($mge_id,$predictions->{$prefix}{'phispy'}{$prediction->{'sequence'}}[$phispy_index]{'id'});
					}
					foreach my $agent_index (@{$prediction->{'agent'}}){
						$mge_agent_stmt->execute($mge_id,$predictions->{$prefix}{'agent'}{$prediction->{'sequence'}}[$agent_index]{'id'});
					}
					foreach my $crispr_index (@{$prediction->{'crispr'}}){
						$mge_crispr_stmt->execute($mge_id,$predictions->{$prefix}{'crispr'}{$prediction->{'sequence'}}[$crispr_index]{'id'});
					}
					foreach my $homology_index (@{$prediction->{'blast'}}){
						$mge_homology_stmt->execute($mge_id,$predictions->{$prefix}{'blast'}{$prediction->{'sequence'}}[$homology_index]{'id'});
					}
					foreach my $rescreen_index (@{$prediction->{'reblast'}}){
						$mge_rescreen_stmt->execute($mge_id,$reblast_predictions->{$prefix}[$rescreen_index]{'id'});
					}
				}
			}
		}
	}
	# Insert Clusters
	VH_helpers->log($params,"\tClusters...",1,0);
	for(my $cluster=0; $cluster<scalar(@{$clusters}); $cluster++){
		$cluster_stmt->execute($run_id,$cluster);
		for(my $sequence=0; $sequence<scalar(@{$clusters->[$cluster]}); $sequence++){
			my $mge = $binned_predictions->{$clusters->[$cluster][$sequence]{'prefix'}}[$clusters->[$cluster][$sequence]{'bin'}][$clusters->[$cluster][$sequence]{'index'}];
			if(exists $mge->{'id'}){
				my $mge_id = $mge->{'id'};
				foreach my $core (@{$clusters->[$cluster][$sequence]{'core'}}){
					if(not exists $core->{'stop'}){
						print $cluster."\n".Dumper($clusters->[$cluster][$sequence]);
					}
					$cluster_core_stmt->execute($run_id,$cluster,$mge_id,$core->{'start'},$core->{'stop'});
				}
			}
		}
	}

}

1;