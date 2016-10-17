#!/usr/bin/env perl
$|++; # Forces stdout flush

# Virus Hunting
# Requirements:
#  Prodigal
#  PhiSpy
#  MUMmer
#  Spine/Agent
#  MCL
#  hmmer
#  Metagene Annotator
#  MUSCLE
#  BLAST
#  BLAST+
#  VirSorter
#  BioPerl

# Inputs:
#  1) input.txt file with genome file prefixes, like:
#    NC_004663
#    GENOMEA
#  2) config.txt:
#    Basic config:
#      input_path=
#        path to the directory containing the input gff/fasta files
#      output_path=
#        path to the base directory where all output should be saved
#    Program paths: (default for all: . or in the case of executables, PATH)
#      genbank_to_seed
#      genbank_to_fasta
#      gff_to_seed
#      prodigal
#      virsorter
#      phiSpy
#      blastn
#      spine
#      agent
#      mcl
#    Program params:
#      phispy_windowsize=
#        size of window to scan genes (set default = 40)
#      phispy_threshold=
#        number of consecutive genes required to call element (set default = 20)
#      spine_percent_input=
#        Number of genome % to include as core (default = 100)
#      spine_max_distance=
#        Max distance between elements (default = 10)
#      spine_agent_min_perc_id=
#        Minimum percent ID for matching regions (default = 85)
#      spine_agent_min_size_core=
#        Minimum core lenght (default = 10)
#      spine_core_file=
#        precomputed core file from separate/previous spine run
#      spacer_fasta_file=
#        Fasta file of spacers
#      known_viral_types=
#        Fasta of known viruses to include clusters
#      virsorter_database=
#        1 for RefseqABVir only, 2 for RefseqABVir + Viromes.
#      masking_file=
#        file of regions to ignore:
#          NC_004663	x	y
#          NC_004663	z	a
#          GENOMEA		m	n
#      missed_element_padding=
#        # of base pairs between two blast hits before they are considered overlapping (default = 5)
use strict;
use Getopt::Long qw(GetOptions);
use YAML qw(LoadFile);
use File::Copy;
use File::Path qw(make_path rmtree);
use File::Spec;

#TODO only used for debugging; remove in production
use Data::Dumper;

# Constants
use define CONVERTED_INPUT_DIR => "Converted_Input_Files";
use constant OUTPUT_DIR => "Output_Files";

use VH_VirSorter;
use VH_PhiSpy;
use VH_SpineAgent;
use VH_CRISPR;
use VH_Blast;
use VH_helpers;
use VH_ReBlast;
use VH_Cluster;

### STEP 1: Parse arguments
print `python -c "import Bio"`;
# quit unless we have the correct number of args
my $num_args = $#ARGV+1;
if($num_args<1 || $ARGV[0] =~ /^-/){
	print "\nUsage: virushunting.pl input.txt [config.txt] [options]\n";
	exit;
}

print VH_helpers->current_time()."Parsing arguments: \n";
# We got at least one argument, the first one is the input text file
my $input_file = $ARGV[0];
print "\tInput file: $input_file\n";

# If we got a second argument, it's the config file
my $config_file = "";
if($num_args>1 && $ARGV[1] =~ /^[^-]/){
	$config_file = $ARGV[1];
	print "\tConfig file: $config_file\n";
}

# The rest of the arguments
# Defaults
my %params = (
	"input_path"=>".",
	"output_path"=>".",
	"genbank_to_seed"=>"genbank_to_seed.py",
	"genbank_to_fasta"=>"genbank_to_fasta.py",
	"gff_to_seed"=>"gff_to_seed.pl",
	"prodigal"=>"prodigal",
	"virsorter"=>"wrapper_phage_contigs_sorter_iPlant.pl",
	"virsorter_data_dir"=>"/data", # This is the virsorter's builtin default
	"phispy"=>"phispy.py",
	"blastn"=>"blastn",
	"makeblastdb"=>"makeblastdb",
	"spine"=>"spine.pl",
	"agent"=>"AGEnt.pl",
	"cut"=>"cut",
	"mcxload"=>"mcxload",
	"mcl"=>"mcl",
	"mcxdump"=>"mcxdump",
	"blast_to_mcl"=>"Blast_to_MCL.1.py",
	"core_genome_cluster_test"=>"core_genome_cluster_test.pl",
	"mcldump2clusters"=>"MCLdump2clusters.pl",
	"phispy_windowsize"=>40,
	"phispy_threshold"=>20,
	"spine_percent_input"=>100,
	"spine_agent_min_perc_id"=>85,
	"spine_agent_min_size_core"=>10,
	"spine_core_file"=>"",
	"spacer_fasta_file"=>"",
	"known_viral_types"=>"",
	"virsorter_database"=>1,
	"masking_file"=>"",
	"mcl_inflation"=>2.0,
	"overhang_threshold"=>100,
	"merge_threshold"=>1500,
	"clustering_parameter"=>"percent_length_aligned",
	"reblast_min_perc_id"=>90,
	"reblast_min_perc_length"=>50,
	"reblast_edge_distance"=>5,
	"reblast_distance"=>5,
	"cluster_core_congruence"=>0.66,
	"percent_id_min_core"=>85,
	"cluster_min_perc_length"=>0.8,
	"cluster_min_length"=>100,
	"cluster_min_bit_score"=>100.0,
	"cluster_core_max_distance"=>10
);

# First, read in the config file, if provided
if($config_file ne ""){
	my $config_hash = LoadFile($config_file);
	foreach my $k (keys %params) {
		$params{$k} = $config_hash->{$k} if exists $config_hash->{$k};
	}
}

# Override with command line arguments, if provided
GetOptions(	"input_path=s"=>\$params{"input_path"},
			"output_path=s"=>\$params{"output_path"},
			"genbank_to_seed=s"=>\$params{"genbank_to_seed"},
			"genbank_to_fasta=s"=>\$params{"genbank_to_fasta"},
			"gff_to_seed=s"=>\$params{"gff_to_seed"},
			"prodigal=s"=>\$params{"prodigal"},
			"virsorter=s"=>\$params{"virsorter"},
			"virsorter_data_dir=s"=>\$params{"virsorter_data_dir"},
			"phispy=s"=>\$params{"phispy"},
			"blastn=s"=>\$params{"blastn"},
			"makeblastdb=s"=>\$params{"makeblastdb"},
			"spine=s"=>\$params{"spine"},
			"agent=s"=>\$params{"agent"},
			"cut=s"=>\$params{"cut"},
			"mcxload=s"=>\$params{"mcxload"},
			"mcl=s"=>\$params{"mcl"},
			"mcxdump=s"=>\$params{"mcxdump"},
			"blast_to_mcl=s"=>\$params{"blast_to_mcl"},
			"core_genome_cluster_test=s"=>\$params{"core_genome_cluster_test"},
			"mcldump2clusters=s"=>\$params{"mcldump2clusters"},
			"phispy_windowsize=i"=>\$params{"phispy_windowsize"},
			"phispy_threshold=i"=>\$params{"phispy_threshold"},
			"spine_percent_input=i"=>\$params{"spine_percent_input"},
			"spine_agent_min_perc_id=i"=>\$params{"spine_agent_min_perc_id"},
			"spine_agent_min_size_core=i"=>\$params{"spine_agent_min_size_core"},
			"spine_core_file=s"=>\$params{"spine_core_file"},
			"spacer_fasta_file=s"=>\$params{"spacer_fasta_file"},
			"known_viral_types=s"=>\$params{"known_viral_types"},
			"virsorter_database=i"=>\$params{"virsorter_database"},
			"masking_file=s"=>\$params{"masking_file"},
			"mcl_inflation=f"=>\$params{"mcl_inflation"},
			"overhang_threshold=i"=>\$params{"overhang_threshold"},
			"merge_threshold=i"=>\$params{"merge_threshold"},
			"clustering_parameter=s"=>\$params{"clustering_parameter"},
			"reblast_min_perc_id=i"=>\$params{"reblast_min_perc_id"},
			"reblast_min_perc_length=i"=>\$params{"reblast_min_perc_length"},
			"reblast_edge_distance=i"=>\$params{"reblast_edge_distance"},
			"reblast_distance=i"=>\$params{"reblast_distance"},
			"cluster_core_congruence=f"=>\$params{"cluster_core_congruence"},
			"percent_id_min_core=i"=>\$params{"percent_id_min_core"},
			"cluster_min_perc_length=f"=>\$params{"cluster_min_perc_length"},
			"cluster_min_length=i"=>\$params{"cluster_min_length"},
			"cluster_min_bit_score=f"=>\$params{"cluster_min_bit_score"},
			"cluster_core_max_distance=i"=>\$params{"cluster_core_max_distance"}
			);

foreach my $k (sort keys %params) {
	print "\t$k: $params{$k}\n";
}
print "\n";

### STEP 2. READ input.txt

#Helper functions
sub file_of_type_exists {
	my $prefix = shift(@_);
	foreach (@_) {
		if (-f $params{"input_path"}."/".$prefix.".".$_) {
			return $prefix.".".$_;
		}
	}
	return undef;
}

# Create output directory, if necessary
make_path($params{"output_path"});

# Save each file prefix to array
my @prefixes;
open(my $fh, '<', $input_file)
	or die "Could not open input file '$input_file' $!";
while (my $row = <$fh>) {
	chomp $row;
	push @prefixes, $row;
}

# Test each file prefix for format
make_path($params{"output_path"}."/".CONVERTED_INPUT_DIR);
print VH_helpers->current_time()."Reading input files from $input_file...\n";
my @valid_prefixes;
my %contigs;
foreach (@prefixes) {
	print VH_helpers->current_time()."\t$_... ";
	my $valid_prefix_found = 0;
	# Skip if already done this step and running again
	if ( -f $params{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna" and 
		-d $params{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_" ){
		print "$_ already converted. Skipping.";
		$valid_prefix_found = 1;
		push @valid_prefixes, $_;
	} elsif (my $gbk_file_name = file_of_type_exists($_,"gbk","gb")){
		## Gbk. Generate fasta+seed
		print "$gbk_file_name found. Converting... ";

		$gbk_file_name = $params{"input_path"}."/$gbk_file_name";
		my $seed_file_name = $params{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_";
		my $fasta_file_name = $params{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna";

		# Convert gbk to seed
		`python "$params{genbank_to_seed}" "$gbk_file_name" "$seed_file_name"`;
		# Convert gbk to fasta
		`$params{genbank_to_fasta} -i "$gbk_file_name" -m genbank -s whole -a accessions -o "$_.fna"`;
		# Move fasta file into new directory
		move $params{"input_path"}."/$_.fna", $fasta_file_name;
		# TODO check formats of fasta definition lines
		$valid_prefix_found = 1;
		push @valid_prefixes, $_;
		print "Done.";
	} elsif ( (my $fasta_file_name = file_of_type_exists($_,"fna","fasta","fa")) && (my $gff_file_name = file_of_type_exists($_,"gff")) ) {
		## Fasta+gff. Generate seed
		print "$fasta_file_name and $gff_file_name found. Converting... ";

		$fasta_file_name = $params{"input_path"}."/$fasta_file_name";
		$gff_file_name = $params{"input_path"}."/$gff_file_name";
		my $seed_file_name = $params{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_";

		# Run gff_to_seed to generate seed file
		`$params{gff_to_seed} $gff_file_name $fasta_file_name $_`;
		# Move seed output into new directory
		rmtree($seed_file_name);
		move "_SEED_$_", $seed_file_name or die "Could not copy _SEED_/$_";
		# Copy fasta file into new directory
		copy $fasta_file_name, $params{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna" or die "Could not copy $_.fna";
		# TODO check formats of fasta definition lines
		$valid_prefix_found = 1;
		push @valid_prefixes, $_;
		print "Done.";
	} elsif (my $fasta_file_name = file_of_type_exists($_,"fna","fasta","fa")){
		## Fasta only. Generate gff, then generate seed
		print "$fasta_file_name found. Converting... ";

		$fasta_file_name = $params{"input_path"}."/$fasta_file_name";
		my $gff_file_name = $params{"input_path"}."/$_.gff";
		my $seed_file_name = $params{"output_path"}."/".CONVERTED_INPUT_DIR."/_SEED_$_";
		
		# Run prodigal to generate gff file
		my $prod_out = `$params{prodigal} -f gff -c -m -i $fasta_file_name -o $gff_file_name 2>&1`;
		if ($? == 0) {
			# Run gff_to_seed to generate seed file
			`$params{gff_to_seed} $gff_file_name $fasta_file_name $_`;
			# Move seed output into new directory
			rmtree $seed_file_name;
			move "_SEED_$_", $seed_file_name or die "Could not copy _SEED_/$_";
			# Copy fasta file into new directory
			copy $fasta_file_name, $params{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna" or die "Could not copy $_.fna";

			# TODO check formats of fasta definition lines
			push @valid_prefixes, $_;
			$valid_prefix_found = 1;
			print "Done."
		} else {
			print "Prodigal returned an error: $?. Skipping.";
			# Prodigal creates the gff file *before* checking to see if the fasta file is valid.
			#  So if it errors out, delete the empty gff file it creates.
			unlink $gff_file_name;
			# TODO save prodigal output to file
		}
		# TODO if cancelled halfway through prodigal run, empty/invalid gff file results
		
	} else {
		# No workable files found.
		print "No genebank or fasta files found. Skipping.";
	}
	if($valid_prefix_found == 1){
		# read fasta file to get sequence borders
		my $fasta_file_name = $params{"output_path"}."/".CONVERTED_INPUT_DIR."/$_.fna";
		open(my $fastafh, '<', $fasta_file_name) or die "Could not open fasta file.";
		my $sequence = "";
		my $seqlength = 0;
		my $seqstart = 1;
		while (my $row = <$fastafh>) {
			chomp $row;
			if(substr($row,0,1) eq '>'){
				if($sequence ne ""){
					# Save previous sequence
					$contigs{$_}{$sequence} = {'start'=>$seqstart,'end'=>$seqstart+$seqlength-1,'length'=>$seqlength};
					$seqstart = $seqstart + $seqlength;
					$seqlength = 0;
				}
				$sequence = substr($row,1);
			} else {
				$seqlength += length($row);
			}
		}
		# Save last sequence
		$contigs{$_}{$sequence} = {'start'=>$seqstart,'end'=>$seqstart+$seqlength-1,'length'=>$seqlength};
	}
	print "\n";
}
print "\n";

### STEP 2B. Read Masking File ###
my %masks;
if($params{'masking_file'} ne '' and -f $params{'masking_file'}){
	print VH_helpers->current_time()."Masking File Found. Parsing... ";
	open(my $masking_fh, '<', $params{'masking_file'});
	while(my $row = <$masking_fh>){
		chomp $row;
		if($row =~ m/^(.+?)\s+(\d+?)\s+(\d+)/m){
			my $prefix = $1;
			my $start = $2;
			my $end = $3;
			if(exists $contigs{$prefix}){
				foreach my $sequence ( keys %{ $contigs{$prefix} } ){
					if($end > $contigs{$prefix}{$sequence}{'start'} and $start < $contigs{$prefix}{$sequence}{'end'}){
						my $maskstart = $start;
						my $maskend = $end;
						if ( $start < $contigs{$prefix}{$sequence}{'start'} ){
							$maskstart = $contigs{$prefix}{$sequence}{'start'};
						}
						if ( $end > $contigs{$prefix}{$sequence}{'end'} ){
							$maskend = $contigs{$prefix}{$sequence}{'end'};
						}
						push @{$masks{$prefix}{$sequence}}, {'start'=>$maskstart, 'end'=>$maskend};
					}
				}
			}
		}
	}
	print "Done.\n";
}
print "\n";

### STEP 3A. Run VirSorter
VH_VirSorter->run(\%params,\@valid_prefixes);

### STEP 3B. Run PhiSpy
VH_PhiSpy->run(\%params,\@valid_prefixes);

### STEP 3C. Run CRISPR (optional)
my $ran_crispr = 0;
if ($params{"spacer_fasta_file"} eq ""){
	print VH_helpers->current_time()."No spacer fasta file given. Skipping CRISPR match.\n\n";
} elsif (! -f $params{"spacer_fasta_file"}) {
	print VH_helpers->current_time()."Spacer file not found. Skipping CRISPR match.\n\n";
} else {
	$ran_crispr = 1;
	VH_CRISPR->run(\%params,\@valid_prefixes);
}

### STEP 3D. Run Spine/Agent
VH_SpineAgent->run(\%params,\@valid_prefixes);

### STEP 3E. Run Blastn (optional)
my $ran_known_types = 0;
if ($params{"known_viral_types"} eq ""){
	print VH_helpers->current_time()."Known viral types not given. Skipping.\n\n";
} else {
	$ran_known_types = 1;
	VH_Blast->run(\%params,\@valid_prefixes);
}

### STEP 4. INITIAL PROCESS FOR GENOME COORDINATES OF OVERLAPPING ELEMENTS

# Finds overlaps between prediction a and list of predictions b; trims overhanging edges in b
sub overlap_exists {
	my ($a,$b) = @_;
	
	my $found_overlap = 0;
	my @predictions_to_add;
	if( not exists $a->{'used'} ){
		foreach my $prediction (@$b){
			if($a->{'start'}<=$prediction->{'start'} and $a->{'end'}>=$prediction->{'end'}){
				# Prediction entirely contained within a
				$found_overlap = 1;
				$prediction->{'used'} = 1;
			} elsif ($a->{'start'}<=$prediction->{'start'} and $a->{'end'}>=$prediction->{'start'}){
				#Prediction overhangs right side of a
				$found_overlap = 1;
				if($prediction->{'end'}-$a->{'end'}>=$params{'overhang_threshold'}){
					push @predictions_to_add, {'start'=>$a->{'end'}+1,'end'=>$prediction->{'end'}};
					$prediction->{'end'} = $a->{'end'};
				}
				$prediction->{'used'} = 1;
			} elsif ($a->{'start'}<=$prediction->{'end'} and $a->{'end'}>=$prediction->{'end'}){
				# Prediction overhangs left side of a
				$found_overlap = 1;
				if($a->{'start'}-$prediction->{'start'}>=$params{'overhang_threshold'}){
					push @predictions_to_add, {'start'=>$prediction->{'start'},'end'=>$a->{'start'}-1};
					$prediction->{'start'} = $a->{'start'};
				}
				$prediction->{'used'} = 1;
			} elsif ($a->{'start'}>$prediction->{'start'} and $a->{'end'}<$prediction->{'end'}){
				# Prediction overhangs a on both sides
				$found_overlap = 1;
				if($a->{'start'}-$prediction->{'start'}>=$params{'overhang_threshold'}){
					push @predictions_to_add, {'start'=>$prediction->{'start'},'end'=>$a->{'start'}-1};
					$prediction->{'start'} = $a->{'start'};
				}
				if($prediction->{'end'}-$a->{'end'}>=$params{'overhang_threshold'}){
					push @predictions_to_add, {'start'=>$a->{'end'}+1,'end'=>$prediction->{'end'}};
					$prediction->{'end'} = $a->{'end'};
				}
				$prediction->{'used'} = 1;
			}
		}
		push @{$b},@predictions_to_add;
	}
	return $found_overlap;
}

# Merges the predictions in the given array as much as possible
# TODO ignore masked-out predictions
sub merge_predictions {
	my $mergeable_predictions = shift;

	# Merge until there ain't no more to merge
	my $found_merges;
	my %merged_predictions;
	foreach my $sequence (keys %{$mergeable_predictions}){
		do {
			# print Dumper($mergeable_predictions{$sequence});
			$merged_predictions{$sequence} = [];
			$found_merges = 0;
			for (my $i = 0; $i < scalar(@{$mergeable_predictions->{$sequence}}); $i++){
				my $found_match = 0; 
				for (my $j=0; $j < scalar(@{$merged_predictions{$sequence}}); $j++){
					# If predictions overlap or are within merge_threshold of each other
					if($found_match == 0 
						and $mergeable_predictions->{$sequence}[$i]{'start'}<=$merged_predictions{$sequence}[$j]{'end'}+$params{'merge_threshold'}
						and $merged_predictions{$sequence}[$j]{'start'}<=$mergeable_predictions->{$sequence}[$i]{'end'}+$params{'merge_threshold'}){
						$found_match = 1;
						$found_merges = 1;
						# Expand bounds of j to contain i
						if($mergeable_predictions->{$sequence}[$i]{'start'}<$merged_predictions{$sequence}[$j]{'start'}){
							$merged_predictions{$sequence}[$j]{'start'} = $mergeable_predictions->{$sequence}[$i]{'start'};
						}
						if($mergeable_predictions->{$sequence}[$i]{'end'}>$merged_predictions{$sequence}[$j]{'end'}){
							$merged_predictions{$sequence}[$j]{'end'} = $mergeable_predictions->{$sequence}[$i]{'end'};
						}
						# TODO merge methods string intelligently
						foreach my $method (@{$mergeable_predictions->{$sequence}[$i]{'methods'}}){
							if(not grep(/^$method$/,@{$merged_predictions{$sequence}[$j]{'methods'}})){
								push @{$merged_predictions{$sequence}[$j]{'methods'}}, $method;
							}
						}
					}
				}
				if($found_match == 0){
					# If this prediction didn't overlap with any others, add it to the array
					push @{$merged_predictions{$sequence}}, $mergeable_predictions->{$sequence}[$i];
				}
			}
			$mergeable_predictions->{$sequence} = [];
			push @{$mergeable_predictions->{$sequence}},@{$merged_predictions{$sequence}};
		} while ($found_merges == 1);
	}

	return %merged_predictions;
}

sub bin_predictions {
	my $merged_predictions = shift;
	my $predictions = shift;
	my $prefix = shift;
	my $methods = shift;

	my @binned_predictions = ([],[],[],[],[]);
	# First pass: compare non-extending predictions with extending predictions
	
	foreach my $sequence (keys %{$merged_predictions}){
		for (my $j = 0; $j < scalar(@{$merged_predictions->{$sequence}}); $j++){
			for(my $i=0; $i<scalar(@{$methods}); $i++){
				if($methods->[$i]{'extend'} == 0){
					if(exists $predictions->{$prefix}{$methods->[$i]{'key'}}{$sequence}){
						if(overlap_exists($merged_predictions->{$sequence}[$j],$predictions->{$prefix}{$methods->[$j]{'key'}}{$sequence}) == 1){
							push @{$merged_predictions->{$sequence}[$j]{'methods'}},$methods->[$i]{'abbr'};
						}
					}
				}
			}
			my $bin;
			if(scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) > 2){
				$bin = 0;
			} elsif(scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) == 2) {
				$bin = 1;
			} elsif (scalar(@{$merged_predictions->{$sequence}[$j]{'methods'}}) == 1){
				for(my $k=0; $k<scalar(@{$methods}); $k++){
					if($methods->[$k]{'abbr'} eq $merged_predictions->{$sequence}[$j]{'methods'}[0]){
						$bin = $methods->[$k]{'bin'};
					}
				}
			}
			push @{$binned_predictions[$bin]}, {'sequence'=>$sequence,'methods'=>join(',',@{$merged_predictions->{$sequence}[$j]{'methods'}}),'start'=>$merged_predictions->{$sequence}[$j]{'start'},'end'=>$merged_predictions->{$sequence}[$j]{'end'}};
		}
	}

	# Second pass: compare non-extending predictions amongst themselves
	for(my $i=0;$i<scalar @{$methods}; $i++) {
		if($methods->[$i]{'extend'}==0){ # Only compare methods we haven't already merged
			foreach my $sequence (keys %{ $predictions->{$prefix}{$methods->[$i]{'key'}}}){
				foreach my $prediction ( @{ $predictions->{$prefix}{$methods->[$i]{'key'}}{$sequence} }){
					if(not exists $prediction->{'used'} and $prediction->{'start'}>0 and $prediction->{'end'}>0){ 
						my $used_methods = $methods->[$i]{'abbr'};
						my $matches = 1;
						for(my $j=$i+1; $j<scalar @{$methods}; $j++){
							if($methods->[$j]{'extend'}==0){
								if(exists $predictions->{$prefix}{$methods->[$j]{'key'}}{$sequence}){
									if(overlap_exists($prediction,$predictions->{$prefix}{$methods->[$j]{'key'}}{$sequence}) == 1){
										$matches++;
										$used_methods .= ",".$methods->[$j]{'abbr'};
									}
								}
							}
						}
						my $bin;
						if($matches > 2){
							$bin = 0;
						} elsif($matches == 2){
							$bin = 1;
						} elsif($matches == 1){
							$bin = $methods->[$i]{'bin'};
						}
						push @{$binned_predictions[$bin]}, {'sequence'=>$sequence,'methods'=>$used_methods,'start'=>$prediction->{'start'},'end'=>$prediction->{'end'}};
					}
				}
			}
		}
	}

	return \@binned_predictions;
}

print VH_helpers->current_time()."Processing output\n";
my %predictions;
my %merged_predictions;
my %binned_predictions;
my @methods = (
	{'key'=>'agent','abbr'=>'A','bin'=>2,'extend'=>1},
	{'key'=>'virsorter','abbr'=>'V','bin'=>2,'extend'=>1},
	{'key'=>'phispy','abbr'=>'P','bin'=>3,'extend'=>0}
);
if($ran_known_types){
	push @methods, {'key'=>'blast','abbr'=>'B','bin'=>2,'extend'=>1};
}
if($ran_crispr){
	push @methods, {'key'=>'crispr','abbr'=>'C','bin'=>4,'extend'=>0};
}
push @methods, {'key'=>'reblast','abbr'=>'R','bin'=>2,'extend'=>1};

foreach my $prefix (@valid_prefixes){
	print VH_helpers->current_time()."\tProcessing $prefix...\n";
	#Process VirSorter output
	print VH_helpers->current_time()."\t\tParsing VirSorter output... ";
	$predictions{$prefix}{'virsorter'} = VH_VirSorter->get_predictions(\%params,$prefix);
	print "Done.\n";
	#Process PhiSpy output
	print VH_helpers->current_time()."\t\tParsing PhiSpy output... ";
	$predictions{$prefix}{'phispy'} = VH_PhiSpy->get_predictions(\%params,$prefix);
	print "Done.\n";
	#Process CRISPR output
	if($ran_crispr == 1){
		print VH_helpers->current_time()."\t\tParsing CRISPR output... ";
		$predictions{$prefix}{'crispr'} = VH_CRISPR->get_predictions(\%params,$prefix);
		print "Done.\n";
	}
	#Process AGEnt output
	print VH_helpers->current_time()."\t\tParsing Agent output... ";
	$predictions{$prefix}{'agent'} = VH_SpineAgent->get_predictions(\%params,$prefix);
	print "Done.\n";
	#Process Known Types Blast output
	if($ran_known_types == 1){
		print VH_helpers->current_time()."\t\tParsing Known Types output... ";
		$predictions{$prefix}{'blast'} = VH_Blast->get_predictions(\%params,$prefix);
		print "Done.\n";
	}

	#Apply masking file
	if(exists $masks{$prefix}){
		print VH_helpers->current_time()."\t\tApplying mask... ";
		foreach my $sequence ( keys %{$masks{$prefix}} ){
			foreach my $mask ( @{$masks{$prefix}{$sequence}} ){
				foreach my $method ( keys %{$predictions{$prefix}} ){
					if(exists $predictions{$prefix}{$method}{$sequence}){
						my @predictions_to_add;
						foreach my $prediction (@{$predictions{$prefix}{$method}{$sequence}}){
							if($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'end'}){
								# Prediction entirely contained within mask; set to ignore prediction
								$prediction->{'start'} = -1; # Ignoring is _much_ easier than removing from the array
								$prediction->{'end'} = -1;
							} elsif ($mask->{'start'}<=$prediction->{'start'} and $mask->{'end'}>=$prediction->{'start'}) {
								# Prediction's left side overlaps mask; trim left side
								$prediction->{'start'} = $mask->{'end'}+1;
							} elsif ($mask->{'start'}<=$prediction->{'end'} and $mask->{'end'}>=$prediction->{'end'}) {
								# Prediction's right side overlaps mask; trim right side
								$prediction->{'end'} = $mask->{'start'}-1;
							} elsif ($mask->{'start'}>$prediction->{'start'} and $mask->{'end'}<$prediction->{'end'} and not exists $prediction->{'used'}) {
								# Prediction contains a masked area within it; set to ignore prediction
								$prediction->{'start'} = -1; # Ignoring is _much_ easier than removing from the array
								$prediction->{'end'} = -1;
								# Prediction contains a masked area within it; split into two predictions. Don't split ignored predictions.
								# push @predictions_to_add, {'start'=>$mask->{'end'}+1,'end'=>$prediction->{'end'}};
								# $prediction->{'end'} = $mask->{'start'}-1;
							}
						}
						push @{$predictions{$prefix}{$method}{$sequence}}, @predictions_to_add;
					}
				}
			}
		}
		print "Done.\n";
	}

	#Merge results
	print VH_helpers->current_time()."\t\tMerging predictions... ";
	my %mergeable_predictions;
	# Dump mergeable predictions into one array
	for(my$i=0; $i<scalar(@methods); $i++){
		if($methods[$i]{'extend'} == 1){ # We only want to merge Agent, Virsorter, and Blast hits
			foreach my $sequence (keys %{$predictions{$prefix}{$methods[$i]{'key'}}}){
				my @predictions_to_add;
				foreach my $prediction (@{$predictions{$prefix}{$methods[$i]{'key'}}{$sequence}}){
					push @predictions_to_add, {"methods"=>[$methods[$i]{'abbr'}],"start"=>$prediction->{'start'},"end"=>$prediction->{'end'}};
				}
				push @{$mergeable_predictions{$sequence}}, @predictions_to_add;
			}
		}
	}
	
	my %merged_predictions = merge_predictions(\%mergeable_predictions);
	print "Done.\n";

	#Compare results
	print VH_helpers->current_time()."\t\tCross-referencing predictions... ";
	$binned_predictions{$prefix} = bin_predictions(\%merged_predictions,\%predictions,$prefix,\@methods);
	print "Done.\n";
}
print "\n";

### STEP 5. RE-SCREEN PREDICTIONS AGAINST OTHER GENOMES FOR MISSED ELEMENTS ###
my $reblast_predictions = VH_ReBlast->run(\%params,\@valid_prefixes,\%binned_predictions,\%masks,\%contigs);
print VH_helpers::current_time()." Final merge...\n";
foreach my $prefix (@valid_prefixes){
	print VH_helpers::current_time()."\t$prefix... ";
	# unset 'used' on all predictions
	foreach my $method (keys %{$predictions{$prefix}}){
		foreach my $sequence (keys %{$predictions{$prefix}{$method}}){
			foreach my $prediction (@{$predictions{$prefix}{$method}{$sequence}}){
				if($prediction->{'start'}>0 and $prediction->{'end'}>0){
					delete $prediction->{'used'};
				}
			}
		}
	}
	# Merge again
	my %mergeable_predictions;
	# Dump mergeable predictions into one array
	for(my$i=0; $i<scalar(@methods); $i++){
		if($methods[$i]{'extend'} == 1){ # We only want to merge Agent, Virsorter, and Blast hits
			foreach my $sequence (keys %{$predictions{$prefix}{$methods[$i]{'key'}}}){
				my @predictions_to_add;
				foreach my $prediction (@{$predictions{$prefix}{$methods[$i]{'key'}}{$sequence}}){
					push @predictions_to_add, {"methods"=>[$methods[$i]{'abbr'}],"start"=>$prediction->{'start'},"end"=>$prediction->{'end'}};
				}
				push @{$mergeable_predictions{$sequence}}, @predictions_to_add;
			}
		}
	}
	foreach my $prediction (@{$reblast_predictions->{$prefix}}){
		if(not exists $prediction->{'masked'}){
			if(not exists $mergeable_predictions{$prediction->{'sequence'}}){
				$mergeable_predictions{$prediction->{'sequence'}} = [];
			}
			push @{$mergeable_predictions{$prediction->{'sequence'}}}, {'methods'=>['R'],'start'=>$prediction->{'start'},'end'=>$prediction->{'end'}};
		}
	}
	my %merged_predictions = merge_predictions(\%mergeable_predictions);

	# Bin again
	$binned_predictions{$prefix} = bin_predictions(\%merged_predictions,\%predictions,$prefix,\@methods);
	# push @{$binned_predictions{$prefix}[2]}, @{$reblast_predictions->{$prefix}};
	print "Done.\n";
}
print "\n";

print VH_helpers::current_time()."Saving Output...\n";
foreach my $prefix (@valid_prefixes){
	#Output to file
	print VH_helpers::current_time()."\t$prefix... ";
	sub output_bin {
		my ($fh, $bin) = @_;
		foreach my $prediction (@$bin){
			if(not exists $prediction->{'masked'}){
				print $fh $prediction->{'sequence'}."\t".$prediction->{'methods'}."\t".$prediction->{'start'}."\t".$prediction->{'end'}."\n";
			}
		}
	}
	make_path($params{'output_path'}."/".OUTPUT_DIR);
	my $output_file_name = $params{'output_path'}."/".OUTPUT_DIR."/$prefix.txt";
	open my $output_fh, '>', $output_file_name;
	print $output_fh "# Type 1: Predicted by >2 methods\n# Sequence\tmethods\tstart\tend\n";
	output_bin($output_fh,$binned_predictions{$prefix}[0]);
	print $output_fh "# Type 2: Predicted by 2 methods\n# Sequence\tmethods\tstart\tend\n";
	output_bin($output_fh,$binned_predictions{$prefix}[1]);
	print $output_fh "# Type 3: Predicted by 1 1° method\n# Sequence\tmethods\tstart\tend\n";
	output_bin($output_fh,$binned_predictions{$prefix}[2]);
	print $output_fh "# Type 4: Predicted by 1 2° method\n# Sequence\tmethods\tstart\tend\n";
	output_bin($output_fh,$binned_predictions{$prefix}[3]);
	print $output_fh "# Type 5: Protospacer match\n# Sequence\tmethods\tstart\tend\n";
	output_bin($output_fh,$binned_predictions{$prefix}[4]);
	close $output_fh;
	print "Done.\n";
}
print "\n";

### STEP 6. COMPARE/CLUSTER PREDICTIONS ###
# TODO this shouldnt be known viral types
# if($params{"known_viral_types"} eq ""){
# 	print VH_helpers->current_time()."No phage database given. Skipping clustering.\n";
# } else {
	VH_Cluster->run(\%params,\@valid_prefixes,\%binned_predictions);
# }
