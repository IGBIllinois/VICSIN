# VICSIN

Viral, Integrative, and Conjugative Sequence Identification and Networking

## Dependencies

VICSIN is built mainly in Perl. Perl version 5.16 and up have been tested and are supported. Older versions may work, but have not been tested. The following Perl packages are required:

* YAML
* Bio::SeqIO

Additionally, some components are written in Python2.7. BioPython is required.

VICSIN uses the following software packages in its pipeline. The pipeline has been tested using the versions given.

* Prodigal >2.6.2
* VirSorter 1.0.3
* PhiSpy 2.3
* BLAST 2.2.26
* BLAST+ >2.3.0
* Spine 0.2.1
* AGEnt 0.2.1
* MCL >12.068

## Usage

vicsin <input.txt> [<config.yml>] [options]

* input.txt: a file listing the names of all input fasta/genbank files to use. One name per line, plain text.
* config.yml: an optional config file, using yaml syntax, which can contain any of the command line options.
* options:
    * input_path: Path in which all the input files listed in input.txt can be found. Default: current working directory.
	* output_path: Path in which all output files will be written. Directory will be created if not found. Default: current working directory.
	* genbank_to_seed: Path to genbank_to_seed.py, if not in PATH
	* genbank_to_fasta: Path to genbank_to_fasta.py, if not in PATH
	* core_genome: Path to core_genome.py, if not in PATH
	* gff_to_seed: Path to gff_to_seed.pl, if not in PATH
	* prodigal: Path to prodigal, if not in PATH
	* virsorter: Path to wrapper_phage_contigs_sorter_iPlant.pl, if not in PATH
	* virsorter_data_dir: Path to VirSorter's data directory (default = "/data")
	* phispy: Path to phispy.py, if not in PATH
	* blastn: Path to blastn, if not in PATH
	* makeblastdb: Path to makeblastdb, if not in PATH
	* spine: Path to spine.pl, if not in PATH
	* agent: Path to AGEnt.pl, if not in PATH
	* cut: Path to cut, if not in PATH
	* mcxload: Path to mcxload, if not in PATH
	* mcl: Path to mcl, if not in PATH
	* mcxdump: Path to mcxdump, if not in PATH
	* blast_to_mcl: Path to Blast_to_MCL.1.py, if not in PATH
	* mcldump2clusters: Path to MCLdump2clusters.pl, if not in PATH
	* num_threads: Number of threads to use when running VirSorter and Spine
	* phispy_windowsize: Size of window to scan genes (default = 40)
	* phispy_threshold: Number of consecutive genes required to call element (default = 20)
	* spine_percent_input: Number of genome % to include as core (default = 100)
	* spine_max_distance: Max distance between elements (default = 10)
	* spine_agent_min_perc_id: Minimum percent ID for matching regions (default = 85)
	* spine_agent_min_size_core: Minimum core length (default = 10)
	* spine_core_file: Path to precomputed core file from separate/previous spine run
	* spacer_fasta_file: Path to fasta file of spacers for CRISPR match
	* known_viral_types: Path to fasta file of known viruses to include in clusters
	* virsorter_database: 1 for RefseqABVir only, 2 for RefseqABVir + Viromes.
	* masking_file: Path to tab delimited file of regions to ignore
	* mcl_inflation: MCL inflation cutoff for isolating clusters (default = 20)
	* overhang_threshold: Threshold for how many bp one prediction must overlap another before they are counted as overlapping a hit from another method (default = 100)
	* merge_threshold: Threshold for how many bp apart to predictions from the same method must be before they are merged into one prediction (default = 1500)
	* reblast_min_perc_id: Minimum percent id threshold during reBLAST step (default = 90)
	* reblast_min_perc_length: Minimum percent length threshold during reBLAST step (default = 50)
	* reblast_edge_distance: Threshold for distance from contig edge before reBLAST hit snaps to edge (default = 5)
	* reblast_distance: Threshold for distance between reBLAST hits before they merge (default = 5)
	* reblast_min_contig_length: Minimum size for reBLAST hits (default = 10000)
	* clustering_parameter: Metric to use during clustering; either total_length_aligned, total_bit_score, or percent_length_aligned (default = percent_length_aligned)
	* cluster_core_congruence: Threshold for percent of matching clustered predictions required to define core genome (default = 0.66)
	* percent_id_min_core: Threshold for percent id in core genome (default = 85)
	* cluster_min_perc_length: Threshold for minimum percent length to be included in cluster, if using percent_length_aligned (default = 0.5)
	* cluster_min_length: Threshold for minimum length to be included in cluster, if using total_length_aligned (default = 100)
	* cluster_min_bit_score: Threshold for minimum bit score to be included in cluster, if using total_bit_score (default = 100.0)
	* cluster_core_max_distance: Distance threshold for defining core genome (default = 10)
	* cluster_size_threshold: Threshold between "small" and "large" predictions when clustering (default = 2000)
	* use_database: Whether or not to log predictions in a MySQL database (default = false)
	* database_host: IP/URL of MySQL database host
	* database_name: MySQL database name
	* database_port: MySQL port (default = 3306)
	* database_user: MySQL username
	* database_pass: MySQL password
	* verbosity: Detail level of output log, 0-2 (default = 0)
	