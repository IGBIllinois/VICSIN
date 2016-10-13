#!/usr/bin/perl
## Algorithm for STEP 7. DEFINE CORE GENOME OF EACH CLUSTER
## core_genome_cluster_test.pl
## P. Degnan 
## 13 June 2016
## v0.1

## Designed for all clusters with ≥2 sequences


## Input and sort phage by size
## Names must match those in blast report

%hash=(10000,'phage1',9000,'phage2',8500,'phage3');
foreach $key (sort {$b <=> $a} keys %hash) { 
	push(@list, $hash{$key}); 
}
## Possible problem if there are phage of identical size - alternate way to get /sort data?
## Sorting by size is probably not essential.

#@list=('phage1','phage2','phage3');

## Set parameter for fraction of genomes in cluster that must have matches. 
$cutoff=2/3;
## Set parameter for merger distance between HSPs 
$distance=10;
## Set perecent identity threshold
$identity=85;
## Blast report
$blast="test.txt";

## Retrieve all relevant lines (K) from a blast report from a all-by-all blastn search
# 	1	2	3	n
# 1	-	K	K	K
# 2		-	K	K
# 3			-	K
# n				-
# Assumes bidirectional search will be the same (possibly over (in query)/under (in subject) count duplicated regions)

foreach $i (0..$#list){
	#print "$list[$i]\n";
	$next=$i+1;
	foreach $j ($next..$#list){
	
		$allmatches=`awk '\$1  ~/$list[$i]/' $blast | awk '\$2  ~/$list[$j]/' `;
		print "[$i][$j]\n$allmatches\n";
		@matches=split(/\n/,$allmatches);
		foreach $m (@matches){
			@cols=split(/\t/,$m);				
			if($cols[2] >= $identity){
				#query [0] [6] [7]
				for $i ($cols[6]..$cols[7]){
					$SEQUENCE{$cols[0]}{$i}++;
				}
				#subject [1] [8] [9]
				if($cols[8] < $cols[9]){
					for $i ($cols[8]..$cols[9]){
						$SEQUENCE{$cols[1]}{$i}++;
					}
				}else{
					for $i ($cols[9]..$cols[8]){
						$SEQUENCE{$cols[1]}{$i}++;
					}				
				}
			}
		}
	}
}

## DETERMINE COORDINATES FOR EACH GENOME
## Currently data are printed to screen. Reporting style is open ended.
## Directly load into MYSQL database?

$total=@list;

foreach $key (sort {$b <=> $a} keys %hash) { 
	$name=$hash{$key};
	print "[$key][$name]\n";
	# Reset placeholder variables for each phage
	$interval_start="";
	$interval_stop="";
	$interspace=0;
	for $i (1..$key){
		
		# print stop position when exceed distance between HSPs
		if($interspace == ($distance+1) && $interval_start ne ""){
			$interval_start="";
			print "$interval_stop\n";
			
		}
		
		# Calculate allowed $cutoff ratio of # genomes w/ shared HSP
		$ratio=($SEQUENCE{$name}{$i} +1 ) / $total;
		#if($name eq "phage3"){print "$i\t$SEQUENCE{$name}{$i}\t$ratio\n";}
		if($ratio >= $cutoff){
			if($interval_start eq ""){
				$interval_start = $i;
				#$interspace=0;
				print "$i\t";	
			}
			$interspace=0;
		}else{	
			#if just switched to from above cutoff to below, set stop position	
			if($interspace==0){
				$interval_stop=$i-1;
			}
			$interspace++;		
		}
	}
	# print out last stop position if not done already
	if($interspace == 0){
		print "$key\n";
	}elsif($interspace <= $distance){
		print "$interval_stop\n";
	}
	
}

__END__

