#! /usr/bin/perl

package VH_helpers;
use POSIX qw(strftime);
use File::Path qw(rmtree);
use Data::Dumper;

sub current_time {
	return strftime("[%Y-%m-%d %H:%M:%S] ",localtime);
}

sub clean_folder {
	shift;
	my $folder_name = shift;
	my $exceptions = shift;

	opendir(my $folderfh, $folder_name);
	my @files = grep !/^\.\.?$/, readdir($folderfh);
	foreach my $f (@files){
		my $filename = $folder_name."/".$f;
		if(not grep /^$filename$/, @{$exceptions}){
			if(-d $filename){
				rmtree($filename);
			} else {
				unlink($filename);
			}
		}
	}
}

sub log {
	shift;
	my $params = shift;
	my $message = shift;
	my $level = shift;
	my $continue = shift;
	if (not defined $level){
		$level = 0;
	}
	if (not defined $continue){
		$continue = 0;
	}
	if($params->{'verbosity'}>=$level){
		print VH_helpers->current_time()." ";
		print $message;
		if($continue == 1){
			print " ";
		} else {
			print "\n";
		}
	}
}

sub log_done {
	shift;
	my $params = shift;
	my $level = shift;
	if($params->{'verbosity'}>=$level){
		print "Done.\n";
	}
}

1;