#! /usr/bin/perl

# Various helper functions for VICSIN Pipeline
# Copyright 2017 University of Illinois at Urbana-Champaign
# Author: Joe Leigh <jleigh@illinois.edu>

package VH_helpers;

use strict;
use POSIX qw(strftime);
use File::Path qw(rmtree);
use Data::Dumper;

use VICSIN;

sub current_time {
	return strftime("[%Y-%m-%d %H:%M:%S] ",localtime);
}

sub clean_folder {
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

sub file_of_type_exists {
	my $prefix = shift;
	foreach (@_) {
		if (-f (VICSIN::param("input_path")."/".$prefix.".".$_) ) {
			return $prefix.".".$_;
		}
	}
	return undef;
}

sub log {
	my $message = shift;
	my $level = shift;
	my $continue = shift;
	if (not defined $level){
		$level = 0;
	}
	if (not defined $continue){
		$continue = 0;
	}
	if(VICSIN::param('verbosity')>=$level){
		print VH_helpers::current_time()." ";
		print $message;
		if($continue == 1){
			print " ";
		} else {
			print "\n";
		}
	}
}

sub run_cmd {
	my $cmd = shift;
	VH_helpers::log("\t\t$cmd",2);
	return `$cmd`;
}

sub log_done {
	my $level = shift;
	if(not defined $level){
		$level = 0;
	}
	if(VICSIN::param('verbosity')>=$level){
		print "Done.\n";
	}
}

1;