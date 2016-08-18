#! /usr/bin/perl

package VH_helpers;
use POSIX qw(strftime);
use File::Path qw(rmtree);
use Data::Dumper;

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

1;