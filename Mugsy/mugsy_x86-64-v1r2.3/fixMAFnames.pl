#!/usr/bin/perl

use strict;


while(my $line=<STDIN>){
    if($line =~ /^s\s+(\S+)\:(\S+):\d+-\d+:\d+:[+-]:\d+/){
	if($1 eq $2){
	    $line =~ s/^s\s+(\S+)\:(\S+):\d+-\d+:\d+:[+-]:\d+/s $1.$1/;
	}
	else{
	    $line =~ s/^s\s+(\S+)\:(\S+):\d+-\d+:\d+:[+-]:\d+/s $1.$2/;
	}
    }
    elsif($line =~ /^s\s+(\S+)\:(\S+):\d+:[+-]:\d+/){
	if($1 eq $2){
	    $line =~ s/^s\s+(\S+)\:(\S+):\d+:[+-]:\d+/s $1.$1/;
	}
	else{
	    $line =~ s/^s\s+(\S+)\:(\S+):\d+:[+-]:\d+/s $1.$2/;
	}
    }
    elsif($line =~ /^s\s+\S+\s+/){
	$line =~ s/^s\s+(\S+)(\s+)/s $1.$1$2/;
    }
    print $line;
}

