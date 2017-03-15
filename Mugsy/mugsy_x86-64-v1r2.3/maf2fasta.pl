#!/usr/bin/perl
#Convert MAF to FASTA
#Optionally only convert blocks that contain label 
#./maf2fasta.pl [label] < maf > fasta

use strict;

my $currscore;
my $currlabel;
my $currcoord;
my $currorient;
my $saveblock=0;
my @matches;

my @blocks;

while(my $line=<STDIN>){
    if($line =~ /a\s+score=([\d\.\-]+)/){
	if($saveblock>0){
	    my @nmatches = @matches;
	    push @blocks,[$currscore,$currlabel,$currorient,$currcoord,\@nmatches];
	}
	($currscore) = ($line =~ /a\s+score=([\d\.\-]+)/);
	($currlabel) = ($line =~ /label=(\d+)/);
	@matches=();
    }
    elsif($line =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)/){
	my $accession = $1;
	my $start = $2;
	my $len = $3;
	my $orientation = $4;
	my $seqlength = $5;
	my $seq = $6;
	if($accession =~ /([^\.]+)\.(\S+)/){
	    
	}
	else{
	    die if($accession =~ /\./);
	    #$accession = "$accession.$accession";
	}
	push @matches,[$accession,$start,$len,$orientation,$seqlength,$seq];
	$saveblock=1;
    }
    else{
	
    }
}
if($saveblock>0){
    my @nmatches = @matches;
    push @blocks,[$currscore,$currlabel,$currorient,$currcoord,\@nmatches];
}
foreach my $block (sort {$a->[3] <=> $b->[3]} @blocks){
    if($ARGV[0]){
	if($block->[1] eq $ARGV[0]){
	    &printFASTA(@$block) ;
	}
    }
    else{
	&printFASTA(@$block) ;
    }
}

sub printFASTA{
    my($score,$label,$orient,$coord,$matches) = @_;
    foreach my $m (@$matches){
	#print ">$m->[0].$label score=$score $m->[1] $m->[2] $m->[3] $m->[4]\n";	
	print ">$m->[0] $m->[1] $m->[2] $m->[3] $m->[4]\n";
	for(my $i=0;$i<length($m->[5]);$i+=60){
	    print substr($m->[5],$i,60),"\n";
	}
    }
    print "=\n";
}
