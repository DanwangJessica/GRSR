#!/usr/bin/perl
#Utility for converting output of Mauve XMFA to MAF format

#USAGE: ./xmfa2maf seqs.len < aln.xmfa > aln.maf

use strict;

my $seqname;
my $start;
my $end;
my $orient;
my $seqinfo = [];
my %lens;
my $blocks = [];
my $usenum = $ARGV[1];
my $idx=0;

if($ARGV[0]){
    open(FILE,$ARGV[0]) or die "Can't open file $ARGV[0] needed for sequence lengths";
    while(my $line=<FILE>){
	chomp $line;
	my($name,$len,$newname) = split(/\s+/,$line);
	if($usenum){
	    $lens{++$idx}->{'len'} = $len;
	    if(length($newname)>0){
		$lens{$idx}->{'name'} = $newname;
	    }
	}
	elsif($name){
	    $lens{$name}->{'len'} = $len;
	    if(length($newname)>0){
		$lens{$name}->{'name'} = $newname;
	    }
	}
    }
    close FILE;
}

print "##maf version=1 scoring=mauve\n";
while(my $line=<STDIN>){
    if($line =~ /^\s*=/){
	if(defined $seqname && $start>0){
	    push @$blocks,[$seqname,$start-1,$end,$orient,$seqinfo];
	}
	if(scalar(@$blocks)>0){
	    #Convert alignment to zero start, interbase coordinates
	    print "a score=1\n";
	    foreach my $l (@$blocks){
		&printMAF(@$l);
	    }
	    print "\n";
	}
	$seqname=undef;
	$start=0;
	$seqinfo=[];
	$blocks = [];
    }
    #Format >id1:start-end orient id2
    elsif(($line =~ /^>\s+\S+\:/ && 
	   $line =~ /^>\s*(\S+)\:(\d+)-(\d+)\s+([\+\-])\s+(\S+)/)
	  || 
	  $line =~ /^>(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\S+)/){
	chomp $line;
	if(defined $seqname && $start>0){
	    push @$blocks,[$seqname,$start-1,$end,$orient,$seqinfo];
	}
	my $seqid;
	if(exists $lens{$1}){
	    $seqid = $1;
	}else{
	    if(exists $lens{$5}){
		$seqid = $5;
	    }
	    else{
		$seqid = $1;
		$lens{$1}->{'len'} = $5;
	    }
	}
	$start = $2;
	if($start>0){
	    $end = $3;
	    #XMFA format start always < end 
	    die "Invalid coordinates $start-$end" if($start>$end);
	    #Relative orientation of the alignment
	    $orient = $4;
	    my $file = $5;
	    $seqname = $seqid;
	    $seqname =~ s/^\/.*\/(\S+)/$1/;
	}
	$seqinfo=[];
    }
    else{
	if($line !~ /\#/){
	    if(defined $seqname){
		chomp $line;
		push @$seqinfo,$line if($line =~ /\S+/);
	    }
	}
    }
}

sub printMAF{
    my($id,$s,$e,$o,$str) = @_;
    die "No length specified for seq $id in $ARGV[0]" if(!exists $lens{$id});
    die "$e<$s" if($e<=$s);
    die if($o ne '+' && $o ne '-');
    my $len = $e-$s;
    $s = ($o eq '-') ? ($lens{$id}->{'len'}-$e) : $s;
    die "Bad coords $s $e $lens{$id}->{'len'}" if($s<0);
    
    my $seqlen = $lens{$id}->{'len'};
    if(exists $lens{$id}->{'name'}){
	$id = $lens{$id}->{'name'};
    }
    print "s $id $s ",$len," $o $seqlen ",join('',@$str),"\n";
}
