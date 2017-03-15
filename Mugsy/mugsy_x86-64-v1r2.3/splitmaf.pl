#!/usr/bin/perl

#Accepts pairwise maf only
#./splitmaf.pl outputprefix < input.maf



my $qfiles = {};

my @seqs;
my @buffer;
my $currscoreline;

my $header = "##maf version=1 scoring=maf_project_simple\n";

while(my $line=<STDIN>){
    if($line =~ /^a/){
	die "Only pairwise seqs accepted" if(scalar(@seqs)>2);
	if(scalar(@seqs)>0){
	    &writemaf(\@seqs,\@buffer);
	}
	$currscoreline=$line;
	@seqs = ();
	@buffer = ();
    }
    elsif($line =~ /^s\s+([^.\s]+)/){
	push @seqs,$1;
    }
    push @buffer,$line;
}
if(scalar(@seqs)>0){
    &writemaf(\@seqs,\@buffer);
}


sub writemaf{
    my($seqs,$buffer) = @_;
    die "Invalid seqs ids $seqs->[0] $seqs->[1]" if(!defined $seqs->[0] || !defined $seqs->[1]);
    my $fh;
    if(! exists $qfiles->{$seqs->[0]}->{$seqs->[1]}){
	open $fh, "+>$ARGV[0]$seqs->[0].$seqs->[1].maf" or die "Can't open file $ARGV[0]$seqs->[0].$seqs->[1].maf: $!";
	print $fh $header;
	$qfiles->{$seqs->[0]}->{$seqs->[1]} = $fh;
	print "$ARGV[0]$seqs->[0].$seqs->[1].maf\n";
    }
    $fh = $qfiles->{$seqs->[0]}->{$seqs->[1]};
    print $fh @$buffer;
}
