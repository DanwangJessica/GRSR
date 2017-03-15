#!/usr/bin/perl

#./plot.pl outputprefix reforganismname ps,xll
#Eg. cat /tmp/testfilter.maf | ./plot.pl /tmp/testfilter genome2 mugsy.out > out.gp
#cat /tmp/plasmidfilter.maf | ./plot.pl /tmp/plasmidfilter AF401292 mugsy.out > out.gp
#gnuplot out.gp
#
#Requires delta files output by mugsy in outputprefix
#

use strict;

my $terminal = ($ARGV[4] =~ /ps/) ? 'postscript' : 'X11';
my($refname) = ($ARGV[1] =~ /^([^:.]+)/);
my $delta = "$ARGV[0].$refname.filt.delta";
print STDERR "Parsing $delta\n";

die "Can't find delta file" if(!-e $delta);

#
#Need to add -R -Q support for specifying the order of draft sequences

my $mummerplotcmd = "/usr/local/projects/angiuoli/developer/sangiuoli/mummer/trunk/MUMmer3.20/mummerplot -p $ARGV[0].$refname \"$delta\"";
`$mummerplotcmd 1> /dev/null 2> /dev/null`;

my $idlenlookup={};

open FILE, "$ARGV[0].$refname.gp" or die "Can't open file $ARGV[0].$refname.gp";

my $savelen;
my @xseqs;
my @yseqs;

while (my $line=<FILE>){
    if($line =~ /^set ytics/){
	$savelen=2;
    }
    if($line =~ /^set xtics/){
	$savelen=1;
    }
    if($line =~ /^set [xy]label "([^\"]+)"/){
	if($1 eq "QRY" || $1 eq "REF"){

	}
	else{
	    my $id = $1;
	    $id =~ /([^\:\.]+)[\:\.]([^\:]+)/;
	    if($1 eq $2){
		$id = $1;
	    }
	    else{
		$id = "$1.$2";
	    }

	    if($line =~ /xlabel/){
		push @xseqs,[$id,0];
	    }
	    elsif($line =~ /ylabel/){
		push @yseqs,[$id,0];
	    }
	}
    }
    if($savelen){
	my($id,$len) = ($line =~ /\"\*?([^\"]+)\"\s+(\d+)\,/);
	$id =~ /([^\:\.]+)[\:\.]([^\:]+)/;
	if($1 eq $2){
	    $id = $1;
	}
	else{
	    $id = "$1.$2";
	}
	if(defined $len && $id ne ""){
	    if($savelen==1){
		push @xseqs,[$id,$len];
	    }
	    elsif($savelen==2){
		push @yseqs,[$id,$len];
	    }
	}
    }
}

my @seqs = (@xseqs,@yseqs);

for(my $i=0;$i<@seqs;$i++){
    my($id,$len) = ($seqs[$i]->[0],$seqs[$i]->[1]);
    $idlenlookup->{$id} = $len;
}

close FILE;

open FILE1,"+>$ARGV[0].$refname.maf.fplot" or die "Can't open plot $ARGV[0].$refname.maf.fplot";
open FILE2,"+>$ARGV[0].$refname.maf.rplot" or die "Can't open plot $ARGV[0].$refname.maf.rplot";

print FILE1 "0 0 0\n";
print FILE1 "0 0 0\n";
print FILE1 "\n\n";

print FILE2 "0 0 0\n";
print FILE2 "0 0 0\n";
print FILE2 "\n\n";

my @accs = `grep ">" $delta`;
&maf2gp(\*FILE1,\*FILE2,$ARGV[1]);

my $synfile = "$ARGV[0].$refname.syn.plot";
my $reportgraphs = {};
my @graphs;
my $varreportgraphs = {};
my @vargraphs;
#Synteny blocks
if($ARGV[2]){
    open FILE,$ARGV[2] or die "Can't open output file $ARGV[2]";
    my $currgraph;
    my $currchain;
    my $name;
    while(my $line=<FILE>){
	chomp $line;
	if($line !~ /^[\s\#]/){
	    my @elts = split(/\s+/,$line);
	    if($name ne $elts[0]){
		$name = "$elts[0]";
	    }
	    my $seq = $elts[1];
	    my $start = $elts[3];
	    my $end = $elts[4];
	    $reportgraphs->{$name}->{$name}->{'seqs'}->{$seq}->{'start'} = $start;
	    $reportgraphs->{$name}->{$name}->{'seqs'}->{$seq}->{'end'} = $end;
	}
    }

    close FILE;
    @graphs = keys %$reportgraphs;
}

#Variants
if(defined $ARGV[3] && -e $ARGV[3]){
    open FILE,$ARGV[3] or die "Can't open output file $ARGV[3]";
    my $currgraph;
    my $currchain;
    my $name;
    while(my $line=<FILE>){
	chomp $line;
	if($line !~ /^[\s\#]/){
	    my @elts = split(/\s+/,$line);
	    if($name ne $elts[0]){
		$name = "$elts[0]";
	    }
	    my $seq = $elts[1];
	    my $start = $elts[3];
	    my $end = $elts[4];
	    $varreportgraphs->{$name}->{$name}->{'seqs'}->{$seq}->{'start'} = $start;
	    $varreportgraphs->{$name}->{$name}->{'seqs'}->{$seq}->{'end'} = $end;
	}
    }

    close FILE;
    @vargraphs = keys %$varreportgraphs;
}

my $first=1;
my @outlabels;
open FILE, "+>$ARGV[0].$refname.syn.plot";
foreach my $graphfile (@graphs){
    chomp $graphfile;
    foreach my $chainname (keys %{$reportgraphs->{$graphfile}}){
	my @labels = keys %{$reportgraphs->{$graphfile}->{'seqs'}};
	foreach my $x (@xseqs){
	    my $xacc = $x->[0];
	    $xacc =~ s/[\.|]/_/g;
	    if(exists $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$xacc}){
		foreach my $y (@yseqs){
		    my $yacc = $y->[0];
		    $yacc =~ s/[\.|]/_/g;
		    if(exists $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$yacc}){
			die "Can't find length for $x->[0]" if(! exists $idlenlookup->{$x->[0]});
			die "Can't find length for $y->[0]" if(! exists $idlenlookup->{$y->[0]});
			my $min0 = $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$xacc}->{'start'} += $idlenlookup->{$x->[0]};
			my $min1 = $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$yacc}->{'start'} += $idlenlookup->{$y->[0]};
			my $max0 = $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$xacc}->{'end'} += $idlenlookup->{$x->[0]};
			my $max1 = $reportgraphs->{$graphfile}->{$chainname}->{'seqs'}->{$yacc}->{'end'} += $idlenlookup->{$y->[0]};
			printf FILE ("%d %d %d #$graphfile\n",$min0,$min1,100);
			printf FILE ("%d %d %d\n",$min0,$max1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$min0,$max1,100);
			printf FILE ("%d %d %d\n",$max0,$max1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$max0,$max1,100);
			printf FILE ("%d %d %d\n",$max0,$min1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$max0,$min1,100);
			printf FILE ("%d %d %d\n",$min0,$min1,100);
			printf FILE ("\n\n");
#			push @outlabels,"set label \"$chainname\" at $min0,",$min1+(($max1-$min1)/2),"\n";
#			push @outlabels,"set label \"$chainname\" at ",$min0+(($max0-$min0)/2),",",$min1+(($max1-$min1)/2),"\n";
			push @outlabels,"set label \"$chainname\" at ",$min0+(($max0-$min0)/2),",",$min1,"\n";
#			push @outlabels,"set label \"$chainname\" at $min0,",$max1,"\n";
		    }
		}
	    }
	}
    }
}
close FILE;

open FILE, "+>$ARGV[0].$refname.var.plot";
foreach my $vargraphfile (@vargraphs){
    chomp $vargraphfile;
    foreach my $chainname (keys %{$varreportgraphs->{$vargraphfile}}){
	my @labels = keys %{$varreportgraphs->{$vargraphfile}->{'seqs'}};
	foreach my $x (@xseqs){
	    my $xacc = $x->[0];
	    $xacc =~ s/[\.|]/_/g;
	    if(exists $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$xacc}){
		foreach my $y (@yseqs){
		    my $yacc = $y->[0];
		    $yacc =~ s/[\.|]/_/g;
		    if(exists $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$yacc}){
			die "Can't find length for $x->[0]" if(! exists $idlenlookup->{$x->[0]});
			die "Can't find length for $y->[0]" if(! exists $idlenlookup->{$y->[0]});
			my $min0 = $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$xacc}->{'start'} += $idlenlookup->{$x->[0]};
			my $min1 = $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$yacc}->{'start'} += $idlenlookup->{$y->[0]};
			my $max0 = $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$xacc}->{'end'} += $idlenlookup->{$x->[0]};
			my $max1 = $varreportgraphs->{$vargraphfile}->{$chainname}->{'seqs'}->{$yacc}->{'end'} += $idlenlookup->{$y->[0]};
			printf FILE ("%d %d %d #$vargraphfile\n",$min0,$min1,100);
			printf FILE ("%d %d %d\n",$min0,$max1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$min0,$max1,100);
			printf FILE ("%d %d %d\n",$max0,$max1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$max0,$max1,100);
			printf FILE ("%d %d %d\n",$max0,$min1,100);
			printf FILE ("\n");
			printf FILE ("%d %d %d\n",$max0,$min1,100);
			printf FILE ("%d %d %d\n",$min0,$min1,100);
			printf FILE ("\n\n");
#			push @outlabels,"set label \"$chainname\" at $min0,",$min1+(($max1-$min1)/2),"\n";
#			push @outlabels,"set label \"$chainname\" at ",$min0+(($max0-$min0)/2),",",$min1+(($max1-$min1)/2),"\n";
			push @outlabels,"set label \"$chainname\" at ",$min0+(($max0-$min0)/2),",",$min1,"\n";
#			push @outlabels,"set label \"$chainname\" at $min0,",$max1,"\n";
		    }
		}
	    }
	}
    }
}
close FILE;

open FILE, "$ARGV[0].$refname.gp" or die "Can't open file $ARGV[0].$refname.gp";

my $inplot=0;
while (my $line=<FILE>){
    if($line =~ /^plot/){
	print join('',@outlabels);

	$inplot++;
    }
    elsif($inplot>0){
	if($line =~ /ls\s+2\s+$/){
	    chomp $line;
	    $line .= ", \\\n";
	}
	$inplot++;
    }
    print $line;
    if($inplot==3){
	print " \"$ARGV[0].$refname.maf.fplot\" title \"MAFFWD\" w lp ls 3, \\\n";
	print " \"$ARGV[0].$refname.maf.rplot\" title \"MAFREV\" w lp ls 4, \\\n";
	print " \"$ARGV[0].$refname.syn.plot\" title \"SYNBLOCKS\" w lp ls 5";	
	if(defined $ARGV[3]){
	    print ", \\\n";
	    print " \"$ARGV[0].$refname.var.plot\" title \"VARBLOCKS\" w lp ls 6 \n";
	}
	else{
	    print "\n";
	}

	$inplot=0;
    }


}


sub maf2gp{
    my($fh1,$fh2,$refacc)=@_;
    my $refline;
    my $x = [];	
    my $lcbnum=0;
    $refacc =~ /([^\:\.]+)[\:\.]([^\:]+)/;
    if($1 && $2 && $1 ne $2){
	$refacc = "$1.$2";
    }
    else{
	$refacc = $1;
    }
    print STDERR "Using accession: $refacc\n";
    while(my $line=<STDIN>){
	if($line =~ /^a\s+/){
	    if(scalar(@$x)>0){
		&printblock($fh1,$fh2,$x,$refacc,$lcbnum);
		$lcbnum++;
		$x = [];
	    }
	}
	else{
	    if($line =~ /^(s.+)\s+\S+/){
		push @$x,$1;
	    }
	}
    }
    &printblock($fh1,$fh2,$x,$refacc,$lcbnum);
}

sub printblock{
    my($fh1, $fh2, $scores,$ref,$lcbnum) = @_;
    my($refa,$refb,$refe,$refo,$reflen);
    my $hasref=0;
    my $refacc;
    foreach my $line (@$scores){
	my($qry) = ($line =~ /s\s+(\S+)/);
	$qry =~ /([^\:\.]+)[\:\.]([^\:]+)/;
	if($1 && $2 && $1 ne $2){
	    $qry = "$1.$2";
	}
	else{
	    $qry = $1;
	}
	if($qry =~ /^$ref/){
	    $refacc = $qry;
	    my $refoffset = $idlenlookup->{$refacc};
	    ($refa,$refb,$refe,$refo,$reflen) = ($line =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)/);
	    $refe = $refb + $refe;
	    $refe += $refoffset;
	    $refb += $refoffset;
	    $hasref=1;
	}
    }
    if($hasref==1){
	foreach my $line (@$scores){
	    my($qry) = ($line =~ /s\s+(\S+)/);
	    $qry =~ /([^\:\.]+)[\:\.]([^\:]+)/;
	    if($1 && $2 && $1 ne $2){
		$qry = "$1.$2";
	    }
	    else{
		$qry = $1;
	    }
	    if($qry ne $refacc){
		my $qryoffset = $idlenlookup->{$qry};
		#print STDERR "$qry $qryoffset\n";
		if(defined $qryoffset){
		    my($qrya,$qryb,$qrye,$qryo,$qrylen) = ($line =~ /s\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+)/);
		    $qrye = $qryb + $qrye;
		    $qryb = $qryb;
		    if($refo eq '+' && $qryo eq '+'){
			#print STDERR "$refa\t$refb\t$refe\t$refo\t$qrya\t$qryb\t$qrye\t$qryo\n";
			$qrye += $qryoffset;
			$qryb += $qryoffset;
			print $fh1 "$refb $qryb 100\n";
			print $fh1 "$refe $qrye 100\n\n\n";
			push @outlabels,"set label \"$lcbnum\" at ",$refe+100,",",$qrye,"\n";
		    }
		    elsif($refo eq '+' && $qryo eq '-'){
			$qrye = ($qrylen - $qrye);
			$qryb = ($qrylen - $qryb);
			#print STDERR "$refa\t$refb\t$refe\t$refo\t$qry\t$qryb\t$qrye\t$qryo\n";
			$qrye += $qryoffset;
			$qryb += $qryoffset;
			print $fh2 "$refe $qrye 100\n";
			print $fh2 "$refb $qryb 100\n\n\n";
			push @outlabels,"set label \"$lcbnum\" at ",$refb+100,",",$qryb,"\n";
			
		    }
		    elsif($refo eq '-' && $qryo eq '+'){
			my $refec = $reflen - $refe;
			my $refbc = $reflen - $refb;
			#print STDERR "$refa\t$refbc\t$refec\t$refo\t$qry\t$qryb\t$qrye\t$qryo\n";
			$qrye += $qryoffset;
			$qryb += $qryoffset;
			print $fh2 "$refec $qrye 100\n";
			print $fh2 "$refbc $qryb 100\n\n\n";
			push @outlabels,"set label \"$lcbnum\" at ",$refbc+100,",",$qryb,"\n";
		    }
		    elsif($refo eq '-' && $qryo eq '-'){
			my $refec = $reflen - $refe;
			my $refbc = $reflen - $refb;
			#print STDERR "$refa\t$refbc\t$refec\t$refo\t$qry\t$qryb\t$qrye\t$qryo\n";
			$qrye = $qryoffset + ($qrylen - $qrye);
			$qryb = $qryoffset + ($qrylen - $qryb);
			print $fh1 "$refec $qrye 100\n";
			print $fh1 "$refbc $qryb 100\n\n\n";
			push @outlabels,"set label \"$lcbnum\" at ",$refbc+100,",",$qryb,"\n";
		    }
		    else{
			die;
		    }
		    #print STDERR "\n";
		}
	    }
	}
    }
}
