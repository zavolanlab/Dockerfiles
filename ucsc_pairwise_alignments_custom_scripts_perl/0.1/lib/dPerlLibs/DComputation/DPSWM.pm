use strict;
use DTable::DTableOps;
package DPSWM;


# reads a list of sequences and returns
# weight matrix filled with probabilities

sub weightMatrixFromSeqs{
    my($seqs)=@_;

    # intitalize weight matrix
    my @pswm;
    for my $i ( 0 .. 3 ) {
	for my $j ( 0 .. length($seqs->[0] )-1) {
	    $pswm[$i][$j]=0;
	}
    }

    # count the occurences of the nucleotide for each position
    for my $i(0..$#{$seqs}){
	my $sig=$seqs->[$i];
    	$sig=~tr/ACGTacgt/012301230/;
	my @nucs = split(//,$sig);
	
	for my $j(0..$#nucs){
	    $pswm[$nucs[$j]][$j]+=1;
	}

	#print "$seqs->[$i] $sig @nucs\n";
    }

    # normalize the occurence

    my ($A,$C,$G,$T);
    for my $j ( 0 .. $#{ $pswm[0] } ) {
	$A=$pswm[0][$j];
	$C=$pswm[1][$j];
	$G=$pswm[2][$j];
	$T=$pswm[3][$j];
	my $sum=$A+$C+$G+$T;
        $pswm[0][$j]=int(100*$A / $sum)/100;
	$pswm[1][$j]=int(100*$C / $sum)/100;
	$pswm[2][$j]=int(100*$G / $sum)/100;
	$pswm[3][$j]=int(100*$T / $sum)/100;

	
    }

    #DTableOps::printTableClean(\@pswm);
    return @pswm;
}


sub applyPSWMtoSeq{
    my($seq,$pswm,$thresh)=@_;
    my @hits;

    $seq =~tr/mswbdrhyvknMSWBDRHYVKN/AGAGGGATGGAAGAGGGATGGA/; 

    my $l=length($seq);
    my $t=$#{ $pswm->[0] }+1;

    my $localSeq;
    for my $i(0..$l-$t){
	$localSeq=substr($seq,$i,$t);
	$localSeq=~tr/ACGTacgt/012301230/;

	my $logOddsRatio=0;
	for my $j(0..length($localSeq)-1){
	    my $probRatio=$pswm->[substr($localSeq,$j,1)][$j]/0.25;
	    my $lnProbRatio;
	    if ($probRatio>0){
		$lnProbRatio=log($probRatio);
	    }else{$lnProbRatio=-100;}

	    $logOddsRatio+=$lnProbRatio;
	    #print "$lnProbRatio\n";
	}
	if ($logOddsRatio>$thresh){
	    push(@hits,[$i+1,$i+$t,0.01*int(100*$logOddsRatio)]);
	}
	#print "$logOddsRatio\n";

    }

    return @hits;
}




1;
