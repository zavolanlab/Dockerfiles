use strict;
use DTable::DTableOps;
use DComputation::DMath;
package DRnaCut;


sub findOptimalCut{
    my ($PPM,$Lmin,$Lmax,$Rmin,$Rmax)=@_;

    ($Lmin,$Lmax,$Rmin,$Rmax)=($Lmin-1,$Lmax-1,$Rmin-1,$Rmax-1); 
    # find the substring ($l..$r) which folds mainly with itself
    
    # algorithm:
    # create four temporary matrices with accums. o(n^2)
    
    my @horAccumL;
    my @horAccumR;
    my @vertAccumUL;
    my @vertAccumUR;
    
    for my $i(0..$#{$PPM}){
	my $tAccum=0;
	for my $j(0..$#{$PPM->[$i]}){
	    $tAccum+=$PPM->[$i][$j];
	    $horAccumL[$i][$j]=$tAccum;
	}
    }

    for (my $i=0; $i<=$#{$PPM}; $i++){
	my $tAccum=0;
	for (my $j=$#{$PPM->[$i]}; $j>=0; $j--){
	    $tAccum+=$PPM->[$i][$j];
	    $horAccumR[$i][$j]=$tAccum;
	}
    }

    for (my $j=0; $j<=$#{$horAccumL[0]}; $j++){
	my $tAccum=0;
	for (my $i=0; $i<=$#horAccumL; $i++){
	    $tAccum+=$horAccumL[$i][$j];
	    $vertAccumUL[$i][$j]=$tAccum;
	}
    }

    for (my $j=$#{$horAccumR[0]}; $j>=0; $j--){
	my $tAccum=0;
	for (my $i=0; $i<=$#horAccumR; $i++){
	    $tAccum+=$horAccumR[$i][$j];
	    $vertAccumUR[$i][$j]=$tAccum;
	}
    }


    my $slider=DMath::min($Lmax-$Lmin,$Rmax-$Rmin);
    # compute the cross edges
    my $bestL=0;
    my $bestR=0;
    my $tempMin=1000000;
    for my $s(0..$slider){
	my $l=$Lmin+$s;
	my $r=$Rmin+$s;
	my $lWin=$vertAccumUL[$r][$l]-$vertAccumUL[$l][$l];
	my $rWin=$vertAccumUR[$r][$r]-$vertAccumUR[$l][$r];
	my $tot=$lWin+$rWin;
	if ($tot<$tempMin){
	   $tempMin=$tot;
	   $bestL=$l;
	   $bestR=$r;
	}
	#print "$l $r $tot\n";
    }


    #DTableOps::printMatrix(\@horAccumL);
    #DTableOps::printMatrix(\@vertAccumUL);
    #DTableOps::printMatrix(\@vertAccumUR);

    #print "$Lmin,$Lmax,$Rmin,$Rmax $slider\n";
    #print "$bestL $bestR \n";
    return ($bestL+1,$bestR+1);
}

1;
