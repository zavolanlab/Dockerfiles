use strict;
package DBisection;

# implementation of the bisection alg.
# usage:
# my $zero=DBisection::findZero($myFunc,-2,3,0.0001); (left,right,prec)
#


sub findZero{
    my($func,$l,$r,$e)=@_;

    my $m=($l+$r)/2.0;

    my $lVal=$func->evaluate($l);
    my $rVal=$func->evaluate($r);

    # check if there is a zero crossing between the starting boundaries
    if (($lVal>=0 and $rVal>=0) or ($lVal<0 and $rVal<0)){
	die "error in bisection.pm. no zero crossing between the given boundaries\n";
    }

    my $mVal=$func->evaluate($m);
 
    my $mValOld=abs($mVal);
    my $iterE=$mValOld;
    while ($iterE>$e){
	if (($mVal>0 and $lVal<0) or ($mVal<0 and $lVal>0)){
	    # go to left
	    $r=$m;
	    $m=($l+$r)/2.0;
	    $rVal=$mVal;
	    $mVal=$func->evaluate($m);
	}else{
	    # go to right
	    $l=$m;
	    $m=($l+$r)/2.0;
	    $lVal=$mVal;
	    $mVal=$func->evaluate($m);
        }	
	$iterE=abs($mVal-$mValOld);
	$mValOld=$mVal;
	#print "$l $m $r\n";
    }

    return $m;
}


1;
