use strict;
use DTable::DSortedTable;
package DSet;


# set operations
# a set is defined as a hash with val 1

sub mergeSets{
    my($set1,$set2)=@_;

    my %resSet;
    for my $k(keys %{$set1}){
	$resSet{$k}=1;
    }
    for my $k(keys %{$set2}){
	$resSet{$k}=1;
    }

    return \%resSet;
}


# returns the three groups after intersection
sub intersect{
    my($set1,$set2)=@_;
    
    my %left;
    my %right;
    my %mid;

    # check 1 in 2
    for my $k1(keys %{$set1}){
	if (defined $set2->{$k1}){
	    # element in intersection
	    $mid{$k1}=1;
	}else{
	    # element is left
	    $left{$k1}=1;
	}

    }
 

    for my $k2(keys %{$set2}){
	if (!defined $mid{$k2}){
	    $right{$k2}=1;
	}
    }

    return (\%left,\%mid,\%right);
}

# returns all the intersecting parts of 3 sets:
# set1
# set2
# set3
# set1+set2
# set2+set3
# set1+set3
# set1+set2+set3
sub intersect3{
    my($set1,$set2,$set3)=@_;

    my ($s1n2,$s1and2,$sn12)=intersect($set1,$set2);
    my ($s2n3,$s2and3,$sn23)=intersect($set2,$set3);
    my ($s1n3,$s1and3,$sn13)=intersect($set1,$set3);

    my ($s2and3n1,$s1and2and3,$x002)=intersect($s2and3,$set1);
    
    my ($s1n2n3,$s1and3n2,$s3n1n2)=intersect($s1n2,$sn23);

    my ($x003,$s1and2n3,$s2n1n3)=intersect($s1n3,$s2n3);


    return ($s1n2n3,$s2n1n3,$s3n1n2,$s1and2n3,$s2and3n1,$s1and3n2,$s1and2and3);


}

# computes a weighted intersection
# both input sets need to have a weight associated to each element
# computed are the expected number of elements:
# n00 that are in none of the two sets
# n10 that are in set1 only
# n01 that are in set2 only
# n11 that are in both sets
#
sub weightedIntersection{
    my ($set1,$set2,$allElements)=@_;

    my($n00,$n10,$n01,$n11)=(0,0,0,0);

    for my $el(keys %{$allElements}){
	my $wSet1=0;
	my $wSet2=0;

	if(defined $set1->{$el}){$wSet1=$set1->{$el};}
	if(defined $set2->{$el}){$wSet2=$set2->{$el};}

	$n00+=(1-$wSet1)*(1-$wSet2);
	$n10+=$wSet1*(1-$wSet2);
	$n01+=(1-$wSet1)*$wSet2;
	$n11+=$wSet1*$wSet2;
    }

    return ($n00,$n01,$n10,$n11);
}



sub signOverlap2SetsConstraint{
    my($set1Ids,$set2Ids,$allIds,$constraint,$rndVers)=@_;

    # find overlap between set1 and set2
    my $inSet1andSet2=0;
    foreach(keys %{$set1Ids}){if(defined $set2Ids->{$_}){$inSet1andSet2++;}}


    # sort $constraint
    my @utrLengthsSorted;
    my @utrLengthsSortedVec;
    foreach(sort {$constraint->{$a} <=> $constraint->{$b}} keys %{$constraint}){
	push(@utrLengthsSorted,[$constraint->{$_},$_]);
	push(@utrLengthsSortedVec,$constraint->{$_});
    }


    # randomize set2
    my $meanTot=0;
    my $varTot=0;
    for my $origId(sort keys %{$set2Ids}){
	if(defined $constraint->{$origId}){
	    my $origUtrLen=$constraint->{$origId};
	
	    
	    # get random genes with approx same utr length
	    # find the utr with the same length as the target utr
	    my ($left,$right)=DSortedTable::findElementInSortedList(\@utrLengthsSortedVec,$origUtrLen);
	    
	    # start at left and get closest utrs in both directions
	    my $bL=$left;
	    my $bR=$left+1;
	    my @closestUTRs;
	    for(my $i=0;$i<$rndVers+1;$i++){
		if(($bL>=0) and ($bR<@utrLengthsSorted)){
		    # both pointers within boundaries
		    my $v1=$utrLengthsSorted[$bL][0];
		    my $v2=$utrLengthsSorted[$bR][0];
		    
		    if(abs($v1-$origUtrLen)<abs($v2-$origUtrLen)){
			# left wins
			if($utrLengthsSorted[$bL][1] ne $origId){
			    push(@closestUTRs,[$utrLengthsSorted[$bL][1],$utrLengthsSorted[$bL][0]]);
			}
			$bL--;
		    }else{
			# right wins
			if($utrLengthsSorted[$bR][1] ne $origId){
			    push(@closestUTRs,[$utrLengthsSorted[$bR][1],$utrLengthsSorted[$bR][0]]);
			}
			$bR++;
		    }
		}elsif(($bL>=0) and ($bR>=@utrLengthsSorted)){
		    # left pointer within boundaries, left wins
		    if($utrLengthsSorted[$bL][1] ne $origId){
			push(@closestUTRs,[$utrLengthsSorted[$bL][1],$utrLengthsSorted[$bL][0]]);
		    }
		    $bL--;
		}elsif(($bL<0) and ($bR<@utrLengthsSorted)){
		    # right pointer within boundaries,right wins
		    if($utrLengthsSorted[$bR][1] ne $origId){
			push(@closestUTRs,[$utrLengthsSorted[$bR][1],$utrLengthsSorted[$bR][0]]);
		    }
		    $bR++;	
		}else{die "not possible 2342\n";}
	    }
	    
	    # compute the fraction of overlap with set1
	    my $pOverl=0;
	    for(my $i=0;$i<$rndVers;$i++){
		if(defined $set1Ids->{$closestUTRs[$i]->[0]}){
		    $pOverl++;
		}
	    }
	    $pOverl/=$rndVers;
	    #print "$origId $origUtrLen $pOverl\n";
	    
	    my $mean=$pOverl;
	    my $var=$pOverl*(1-$pOverl);
	    
	    $meanTot+=$mean;
	    $varTot+=$var;
	    
	}else{die "error in DSet library. no utr annotation for $origId\n";}
	
    }
    my $stdevTot=sqrt($varTot);

    # real overlap, expected overlap, standard dev
    return ($inSet1andSet2,$meanTot,$stdevTot);
}


1;
