use strict;
package DMath;


#simple math ops


# returns the maximum of a list with its corresponding index
# usage:
# my ($max,$ind)=DMath::maximum((3,4,6,2,13,4));
sub maximum{
    my (@list)=@_;

    my $max=$list[0];
    my $maxIndex=0;
    for my $i (0..$#list){
	if ($list[$i]>$max){
	    $max=$list[$i];
	    $maxIndex=$i;
	}
    }
    return ($max,$maxIndex)
}


sub minimum{
    my (@list)=@_;

    my $min=$list[0];
    my $minIndex=0;
    for my $i (0..$#list){
	if ($list[$i]<$min){
	    $min=$list[$i];
	    $minIndex=$i;
	}
    }
    return ($min,$minIndex)
}

sub max{
    my (@list)=@_;
    my ($max,$ind)=DMath::maximum(@list);
    return $max;
}

sub min{
    my (@list)=@_;
    my ($min,$ind)=DMath::minimum(@list);
    return $min;
}


sub choose {
  my($n, $k) = @_;

  die "choose($n, $k) not defined" if $n < $k;

  my $tempAccum=1;
  for my $i(0..$k-1){
      $tempAccum*=($n-$i)/($k-$i);

  }
  
  return $tempAccum;
}


sub sortIndex{
    my($list)=@_;

    my @sortedIndices = sort { $list->[$a] <=> $list->[$b] } 0..$#{$list};
    my @sortedValues= @{$list}[@sortedIndices];

    return (\@sortedValues,\@sortedIndices);
}



sub test{
    print "hello\n";
    print "@INC\n";

}


1;


