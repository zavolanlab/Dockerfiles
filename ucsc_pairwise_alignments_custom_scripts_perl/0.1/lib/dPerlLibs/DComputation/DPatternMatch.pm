use strict;
package DPatternMatch;

# usage
# my $cs1='TTGAC[AT].{17}TATA[AC]T';
# my @hits=DPatternMatch::matchSeqToRegexp($seq,$cs2);

sub matchSeqToRegexp{
    my($seq,$regexp)=@_;
    my @tab;

    my($start,$end,$hit);
    while ($seq =~ /($regexp)/g) {
        $end=pos $seq;
	$start=$end-length($1)+1;
	$hit=$1;

	push(@tab,[$start,$end,$hit]);
    }  

    return @tab;

}

# returns only the number of hits
sub matchSeqToRegexpHits{
    my($seq,$regexp)=@_;

    my @tab=matchSeqToRegexp($seq,$regexp);
    return ($#tab+1);
}


sub matchSeqToRegexpOverlap{
    my($seq,$regexp)=@_;
    my @tab;

    my($start,$end,$hit);
    while ($seq =~ /(?=($regexp))/g) {
        $start=(pos $seq) +1;
	$end=$start+length($1)-1;
	$hit=$1;

	push(@tab,[$start,$end,$hit]);
    }  

    return @tab;

}



1;
