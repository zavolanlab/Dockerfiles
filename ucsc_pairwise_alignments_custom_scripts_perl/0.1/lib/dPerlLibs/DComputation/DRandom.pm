use strict;
package DRandom;




# computes a random permutation and returns an array of indices
# input: number of elements
# output: indices (0..$elements-1) !!!
sub perm{
    my ($elements)=@_;

    # fisher_yates_shuffle

    my $array = [0..$elements-1];
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }

    return @{$array};
}


sub select{
    my ($array,$selector)=@_;
    my @outArray;

    for my $i(0..$#{$selector}){
	push(@outArray,$array->[$selector->[$i]]);
    }

    return @outArray;
}


1;
