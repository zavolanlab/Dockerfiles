use strict;
use DTable::DTableOps;
package DSortedTable;


# all of the functionality here can be used for sorted tables. 
# since these are sorted, the functions will run much faster.
# they need to be sorted according to the first col
#
# example for a sorted table
#        1  1524     
#     2052  3260  	
#     3280  4437    	
#     4434  4997   
#     5123  7267   
#     7302  9818    
#     9914  10828    	
#    11874  12311   	
#    12468  13016   
#--------------------------------------------------------------


# returns all rows which are in the searched range
# my @tabOut=DSortedTable::getRowsWithinRange(\@tabIn,3000,6000);
# only the first col from the input is used
#
#     result
#     3280  4437    	
#     4434  4997   
#     5123  7267 


sub getRowsWithinRange{
     my ($tabRef,$min,$max)=@_;
     my @outTab;

     if ($min>$max){my $tt=$min;$min=$max;$max=$tt;}

     if (($min<$tabRef->[0][0]) and ($max<$tabRef->[0][0]) ) {return @outTab;}
     if (($min>$tabRef->[-1][0])and ($max>$tabRef->[-1][0])){return @outTab;}


     #my($s1,$e1)=findElement($tabRef,$min);
     #my($s2,$e2)=findElement($tabRef,$max);
    
    
     #for my $i($e1..$s2){
	#push(@outTab,$tabRef->[$i]);
     #}

     return @outTab;

}

# returns the indices of the regions that are completely within the range provided with min, max

sub getRegionsWhichAreInsideRange{
    my ($featTab,$min,$max)=@_;



}



# finds an element in a sorted list
# returns the indices that surround the query
# my @l=(1,2052,3280,4434,5123,7302,9914,11874,12468);
# DSortedTable::findElementInSortedList(\@l,5130);
# returns: (4,5);
sub findElementInSortedList{
    my ($list,$query)=@_;
 
    my $l=@{$list};
    if ($query<$list->[0]){
	return (0,0);
    }
    if ($query>$list->[-1]){
        return ($l-1,$l-1);
    }

    my ($left,$mid,$right)=(0,int( ($l+1)/2 -1),$l-1);

    while ($left <$right-1){
	if ($query<$list->[$mid]){
	    ($left,$mid,$right)=($left, int( ($left+$mid)/2 ) ,$mid);
	}
	else{
	    ($left,$mid,$right)=($mid, int(($mid+$right)/2) ,$right);
	}
	#print "$left,$mid,$right\n";
    }
    return ($left,$right);
}




sub getQueryOverlapWithRegions{
    my ($tab,$from,$to)=@_;

    



}

1;
