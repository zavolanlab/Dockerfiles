use strict;
use DTable::DTableOps;
package DBracket;


# usage
# my @pairs=DBracket::getPairTable("((())((())))");
# 2 3
# 1 4
# 7 8
# 6 9
# 5 10
# 0 11
#
sub getPairTable{
   my($f)=@_;
   my@UNPAIRED=();
   my@Pairs;

   for(my$j=0;$j<length($f);$j++){
      my$s=substr($f,$j,1);
      if($s eq "("){
         push(@UNPAIRED,$j);
      };
      if($s eq ")"){
         my$partner=$UNPAIRED[$#UNPAIRED]; #finds the pairing partner
         pop(@UNPAIRED);                   #removes it from the list of singles
         my$p1=$partner;
         my$p2=$j;
         push(@Pairs,[$p1,$p2]);
      }
   }

   return @Pairs;
}

#my %pairs=DBracket::getPairHash("((())((())))");
#
# note that this function returns pairs in both directions (double length of list)
# 7  -> 8
# 9  -> 6
# 11 -> 0
# 2  -> 3
# 8  -> 7
# 1  -> 4
# 4  -> 1
# 0  -> 11
# 5  -> 10
# 6  -> 9
# 10 -> 5
# 3  -> 2
#
sub getPairHash{
    my($bracket)=@_;
    my %pairHash;

    my @pairTab=getPairTable($bracket);
    
    for my $i(0..$#pairTab){
	$pairHash{$pairTab[$i][0]} = $pairTab[$i][1];
	$pairHash{$pairTab[$i][1]} = $pairTab[$i][0];
    }
    return %pairHash;
}

sub checkBracketIntegrity{
    

}

1;
