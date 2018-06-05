use strict;
use DTable::DTableOps;
package DOverlappingSequences;


# compute the clusters of overlapping sequences
# every sequence is determined by its start and end point
# if two or more sequences overlap, then a cluster is built
# 2 input parameters:
# 1. table with the ranges of all the seqences, e.g 
# ([ 1,30],
#  [80,99],
#  [69,90],
#  [45,80])
# 2. match function callback. The function which decides
#    which regions should merge and which should not.
#    these may be defined by the user. see 
#    - doRegionsOverlap($r1Start,$r1End,$r2Start,$r2End);
#    - isOneRegionContainedInTheOther($r1Start,$r1End,$r2Start,$r2End);
#
# 2 output parameters:
#
# 1. table with all the ranges of the clusters and the ids of sequences they contain
# ([ 1,30, 0],
#  [45,99, 1 , 2 ,3]) 
#
# 2. table with the original ranges and the correspondig cluster number
# ([ 1,30, 0],
#  [80,99, 1],
#  [69,90, 1],
#  [45,80, 1])
#
# usage: 
# my($clusterTableRef,$seqsClusterTableRef)=DOverlappingSequences::mapSeqsIntoClusters(\@concatRegionIndices,
#    \&DOverlappingSequences::isOneRegionContainedInTheOther);
#
# the complexity of the algorithm is o(n log n) where n is the number of sequences

sub mapSeqsIntoClusters{
    my ($rangesTabRef,$callback)=@_;
    my @rangesTab=@{$rangesTabRef};

    my @clusterTable;
    my @seqsClusterTable;

    # if input is empty return
    if ($#{$rangesTabRef}==-1){return(\@clusterTable,\@seqsClusterTable);}

    # sort the seqences according to the start site

    my ($sortedTableRef,$sortedIndicesRef)=DTableOps::sortTableForColumn(\@rangesTab,0); 
    my @sortedTable=@{$sortedTableRef};
    my @sortedIndices=@{$sortedIndicesRef};


    # start building up a cluster table:
    # iterate linearly through the sequences and look for overlappings with the last cluster
    # if overlap, then put the sequence into this cluster
    # else generate a new cluster and push the sequence
    
    # initialize the cluster table with one cluster containing the first sequence
    my @rangeOfFirstSequence=($sortedTable[0][0],$sortedTable[0][1] );   #(rStart,rEnd)
    $clusterTable[0]=\@rangeOfFirstSequence;

    my $seqCounter=0;
    foreach (@sortedTable){

	# check if the sequence overlaps with the boundries of the last cluster 
	if (&$callback(@{$_}[0], @{$_}[1], ${$clusterTable[-1]}[0], ${$clusterTable[-1]}[1] )){
	    # insert the sequence into this cluster and update the boundaries of the cluster

	    my ($cluStart,$cluEnd)=spanOfTwoRegions( @{$_}[0], @{$_}[1], ${$clusterTable[-1]}[0], ${$clusterTable[-1]}[1] );

	    # update the boundaries of the current cluster
	    ${$clusterTable[-1]}[0]=$cluStart;
	    ${$clusterTable[-1]}[1]=$cluEnd;
	    
	    push( @{$clusterTable[-1]},$seqCounter);
	    #print "cluster collision with cluster nr $#clusterTable at @{$_}  \n";

	}else{
	    # create a new cluster for this sequence
	    my @newCluster=( @{$_}[0], @{$_}[1],$seqCounter);
	    push (@clusterTable,\@newCluster); 
	}

	# insert seqence into the seqsClusterTable (2. output parameter)
	my @seqClusterline=( @{$_}[0], @{$_}[1], $#clusterTable);
	push (@seqsClusterTable,\@seqClusterline);
	$seqCounter++;

    }
 
    # transform back the coordinates
    for my $i ( 0 .. $#clusterTable ) {
	for my $j ( 2 .. $#{ $clusterTable[$i] } ) {
	    $clusterTable[$i][$j]=$sortedIndices[$clusterTable[$i][$j]];
	}
    }  

    return (\@clusterTable,\@seqsClusterTable);

}

# check if two sequences determined by their start and end positions overlap
# there is no need for the regions to be sorted
sub doRegionsOverlap{
    my($r1Start,$r1End,$r2Start,$r2End)=@_;

    # in case that the indices are not in the correct direction, sort
    my $tempVar;
    if ($r1Start>$r1End){$tempVar=$r1Start;$r1Start=$r1End;$r1End=$tempVar;}
    if ($r2Start>$r2End){$tempVar=$r2Start;$r2Start=$r2End;$r2End=$tempVar;}

    # check if they overlap
 
    if (($r1End<$r2Start-2) or ($r1Start-2>$r2End))  {
	return 0;
    }else{
	return 1;
    }
}
# check if one region is completely contained in the other
sub isOneRegionContainedInTheOther{
    my($r1Start,$r1End,$r2Start,$r2End)=@_;

    # in case that the indices are not in the correct direction, sort
    my $tempVar;
    if ($r1Start>$r1End){$tempVar=$r1Start;$r1Start=$r1End;$r1End=$tempVar;}
    if ($r2Start>$r2End){$tempVar=$r2Start;$r2Start=$r2End;$r2End=$tempVar;}

    if ((($r1Start<=$r2Start) and ($r1End>=$r2End)) or (($r2Start<=$r1Start) and ($r2End>=$r1End)))  {
	return 1;
    }else{
	return 0;
    }
}



# takes two regions as an input and computes the spanning region
sub spanOfTwoRegions{
    my($r1Start,$r1End,$r2Start,$r2End)=@_;
     # in case that the indices are not in the correct direction, sort
    my $tempVar;
    if ($r1Start>$r1End){$tempVar=$r1Start;$r1Start=$r1End;$r1End=$tempVar;}
    if ($r2Start>$r2End){$tempVar=$r2Start;$r2Start=$r2End;$r2End=$tempVar;}

    my $spanStart;
    my $spanEnd;

    if ($r1Start<$r2Start){$spanStart=$r1Start;}else{$spanStart=$r2Start;}
    if ($r1End>$r2End){$spanEnd=$r1End;}else{$spanEnd=$r2End;}

    return ($spanStart,$spanEnd);
}

# generates the complement to a given set of sequences determined by its start and end points
# 3 input parameters
# 1. table with all the regions sorted and non-overlapping
# ([  1  60],
#  [ 65 120], 
#  [145 170],
#  [175 180])
#
# 2+3  start and end index which determine the window for which 
#      the complement is built (see drawing)
#      "where to start and where to stop"
#
#  regions     ---  ----     -----    ---    ---   
#  window      |start=1                              |end=39
#  complement     --    -----     ----   ----   ------
#
# 1 output parameter
# table with all the regions of the complement (start=1, end=213)
# ([ 61  64],
#  [121 144]
#  [171 174]
#  [181 213]
#
# usage:
# my @intergenicRegionsIndices=DOverlappingSequences::complementRegions(\@toCoReg,1,3267489);
#
sub complementRegions{
    my($tabRegionsRef)=@_;
    my @tabRegions=@{$tabRegionsRef};


    my $globalStart=$tabRegions[0][0];
    my $globalStop=$tabRegions[-1][1];
    
    my $tempStart=$globalStart;
    my $tempEnd;

    #DTableOps::printTable(\@tabRegions);

    my @tabComplement;
    foreach (@tabRegions){
	
	$tempEnd=${$_}[0]-1;
	# built a range from tempStart to tempEnd
	# and check if it is inside globalStart and globalStop 
	
	if ($tempEnd>$tempStart){
	    # the region has at least two elements, do further tests
	    
	    if (($tempStart>$globalStart-1) and ($tempEnd<$globalStop+1)) {
		# the region is completely within the window, so take it
		my @compLine=($tempStart,$tempEnd);
		push(@tabComplement,\@compLine); 
	    }
	    elsif( ($tempStart<$globalStart) and ($tempEnd>$globalStart) ){
		# region is chopped by the globalStart
		# find new boundaries and insert the region
		my @compLine=($globalStart,$tempEnd);
		push(@tabComplement,\@compLine);
	    }
	    elsif( ($tempStart<$globalStop) and ($tempEnd>$globalStop) ){
		#region is chopped by the globalStart
		# find new boundaries and insert the region
		my @compLine=($tempStart,$globalStop);
		push(@tabComplement,\@compLine);
	    }
		
	}

	# set a new tempStart for the next iteration
	$tempStart=${$_}[1]+1;
    }

    # if globalStop is beyond the last sequence, add another one   
    if ($tempStart<$globalStop){
	my @compLine=($tempStart,$globalStop);
	push(@tabComplement,\@compLine); 

    }


    return @tabComplement;

}

# converts the regions to real sequences
# 2 input parameters
# 1. table containing the start and end sites of the sequences
# ([ 3  7],
#  [ 9 13], 
#  [23 29],
#  [17 19])
#
# 2. a string containing e.g. the genome
# 
# genomoe: ATGTTTGTACCGCACGCCAAAAAGCCCGAAATTTA
# output:    GTTTG
#       :          ACCGC
#       :                        AGCCCGA 
#       :                  CCA                   
#
# note that the indices are converted (+1) since genome information 
# is usually (1..n) and not (0..n-1)
#
# 1 output parameter
# array of strings containing all the extracted seqs.

sub convertRegionsToSeqs{
    my($regionsTabRef,$sequence)=@_;
    my @regionsTab=@{$regionsTabRef};
    my @seqList;
    
    #DTableOps::printTable(\@regionsTab);

    foreach (@regionsTab){
	my $tempStart=${$_}[0];
	my $tempEnd=${$_}[1];
	
	my $locSequence=substr($sequence,$tempStart-1,$tempEnd-$tempStart+1);
	push(@seqList,$locSequence);

    }

    return @seqList;
}



1;



