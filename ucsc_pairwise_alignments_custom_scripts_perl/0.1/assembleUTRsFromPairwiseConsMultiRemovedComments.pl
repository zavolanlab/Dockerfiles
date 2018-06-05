#!/usr/bin/env perl

BEGIN {
	use File::Basename;
	my $path_to_script = dirname(__FILE__);
	my $dtable_lib_path = "/app/lib/dPerlLibs/";
	unshift @INC, $dtable_lib_path;
}

use strict;
use DTable::DTableOps;
use DComputation::DMath;

if(@ARGV != 5){
    die "USAGE: ./assembleUTRsFromPairwiseConsMulti.pl dir hg17 \"panTro1,rheMac1,canFam1,bosTau1,mm6,rn3,monDom1,galGal2\" hg_TR-to-hg17.match.3UTR.tab .\n";
}

my $UNK="N";

my ($dir,$anchorOrg,$alnOrgsSt,$regFile,$outDir)=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
my @alnOrgs=split('\,',$alnOrgsSt);

my $counter=1;
for my $alnOrg(@alnOrgs){
    my $alnDir="$dir/$anchorOrg\_to\_$alnOrg";
    my $outFile= "$outDir/$anchorOrg\_Reg-to-$alnOrg\_a$counter.aln";
    print "...generating $outFile\n";
    assembleAlignment($alnDir,$regFile,$UNK,$outFile);
    $counter++;
}


sub assembleAlignment{
    my($alnDir,$regFile,$UNK,$outFile)=@_;

    open(OUT, ">".$outFile) or die "cannot write to output file $outFile\n";

    #read the regions file
    my $regLines=DTableOps::loadFileRef($regFile,'');
    my $regions=DTableOps::convertLinesToTableRef($regLines,'\s+');

    # split the file according to chromosome
    # and create hash: id->strand
    my %regionsHash;
    my %trStrandHash;
    foreach(@{$regions}){
	my $key="$_->[1]";

	if (defined $regionsHash{$key}){
	    push(@{$regionsHash{$key}},$_);
	}else{
	    $regionsHash{$key}=[$_];
	}

	$_->[0]=~/(\S+)\.\d+$/;
	my $trIDlong=$1;

	my $trIDshort=$trIDlong;
	if ($_->[0]=~/(\S+)\|TR/){
	    $trIDshort=$1;
	}

	my $strand=$_->[4];

	$trStrandHash{$trIDshort}=[$trIDlong,$strand];
    }


    # extract the alignments
    foreach(sort keys %regionsHash){ # for each chromosome
	#print "$_\n";
	my ($coordToTR,$minCoord,$maxCoord)=getCoordToTR($regionsHash{$_});

	# get all the regions for the actual chromosome and sort
	my @regionsSummary=DTableOps::selectColumns($regionsHash{$_},[2,3]);
	my ($regionsSummarySorted,$regionsSummarySortedInd)=DTableOps::sortTableForColumn(\@regionsSummary,0);

	#DTableOps::printTable($regionsSummarySorted);

	# create hash with the alignments which will grow during the time
	my %trAlnHash;
	foreach(keys %{$coordToTR}){$trAlnHash{$coordToTR->{$_}[0]}=["",""];}

	my $lastCounter=$minCoord-1;
	my $chrFile="$alnDir/$_.axt";

	if(!(open(F,$chrFile))){
	    print STDERR "The alignment file $chrFile does not exist\n";
	}else{

	    while(<F>){
        	next if $_ =~ /^\#/;
		if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
		    my ($nr,$chrQ,$fromQ,$toQ,$chrT,$fromT,$toT,$strandT,$score)=($1,$2,$3,$4,$5,$6,$7,$8,$9);
		    my $seqQ=<F>;
		    $seqQ=~s/\n//g;
		    my $seqT=<F>;
		    $seqT=~s/\n//g;

		    # remove all regions which have been passed in respect to the intermediate piece
		    my $cou=0;
		    while(($cou<=$#{$regionsSummarySorted}) and ($regionsSummarySorted->[$cou][1]<$lastCounter)){
			shift(@{$regionsSummarySorted});
		    }

		    # check if the intermediate piece overlaps with anything in our set $regionsSummarySorted
		    my $overlapINT=0;
		    $cou=0;
		    while(($cou<=$#{$regionsSummarySorted}) and ($regionsSummarySorted->[$cou][0]<=$fromQ-1)){
			if ( ($regionsSummarySorted->[$cou][1]>=$lastCounter+1) and ($regionsSummarySorted->[$cou][0]<=$fromQ-1)){
			    $overlapINT=1;
			}
			$cou++;
		    }

		    # remove all regions which have been passed in respect to the aligned piece
		    $cou=0;
		    while(($cou<=$#{$regionsSummarySorted}) and ($regionsSummarySorted->[$cou][1]<$fromQ)){
			shift(@{$regionsSummarySorted});
		    }

		    # check if this piece overlaps with anything in our set $regionsSummarySorted
		    $cou=0;
		    my $overlapALN=0;
		    while(($cou<=$#{$regionsSummarySorted}) and ($regionsSummarySorted->[$cou][0]<=$toQ)){
			if ( ($regionsSummarySorted->[$cou][1]>=$fromQ) and ($regionsSummarySorted->[$cou][0]<=$toQ)){
			    $overlapALN=1;
			}
			$cou++;
		    }

		    # check if some positions fall between the last block($lastCounter) and the new one.
		    if ($overlapINT){
			for (my $j=$lastCounter+1;$j<=$fromQ-1;$j++){
			    if (defined $coordToTR->{$j}){
				growALN($coordToTR,\%trAlnHash,$j,$UNK,$UNK);
			    }
			}

		    }
		    if ($overlapALN){
			#traverse the alignment
			my $tempCounter=$fromQ-1;

			my $cQ;
			my $cT;
			my $seqQlen=length($seqQ);
			for (my $i=0;$i<$seqQlen;$i++){
			    $cQ=substr($seqQ,$i,1);
			    $cT=substr($seqT,$i,1);
			    if ($cQ ne "-"){$tempCounter++};
			    if ((defined $coordToTR->{$tempCounter}) and ($tempCounter>$lastCounter)){ # 2. term for the case where alns overlap (in axt)
				growALN($coordToTR,\%trAlnHash,$tempCounter,$cQ,$cT);
			    }
			}
		    }

		    $lastCounter=DMath::max($toQ,$lastCounter); # $lastCounter in the case that the last block was even longer than this one
		    #print "$nr\n";
		}
	    }

	    close(F);
	}

	# finally check if there are still regions after the last alignment block

	if ($#{$regionsSummarySorted}>=0){
	    # find maximum in remaining
	    my $maxInRem=0;
	    for my $z(@{$regionsSummarySorted}){$maxInRem=DMath::max($maxInRem,$z->[1]);}

	    #for (my $j=$lastCounter+1;$j<=$regionsSummarySorted->[$#{$regionsSummarySorted}][1];$j++){
	    for (my $j=$lastCounter+1;$j<=$maxInRem;$j++){
		if (defined $coordToTR->{$j}){
		    growALN($coordToTR,\%trAlnHash,$j,$UNK,$UNK);
		    #DTableOps::printTable($regionsSummarySorted);
		}
	    }
	}


	# print the resulting utrs
	foreach(sort keys %trAlnHash){

	    my$trId=$trStrandHash{$_}->[0];
	    my$strand=$trStrandHash{$_}->[1];
	    my$seqQ=$trAlnHash{$_}->[0];
	    my$seqT=$trAlnHash{$_}->[1];

	    # reverse complement
	    if ($strand eq "-"){
		$seqQ =~ tr/[ACGTNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCANYRWSKMVDHBtgcanyrwskmvdhb]/;
		$seqT =~ tr/[ACGTNRYWSMKBDHVacgtnrywsmkbdhv]/[TGCANYRWSKMVDHBtgcanyrwskmvdhb]/;
		$seqQ = reverse $seqQ;
		$seqT = reverse $seqT;
	    }

	    print OUT ">>$trId\n";
	    print OUT $seqQ,"\n";
	    print OUT $seqT,"\n";
	}

    }
    close(OUT);
}

sub growALN{
    my($coordToTR,$trAlnHash,$coord,$alnQ,$alnT)=@_;

    # get the mRNAs which fall into the $coord
    # grow the alignment by $alnQ and $alnT
    foreach(@{$coordToTR->{$coord}}){
	$trAlnHash->{$_}[0].=$alnQ;
	$trAlnHash->{$_}[1].=$alnT;
    }
    #print "$coord\n";

}



sub getCoordToTR{
    my($chrRegTabRef)=@_;

    my @chrRegTab=@{$chrRegTabRef};
    my %trRegions; # contains all the regions for one transcript
    foreach(@chrRegTab){

	$_->[0]=~/(\S+)\.\d+$/;
	my $tID=$1;
	if ($_->[0]=~/(\S+)\|TR/){$tID=$1;}

	if (defined $trRegions{$tID}){
	    push(@{$trRegions{$tID}},$_);
	}else{
	    $trRegions{$tID}=[$_];
	}
    }

    # enumarate all the coordinates
    my %coordToTR;
    my $minCoord=1000000000;
    my $maxCoord=1;
    foreach(keys %trRegions){
	my $idTR=$_;
	my @tab=@{$trRegions{$_}};
	my @trEnumCoords;
	for (my $i=0;$i<=$#tab;$i++){
	    my $from=$tab[$i][2];
	    my $to=$tab[$i][3];
	    my @exonEnum=($from..$to);
	    @trEnumCoords=(@trEnumCoords,@exonEnum);
	    # keep the min and the max for later
	    $minCoord=DMath::min($minCoord,$from);
	    $maxCoord=DMath::max($maxCoord,$to);
	}
	#print "@trEnumCoords\n";
	#DTableOps::printTable($trRegions{$_});

	foreach(@trEnumCoords){
	    if (defined $coordToTR{$_}){
		push(@{$coordToTR{$_}},$idTR);
	    }else{
		$coordToTR{$_}=[$idTR];
	    }
	}
    }

    return (\%coordToTR,$minCoord,$maxCoord);

}
