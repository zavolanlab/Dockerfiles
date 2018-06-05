use strict;
use DTable::DTableOps;
use DAppWrapper::DRnaFold;
use DAppWrapper::DRnaHybrid;
use DComputation::DMath;
use DComputation::DRnaCut;
use DStat::StatsLite;
use DMol::DHybrid;
use DComputation::DBracket;
use DAppWrapper::DRnaEval;
package DOpenSecStructure;


sub openLocalFixedInv{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my @meanTransform=(0.2638,1.0864);
    my @sigmaTransform=(0.0794,1.5622);

    my $targetSiteLength=$hybrid->getTargetBindingSiteLength();
    my $mean=$meanTransform[0]*$targetSiteLength+$meanTransform[1];
    my $sigma=$sigmaTransform[0]*$targetSiteLength+$sigmaTransform[1];

    my $dE=openLocalFixed($target,$miRNA,$pos,$hybrid);
    #print "dE=$dE targetSiteLength=$targetSiteLength mean=$mean sigma=$sigma\n";
    my $p=StatsLite::normcdf($dE,$mean,$sigma);

    return ($dE,$p);
}


sub openLocalFixedHum{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my @meanTransform=(0.3373,1.4778);
    my @sigmaTransform=(0.1075,1.7999);

    my $targetSiteLength=$hybrid->getTargetBindingSiteLength();
    my $mean=$meanTransform[0]*$targetSiteLength+$meanTransform[1];
    my $sigma=$sigmaTransform[0]*$targetSiteLength+$sigmaTransform[1];

    my $dE=openLocalFixed($target,$miRNA,$pos,$hybrid);
    #print "$dE $mean $sigma\n";
    my $p=StatsLite::normcdf($dE,$mean,$sigma);

    return ($dE,$p);
}




sub openLocalFixed{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my $windowSize=100;

    my $targetSiteFrom=$pos;
    my $targetSiteTo=$hybrid->getTargetBindingSiteLength()+$pos-1;

    #print "$target,$pos,$hybrid\n";

    # extract local sequence around the binding Site
    my $limFrom=DMath::max($targetSiteFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteTo+int($windowSize/2),length($target));
    my $L1=substr($target,$limFrom-1,$targetSiteFrom-$limFrom);
    my $L3=substr($target,$targetSiteTo,$limTo-$targetSiteTo);  
    my $L2=nucChain("N",$hybrid->getTargetBindingSiteLength());
    my $unBoundLocalTarget=$L1.$L2.$L3;
    my $boundLocalTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);

    my $rnaFoldBound = new DRnaFold();
    $rnaFoldBound->fold($boundLocalTarget);

    my $rnaFoldUnbound = new DRnaFold();
    $rnaFoldUnbound->fold($unBoundLocalTarget);

    my $boundEE=$rnaFoldBound->getEnergyEnsemble();
    my $unboundEE=$rnaFoldUnbound->getEnergyEnsemble();

    my $dE=$unboundEE-$boundEE;

    #print "$targetSiteFrom $targetSiteTo $limFrom $limTo\n";
    #print "$unBoundLocalTarget\n";
    #print "$boundLocalTarget\n";
 
    my $sK=$boundEE-5*$dE;
   
    print $dE,"\n";

    return $dE;

}

sub openLocalVariable{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my $windowSize=200;

    my $targetSiteFrom=$pos;
    my $targetSiteTo=$hybrid->getTargetBindingSiteLength()+$pos-1;
    my $limFrom=DMath::max($targetSiteFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteTo+int($windowSize/2),length($target));
    my $seedLocalTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);
    my $seedTargetSiteFrom=$pos-$limFrom+1;
    my $seedTargetSiteTo=$hybrid->getTargetBindingSiteLength()+$seedTargetSiteFrom-1;


    # compute the PPM for the complete Window
    my $rnaFoldSeed = new DRnaFold();
    $rnaFoldSeed->fold($seedLocalTarget);
    my @PPM=$rnaFoldSeed->getPairingProbMatrix();

    # find the best place to cut the rna
    my $Lmin=1;
    my $Lmax=$targetSiteFrom-$limFrom;
    my $Rmin=$Lmax+$hybrid->getTargetBindingSiteLength()+1;
    my $Rmax=$Rmin+$limTo-$targetSiteTo-1;
    my ($bestL,$bestR)=DRnaCut::findOptimalCut(\@PPM,$Lmin,$Lmax,$Rmin,$Rmax);
    

    # extract local sequence around the binding Site
    my $L1=substr($seedLocalTarget,$bestL-1,$seedTargetSiteFrom-$bestL);
    my $L3=substr($seedLocalTarget,$seedTargetSiteTo,$bestR-$seedTargetSiteTo);  
    my $L2=nucChain("N",$hybrid->getTargetBindingSiteLength());
    my $unBoundLocalTarget=$L1.$L2.$L3;
    my $boundLocalTarget=substr($seedLocalTarget,$bestL-1,$bestR-$bestL+1);

    my $rnaFoldBound = new DRnaFold();
    $rnaFoldBound->fold($boundLocalTarget);

    my $rnaFoldUnbound = new DRnaFold();
    $rnaFoldUnbound->fold($unBoundLocalTarget);

    my $boundEE=$rnaFoldBound->getEnergyEnsemble();
    my $unboundEE=$rnaFoldUnbound->getEnergyEnsemble();

    my $dE=$unboundEE-$boundEE;

    #print "$seedLocalTarget\n";
    #print "$bestL $bestR\n";
    #print "$targetSiteFrom $targetSiteTo $limFrom $limTo\n";
    #print "$unBoundLocalTarget\n";
    #print "$boundLocalTarget\n";

    return $dE;
   

}

sub openSeedThenRest{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my $windowSize=100;
    my $seedSize=9;

    my $targetSiteGlobalFrom=$pos;
    my $targetSiteGlobalTo=$hybrid->getTargetBindingSiteLength()+$pos-1;

    # extract local sequence around the binding Site
    my $limFrom=DMath::max($targetSiteGlobalFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteGlobalTo+int($windowSize/2),length($target));
    my $localTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);
    my $targetSiteLocalFrom=$targetSiteGlobalFrom-$limFrom+1;
    my $targetSiteLocalTo=$targetSiteGlobalTo-$limFrom+1;
    my $targetSiteLocalSeedEnd=$targetSiteLocalTo-$seedSize+1;

    my $stOpen=substituteRegionWithN($localTarget,$targetSiteLocalFrom,$targetSiteLocalTo);
    my $stSeed=substituteRegionWithN($localTarget,$targetSiteLocalSeedEnd,$targetSiteLocalTo);

    # run rnafold
    my $rnaFoldStAll = new DRnaFold();
    $rnaFoldStAll->fold($localTarget);

    my $rnaFoldStOpen = new DRnaFold();
    $rnaFoldStOpen->fold($stOpen);

    my $rnaFoldStSeed = new DRnaFold();
    $rnaFoldStSeed->fold($stSeed);

    my $F_All=$rnaFoldStAll->getEnergyEnsemble();
    my $F_Open=$rnaFoldStOpen->getEnergyEnsemble();
    my $F_Seed=$rnaFoldStSeed->getEnergyEnsemble();

    # run rnahybrid
    my ($s1seed,$s2seed)=$hybrid->get5pHelix(1,9);
    my $rnaHybrid = new DRnaHybrid();
    $rnaHybrid->hybridize($s1seed,$s2seed);
    my $E_Seed=$rnaHybrid->getEnergy();

    #print "$s1seed\n$s2seed\n";

    my ($s1open,$s2open)=$hybrid->getTargetAndQuery();
    $rnaHybrid->hybridize($s1open,$s2open);
    my $E_Open=$rnaHybrid->getEnergy();


    #print "$p_seed $p_bound_seedBound\n";
    #print "F_All=$F_All  F_Open=$F_Open E_Open=$E_Open F_Seed=$F_Seed E_Seed=$E_Seed\n";
    #print "$s1open $s2open $E_Open\n";
    #print "$target\n $pos\n";
    #print "$localTarget\n$stOpen\n$stSeed\n";

    return "F_All=$F_All;F_Open=$F_Open;E_Open=$E_Open;F_Seed=$F_Seed;E_Seed=$E_Seed;";

}



sub open5pEnd{
    my($target,$miRNA,$pos,$hybrid)=@_;

    my $windowSize=200;
    my $seed=8;

    my $targetSiteTo=$hybrid->getTargetBindingSiteLength()+$pos-1;
    my $targetSiteFrom=$targetSiteTo-$seed+1;

    # extract local sequence around the binding Site
    my $limFrom=DMath::max($targetSiteFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteTo+int($windowSize/2),length($target));
    my $localTargetSiteFrom=$targetSiteFrom-$limFrom+1;
    my $localTargetSiteTo=$seed+$localTargetSiteFrom-1;

    my $boundLocalTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);

    my $rnaFoldBound = new DRnaFold();
    $rnaFoldBound->fold($boundLocalTarget);
    my @puv=$rnaFoldBound->getUnboundProbVector();
    my @puvLocal=@puv[$localTargetSiteFrom-1..$localTargetSiteTo-1];

    my ($seedAln2,$seedAln3)=$hybrid->get5pSeed($seed);

    my %hybridEnergies;
    $hybridEnergies{"GC"}=-3;
    $hybridEnergies{"CG"}=-3;
    $hybridEnergies{"AU"}=-2;
    $hybridEnergies{"UA"}=-2;
    $hybridEnergies{"GU"}=-1;
    $hybridEnergies{"UG"}=-1; 
    $hybridEnergies{"  "}=0;
    

    my $hybScore=0;
    for my $i(0..length($seedAln2)-1){
	my $c1=substr($seedAln2,$i,1);
	my $c2=substr($seedAln3,$i,1);
	my $dE=$hybridEnergies{"$c1$c2"}*$puvLocal[$i];
	$hybScore+=$dE;
    }

    #print "$seedAln2\n$seedAln3\n";
    #print "$freeScore\n";
    #print "@puvLocal\n";
    #print "$target\n$boundLocalTarget\n$localTargetSiteFrom\n$localTargetSiteTo\n";

    return ($hybScore,1);

}


sub openPartOfSeed{
    my($target,$miRNA,$pos,$hybrid,$seedStart,$seedLength)=@_;

    my $windowSize=200;

    my $targetSiteFrom=$pos+$hybrid->getTargetBindingSiteLength()-$seedLength-$seedStart+1;
    my $targetBindingSiteLength=$seedLength;
    my $targetSiteTo=$targetBindingSiteLength+$targetSiteFrom-1;


    #print "$target,$pos,$hybrid\n";

    # extract local sequence around the binding Site
    my $limFrom=DMath::max($targetSiteFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteTo+int($windowSize/2),length($target));
    my $L1=substr($target,$limFrom-1,$targetSiteFrom-$limFrom);
    my $L3=substr($target,$targetSiteTo,$limTo-$targetSiteTo);  
    my $L2=nucChain("N",$targetBindingSiteLength);
    my $unBoundLocalTarget=$L1.$L2.$L3;
    my $boundLocalTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);

    my $rnaFoldBound = new DRnaFold();
    $rnaFoldBound->fold($boundLocalTarget);

    my $rnaFoldUnbound = new DRnaFold();
    $rnaFoldUnbound->fold($unBoundLocalTarget);

    my $boundEE=$rnaFoldBound->getEnergyEnsemble();
    my $unboundEE=$rnaFoldUnbound->getEnergyEnsemble();

    my $dE=$unboundEE-$boundEE;

    #print "$targetSiteFrom $targetSiteTo $limFrom $limTo\n";
    #print "$unBoundLocalTarget\n";
    #print "$boundLocalTarget\n";

    return ($dE,1);

}

sub dGtoMakeAccessibleUnzipTarget{
    my($target,$miRNA,$pos,$hybrid)=@_;
 

    my $windowSize=100;
    my $seedSize=4;

    my $targetSiteGlobalFrom=$pos;
    my $targetSiteGlobalTo=$hybrid->getTargetBindingSiteLength()+$pos-1;

    # extract local sequence around the binding Site
    my $limFrom=DMath::max($targetSiteGlobalFrom-int($windowSize/2),1);
    my $limTo=DMath::min($targetSiteGlobalTo+int($windowSize/2),length($target));
    my $localTarget=substr($target,$limFrom-1,$limTo-$limFrom+1);
    my $targetSiteLocalFrom=$targetSiteGlobalFrom-$limFrom+1+1;
    my $targetSiteLocalTo=$targetSiteGlobalTo-$limFrom+1-1;
    my $targetSiteLocalSeedEnd=$targetSiteLocalTo-$seedSize+1;

    my $targetSiteBindingLength=$hybrid->getTargetBindingSiteLength();

    # fold the complete target
    my ($targetEnergy,$targetBracketNot)=DRnaFold::invokeRNAFoldFast($localTarget);

    my $modBracket=removeBracketPairsInRegion($targetBracketNot,$targetSiteLocalSeedEnd-1,$targetSiteLocalTo-1);
    my $unzippedEnergy=DRnaEval::evaluateSeqWithBracket($localTarget,$modBracket);

    my $dE=$unzippedEnergy-$targetEnergy; 
    print "$dE\n";

    return ($dE,1);

}

sub removeBracketPairsInRegion{
    my ($bracket,$from,$to)=@_;  # bis und mit!!

    my $outBracket;
    my @bracketInList = split(//,$bracket);

    my %pairs=DBracket::getPairHash($bracket);

    for my $i($from..$to){
	# check if ind points to an empty place or to a bracket
	if (substr($bracket,$i,1) ne "."){
	    my $partner=$pairs{$i};
	    $bracketInList[$i]=".";
	    $bracketInList[$partner]=".";
	}
    }


    $outBracket=join("",@bracketInList);

    return $outBracket;

}


sub substituteRegionWithN{
    my ($s,$from,$to)=@_;

    my $L1=substr($s,0,$from-1);
    my $L3=substr($s,$to,length($s)-$to);  
    my $L2=nucChain("N",$to-$from+1);
    my $sOut=$L1.$L2.$L3;

    return $sOut;
}


sub nucChain{
    my ($nuc,$nrOfNucs)=@_;
	
    my $nucsLine=""; 
    for my $i(0..$nrOfNucs-1){
	$nucsLine.=$nuc;
    }
    return $nucsLine;
}

1;
