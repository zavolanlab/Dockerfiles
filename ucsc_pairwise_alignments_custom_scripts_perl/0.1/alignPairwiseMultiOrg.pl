#!/usr/bin/env perl

BEGIN {
	use File::Basename;
	my $path_to_script = dirname(__FILE__);
	my $dtable_lib_path = "/app/lib/dPerlLibs/";
	unshift @INC, $dtable_lib_path;
}

use strict;
use DTable::DTableOps;
use DFileFormat::DFastaFileFormat;


if(@ARGV != 3){
    die "USAGE: ./alignPairwiseMultiOrg.pl hg17_3UTR-to-MAM.mln hg_3UTR_GMAP.fa hg17\n";
}

my $anchorOrg=$ARGV[2];

# load the original transcipts
my $trs=DFastaFileFormat::loadFileFast($ARGV[1]);
my $trsHash=DFastaFileFormat::convertSeqsToHash($trs,'\S+');
#foreach(keys %{$trsHash}){print "$_\n";}

# read the file with the alignments
open(F,$ARGV[0]) or die "Cannot open $ARGV[0] \n";
while(<F>){
    if (/\>(\S+)\s+\S+/){
	my $trId=$1;
	my @mAln;
	while($_ =~ /\S/) {
	    /\>\S+\s+(\S+)/;
	    my $org=$1;
	    my $seqQ= <F>;
	    my $seqT= <F>;


	    chomp($seqQ);
	    chomp($seqT);

	    push(@mAln,[$seqQ,$seqT,$org]);
	    $_ = <F>;
	}

	my $aln=alignMultiUTR(\@mAln);

	my $originalAnchorSeq="";
	if (defined $trsHash->{$trId}){
	    $originalAnchorSeq=$trsHash->{$trId};
	}else{die "the sequence $trId is not in the $ARGV[1] \n";}


	my $replAlnAnchor=replaceAlnAnchorByOriginal($aln->[0],$trsHash->{$trId});

	print ">>$anchorOrg\_$trId\n";
	print "$replAlnAnchor\n";
	for (my $i=1;$i<@{$aln};$i++){
	    print ">$mAln[$i-1][2]\_$trId\n";
	    print "$aln->[$i]\n";
	}
	print "\\\\\n";
    }

}
close(F);



sub alignMultiUTR{
    my($mAln)=@_;

    #DTableOps::printTable($mAln);

    my $orgNr=$#{$mAln}+1;

    my @aln;
    # intialize position pointers
    my @posPointers;
    for (my $i=0;$i<$orgNr;$i++){push(@posPointers,0);}

    # traverse the assembled utrs and create alignment
    while($posPointers[0]<length($mAln->[0][0])){
	my @chars;
	for (my $i=0;$i<$orgNr;$i++){
	    push(@chars,substr($mAln->[$i][0],$posPointers[$i],1));
	}

	# find the chars that are "-"
	my %gapNrs;
	for (my $i=0;$i<@chars;$i++){
	    if ($chars[$i] eq "-"){
		$gapNrs{$i}=1;
	    }
	}

	#for (my $i=0;$i<$orgNr;$i++){
	#    print substr($mAln->[$i][0],$posPointers[$i],1)," ",$posPointers[$i]," ";
	#}
	#print "\n";


	# if there is not a single gap: save all
        if (keys(%gapNrs) == 0){
	    $aln[0].=substr($mAln->[0][0],$posPointers[0],1);
	    for (my $i=0;$i<$orgNr;$i++){
		$aln[$i+1].=substr($mAln->[$i][1],$posPointers[$i],1);
		$posPointers[$i]++;
	    }

	# deal with gaps
	}else{
	    if ($gapNrs{0}){
		$aln[0].=substr($mAln->[0][0],$posPointers[0],1);
	    }else{
		$aln[0].="-";
	    }

	    for (my $i=0;$i<$orgNr;$i++){
		if ($gapNrs{$i}){
		    $aln[$i+1].=substr($mAln->[$i][1],$posPointers[$i],1);
		    $posPointers[$i]++;
		}else{
		    $aln[$i+1].="-";
		}
	    }
	}


    }
    #foreach(@aln){print "$_\n";}
    #print "\n\n\n\n\n";

    return \@aln;
}


sub replaceAlnAnchorByOriginal{
    my($alnAnchor,$originalAnchor)=@_;

    my $replAlnAnchor=$alnAnchor;
    my $counter=-1;
    for (my $i=0;$i<length($alnAnchor);$i++){
	my $aA=substr($alnAnchor,$i,1);
	if ($aA ne "-"){
	    $counter++;

	    if(!defined substr($originalAnchor,$counter,1)){
		print STDERR $originalAnchor."\t".$counter."\t".length($originalAnchor)."\n";
            }
	    substr($replAlnAnchor,$i,1)=substr($originalAnchor,$counter,1);
	}

    }

    #print "$alnAnchor \n$replAlnAnchor\n\n";

    return $replAlnAnchor;
}
