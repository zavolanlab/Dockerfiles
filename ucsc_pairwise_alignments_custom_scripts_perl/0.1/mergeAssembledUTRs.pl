#!/usr/bin/env perl

BEGIN {
	use File::Basename;
	my $path_to_script = dirname(__FILE__);
	my $dtable_lib_path = "/app/lib/dPerlLibs/";
	unshift @INC, $dtable_lib_path;
}

use strict;
use DTable::DTableOps;
use DComputation::DPatternMatch;
use DComputation::DMath;
use DFileFormat::DFastaFileFormat;

# mutiplex all pairwise utr alignments in one file
# the aln files need to have a tag corresponding to the number in the final msa (aX)
# e.g. hg17_Reg-to-rn3_a6.aln


if(@ARGV < 1){
    die "USAGE: ./mergeAssembledUTRs.pl file1 file2 file3 ... \n";
}

# load the transcript fasta files
#my $trsTab=DFastaFileFormat::loadFileFast($ARGV[1]);
#my %trsHash=DFastaFileFormat::convertSeqsToHash($trsTab,'^\S+');

my @orgList;
my $cnt=1;
foreach (@ARGV){
    if ($_=~/\.aln$/){
	# check if the org name is encoded in the filename
	my $org;
	my $orgAlnNr;
	if ($_=~/\-to\-(\S+)\_a(\d+)\.aln/){
	    $org=$1;
	    $orgAlnNr=$2;
	}else{
	    $org="org";
	    $orgAlnNr=$cnt;$cnt++;
	}
	$orgList[$orgAlnNr-1]=[$_,$org,$orgAlnNr];
	#print "$org $orgAlnNr\n";
    }
}


# open AllFileHandles
my @fhhs;
foreach(@orgList){
    open(my $in, $_->[0] )   or die "Couldt read $_->[0]\n";
    push(@fhhs,$in);
}


# read all the files simultaneously

my $readFlag=1;
while  ($readFlag) {
    for (my $i=0;$i<=$#fhhs;$i++){
	my $fhl=$fhhs[$i];
	if (my $line=<$fhl>){
	    #$line=~/\>\>(\S+\|TR\(\d+\.\.\d+\)CDS\(\d+\.\.\d+\))/;
	    $line=~/\>\>(\S+)/;
	    my $id=$1;
	    my $seqQ=<$fhl>;
	    my $seqT=<$fhl>;
	    chomp($seqQ);
	    chomp($seqT);

	    #my $seqQsubst=substNsInAnchor($id,$seqQ,\%trsHash);

	    print ">$id $orgList[$i][1]\n";
	    print "$seqQ\n";
	    print "$seqT\n";

	}else{$readFlag=0;}
    }
    print "\n";
}

# close AllFileHandles
for (my$i=0;$i<=$#fhhs;$i++){
    my $ff=$fhhs[$i];
    close $ff;
}



sub substNsInAnchor{
    my ($id,$seqQ,$trsHash)=@_;


    if (defined $trsHash->{$id}){
	my $origTR=$trsHash->{$id};
	$id=~/\S+\|TR\(\d+\.\.(\d+)\)CDS\(\d+\.\.(\d+)\)/;
	my $trLen=$1;
	my $cdsTo=$2;

	# find all the stretches with N's
	my @matches=DPatternMatch::matchSeqToRegexp($seqQ,'N+');

	foreach(@matches){
	    my ($from,$to)=($_->[0],$_->[1]);

	    $seqQ=substr($seqQ,$from-1,$to-$from+1)=substr($origTR,$cdsTo+$from-1,$to-$from+1);
	    #$seqQ=substr($seqQ,$from-1,$to-$from+1);
	}


	DTableOps::printTable(\@matches);




    }else{die "$id is not known\n";}

}
