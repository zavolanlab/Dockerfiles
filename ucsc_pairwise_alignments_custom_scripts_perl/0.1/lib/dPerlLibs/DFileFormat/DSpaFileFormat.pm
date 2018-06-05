use strict;
use DTable::DTableOps;
use DComputation::DMath;
use DString::DStringOps;
package DSpaFileFormat;

sub new {
    my $class = shift;
    my $self = {};
    bless $self;
    my $filename = shift;

    my %coordsHash;
    open(F,$filename) or die "Cannot open $filename \n";

    while(<F>){
	next if $_ =~ /^\#/;
	my @lines;
	while(!($_ =~ /\/\//)) {
	    push @lines, $_;
	    $_ = <F>;
	}
	my @aln=DTableOps::convertLinesToTable(\@lines,'\t');

	my $mrna=$aln[0][0];
	my $chr=$aln[0][3];

	if (defined $coordsHash{$mrna}){die "error in DSpaFileFormat. $mrna is two times in $filename";}

	$coordsHash{$mrna}=[$chr,\@aln];

    }

    $self->{coordsHash}=\%coordsHash;
    return $self;
}

sub new22 {
    my $class = shift;
    my $self = {};
    bless $self;
    my $filename = shift;

    my %coordsHash;
    open(FHANDLE,$filename) or die "Cannot open $filename \n";

    while(<FHANDLE>){
	if ($_=~/seq1 \= (\S+)\, /){
	    my $mrna=$1;
	    my $genLocus=<FHANDLE>;
	    $genLocus=~/seq2 \= (\S+)\, /;
	    my $chr=$1;

	    $_=<FHANDLE>;
	    $_=<FHANDLE>;
	    $_=<FHANDLE>;

	    # extract the coords
	    my @ll;
	    while($_ =~ /\S/) {
		push @ll, $_;
		$_ = <FHANDLE>;
	    }
	    # convert to table
	    my @genLoci=DTableOps::convertLinesToTable(\@ll,'\t');

	    # extract the alignment
	    




	    # save into a hash of array of array
	    if (defined $coordsHash{$mrna}){die "error in DSpaFileFormat. $mrna is two times in $filename";}

	    $coordsHash{$mrna}=[$chr,\@genLoci];
	    
	    #print "$mrna $chr \n";
	    #DTableOps::printTable(\@genLoci);
	}
    }

    $self->{coordsHash}=\%coordsHash;

    return $self;
}

sub getGeneticLocus{
    my $self = shift;
    my $mrna = shift;	
    my $from = shift;
    my $to   = shift;

    if (!defined ${$self->{coordsHash}}{$mrna}){return ("",0,0)};

    my $chr=${$self->{coordsHash}}{$mrna}[0];
    my $genLociRef=${$self->{coordsHash}}{$mrna}[1];

    my $locusStart;
    my $locusEnd;
    my $locusStrand;

    #DTableOps::printTable($genLociRef);

    # find the genetic locus in the correct exon 
    for my $i(1..$#{$genLociRef}){
	my ($exonStart,$exonEnd,$strand)=($genLociRef->[$i][2],$genLociRef->[$i][3],$genLociRef->[$i][8]);
	$locusStrand=$strand;
 
	my $overlap=doRegionsOverlap($exonStart,$exonEnd,$from,$to);
	if ($overlap==1){
	    # cut away region which is outside of the exon
	    my $fromN=DMath::max($exonStart,$from)-$exonStart+1;
	    my $toN=DMath::min($exonEnd,$to)-$exonStart+1;

	    my $cloneSeq=$genLociRef->[$i][4];
	    my $exonSeq=$genLociRef->[$i][7];
	    
	    my $exonGenLocStart=$genLociRef->[$i][5];
	    my $exonGenLocEnd=$genLociRef->[$i][6];	

	    if ($strand eq "+"){
		my $aA=DStringOps::goToPosMatch($cloneSeq,$fromN);
		my $aE=DStringOps::goToPosMatch($cloneSeq,$toN);
		my $bA=DStringOps::goToPosMatchRev($exonSeq,$aA);
		my $bE=DStringOps::goToPosMatchRev($exonSeq,$aE);

		$locusStart=$exonGenLocStart+$bA-1;
		$locusEnd=$exonGenLocStart+$bE-1;
	    }else{
		#$cloneSeq=reverse $cloneSeq;
		#$exonSeq=reverse $exonSeq;
		
		my $aA=DStringOps::goToPosMatch($cloneSeq,$fromN);
		my $aE=DStringOps::goToPosMatch($cloneSeq,$toN);
		my $bA=DStringOps::goToPosMatchRev($exonSeq,$aA);
		my $bE=DStringOps::goToPosMatchRev($exonSeq,$aE);

		$locusEnd=$exonGenLocStart-$bA+1;
		$locusStart=$exonGenLocStart-$bE+1;


	    }


	    my $a=DStringOps::goToPosMatchRev("AGTGCT--GTA---ATG",9);

	    #print "$exonStart $exonEnd $strand\n";
	    #print "$fromN $toN\n";
	    #print "$exonGenLocStart $exonGenLocEnd\n";
	    #print "$cloneSeq\n$exonSeq\n";
	}
    }
    

    #print "$chr\n";
    #DTableOps::printTable($genLociRef);
    #print "$locusStart $locusEnd\n";

    return ($chr,$locusStart,$locusEnd,$locusStrand);
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
 
    if (($r1End<$r2Start) or ($r1Start>$r2End))  {
	return 0;
    }else{
	return 1;
    }
}

1;
