use strict;
use DTable::DTableOps;
use DConst::DBioConst;
use DComputation::DMath;
package PositionScoreAlignment;


# ----------------------------------------------------------------------
# Global constants

my $matrix=\@DBioConst::rnaSimpleHybrid;

# The traceback matrices are indexed by (direction, row, column).

my @DIR = (1, 2, 3);
my $STOP = 0;

# Directions in the linear (2D) traceback: 
# 0=stop; 1=from North (above); 2=from Northwest; 3=from West (left)
my ($FROMN, $FROMNW, $FROMW) = @DIR;

# Directions in the affine (3D) traceback: 
my ($FROMM, $FROMIX, $FROMIY) = @DIR;

# Minus infinity

my $minusInf = -2111111111;     # Representable in 32 bits 2's compl.

# Color codes for the traceback
my ($RED, $BLUE, $GREEN) = (1, 2, 3);





sub dimosAlign{
    my ($x,$posOffsetx,$y,$posOffsety)=@_;      # By val, ref, val, ref, ref, val

    #check if the inputs have the correct dimensionality
    if ( (length($x)!= $#{$posOffsetx}+1) or (length($y)!= $#{$posOffsety}+1) ){
       die "Position Score Alignment error, wrong dimensionality";
    }
    # convert DNA to RNA and convert to uppercase
    $x=formatSeq($x);
    $y=formatSeq($y);
  
    # reverse the SRNA so it can bind to the target
    $x= reverse $x;

    my ($n, $m) = (length($x), length($y));


    #initialize dp and backtrack matrices
    my @F;
    my @Fback;

    for (my $j=0; $j<=$m; $j++) {
	 $F[0][$j] = [(0) x ($n+1)];
	 $F[1][$j] = [(0) x ($n+1)];
	 $Fback[0][$j] = [(0) x ($n+1)];
	 $Fback[1][$j] = [(0) x ($n+1)];
    }

    # direction codes:
    #
    # low       | 1| 3|  high       | 2| 4|
    # (match)   | 5|xl|  (no match) | 6|xh|

    my $mOpenP=6;  # missmatch opening penalty
    my $mExtP=0.2; # missmatch extend penalty
    my $gOpenP=4;  # gap opening penalty
    my $gExtP=0.4; # gap extend penalty


    my ($Nl,$NWl,$Wl,$Nh,$NWh,$Wh);

    # start dp
    for (my $i=1; $i<=$n; $i++) {
	for (my $j=1; $j<=$m; $j++) {

	     # form the scoring matrix
	     my $s = DBioConst::nucScore($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));

	     # position score
	     #my $pS= $s+${$posOffsetx}[$i-1]+ ${$posOffsety}[$j-1];
	     my $pS= $s*${$posOffsetx}[$i-1]* ${$posOffsety}[$j-1];
	     if ($s<0){$pS=-9};
	     
	     $Nl =$F[0][$j-1][$i]; $NWl=$F[0][$j-1][$i-1]; $Wl =$F[0][$j][$i-1];
	     $Nh =$F[1][$j-1][$i]; $NWh=$F[1][$j-1][$i-1]; $Wh =$F[1][$j][$i-1];

	     # the maximum returns the winning index. 0 for new start
	     my ($xl,$xlInd)=DMath::maximum(0,            #(0)
			                    $NWl+$pS,     #(1)
			                    $NWh+$pS);    #(2) 
				   
	     my ($xh,$xhInd)=DMath::maximum(0,            #(0)
				            $NWl-$mOpenP, #(1)
				            $NWh-$mExtP,  #(2)
					    $Nl-$gOpenP,  #(3)
					    $Nh-$gExtP,   #(4)
					    $Wl-$gOpenP,  #(5)
					    $Wh-$gExtP);  #(6)

	     # write the max scores into the dp matrices
	     $F[0][$j][$i] = $xl;
	     $F[1][$j][$i] = $xh;

	     #save the two winners for low and high
	     $Fback[0][$j][$i]=$xlInd;
	     $Fback[1][$j][$i]=$xhInd;
	}
    }

    # Find maximal score in matrix
    my $sMax = $F[0][0][0];
    my ($sMaxI,$sMaxJ,$sMaxH);  # vMaxH: 0="low",  1="high"
    for (my $i=1; $i<=$n; $i++) {
	for (my $j=1; $j<=$m; $j++) {
	    my ($tMax,$tMaxInd)=DMath::maximum($F[0][$j][$i],$F[1][$j][$i],$sMax);
	    if ($tMaxInd != 2){
		$sMax=$tMax;
		$sMaxI=$i; 
		$sMaxJ=$j;
		$sMaxH=$tMaxInd;
	    }
	}
    }  
    # the winner must be in the low matrix since it is a local alignment
    # transition from low to high gives penalty
    if ($sMaxH==1){ die "error in Position Score Alignment: winner is in the high-matrix\n" }
    

    # backtrack to find solution
    my ($btScore,$btI,$btJ,$btH)=($sMax,$sMaxI,$sMaxJ,$sMaxH);

    # build the alignment
    # seqA  AGUAGGUA--GAAA
    # seqM  |||   ||  ||||
    # seqB  UCAAG-AUAUCUUU
    my($seqA,$seqB,$seqM)=("","","");

    while ($btScore != 0){

	# direction codes:
	#
	# low       | 1| 3|  high       | 2| 4|
	# (match)   | 5|xl|  (no match) | 6|xh|
	# 0 : new start, stop backtracking

	my $dir=$Fback[$btH][$btJ][$btI];

	if ($btH==0){

	    if   ($dir==0){$btScore=0;
	    }elsif ($dir==1){
		$seqA.=substr($x,$btI-1,1);
		$seqB.=substr($y,$btJ-1,1);
		$seqM.="|";

		$btI--;$btJ--;$btH=0;
	
            }elsif ($dir==2){
		$seqA.=substr($x,$btI-1,1);
		$seqB.=substr($y,$btJ-1,1);
		$seqM.="|";

		$btI--;$btJ--;$btH=1;
	    }

	}else{

	    if     ($dir==0){$btScore=0;
	    }elsif ($dir==1){
		$seqA.=substr($x,$btI-1,1);
		$seqB.=substr($y,$btJ-1,1);
	        $seqM.=" ";

		$btI--;$btJ--;$btH=0;

	    }elsif ($dir==2){
		$seqA.=substr($x,$btI-1,1);
		$seqB.=substr($y,$btJ-1,1);
	        $seqM.=" ";

		$btI--;$btJ--;$btH=1;
		
	    }elsif ($dir==3){
		$seqA.="-";
		$seqB.=substr($y,$btJ-1,1);
	        $seqM.=" ";

		$btJ--;$btH=0;

	    }elsif ($dir==4){
		$seqA.="-";
		$seqB.=substr($y,$btJ-1,1);
	        $seqM.=" ";

		$btJ--;$btH=1;
	    }elsif ($dir==5){
		$seqA.=substr($x,$btI-1,1);
		$seqB.="-";
	        $seqM.=" ";

		$btI--;$btH=0;

	    }elsif ($dir==6){
		$seqA.=substr($x,$btI-1,1);
		$seqB.="-";
	        $seqM.=" ";

		$btI--;$btH=1;
	    }
       }
	    
    }
    $seqA=reverse $seqA;
    $seqB=reverse $seqB;
    $seqM=reverse $seqM;
    my $alignment="   $seqA\n   $seqM\n   $seqB\n";

    #DTableOps::printMatrix($F[0]);
    #DTableOps::printMatrix($Fback[1]);

    #print "$sMax $sMaxI,$sMaxJ \n";
    #print "$alignment\n";
    
    my $position=$btJ+1;

    return (-$sMax,$position,$alignment);
}







# ----------------------------------------------------------------------
# The Smith-Waterman local alignment algorithm, linear gap costs
# Input: The sequences $x @positionOffsetx and $y @postitionOffsety to align
# Output: references to the F and B matrices, and the aligned sequences

# This is a special alignment. It takes not only the 2 sequences as input, 
# but it also takes two vectors of floats. Each one of them is associated with
# a seqence 
# The weights of these vectors are added nucleotidewise to the scoring matrix
# during the alignment.

sub localAlignLinear {
  my ($x,$posOffsetx,$y,$posOffsety,$e)=@_;      # By val, ref, val, ref, ref, val

  #check if the inputs have the correct dimensionality
  if ( (length($x)!= $#{$posOffsetx}+1) or (length($y)!= $#{$posOffsety}+1) ){
      die "Position Score Alignment error, wrong dimensionality";
  }
  # convert DNA to RNA and convert to uppercase
  $x=formatSeq($x);
  $y=formatSeq($y);
  
  # reverse the SRNA so it can bind to the target
  $x= reverse $x;


  my ($n, $m) = (length($x), length($y));
  # The dynamic programming matrix; also correctly initializes borders
  my @F; 
  for (my $j=0; $j<=$m; $j++) {
    $F[$j] = [(0) x ($n+1)];
  }
  

  # The traceback matrix; also correctly initializes borders
  my @B; 
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      $B[$dir][$j] = [($STOP) x ($n+1)];
    }
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s = DBioConst::nucScore($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      my $val = DBioConst::max(0, 
                     $F[$j-1][$i-1]+$s+ ${$posOffsetx}[$i-1]+ ${$posOffsety}[$j-1] , 
                     $F[$j][$i-1]-$e, 
                     $F[$j-1][$i]-$e);
      $F[$j][$i] = $val;
      # Record all traceback directions (none if we restart at score 0):
      if ($val == $F[$j-1][$i-1]+$s+ ${$posOffsetx}[$i-1]+ ${$posOffsety}[$j-1]) {
        $B[$FROMNW][$j][$i] = $RED;
      } 
      if ($val == $F[$j][$i-1]-$e) {
        $B[$FROMW][$j][$i] = $RED;
      } 
      if ($val == $F[$j-1][$i]-$e) {
        $B[$FROMN][$j][$i] = $RED;
      } 
    }
  }
  # Find maximal score in matrix
  my $vmax = 0;
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      $vmax = DBioConst::max($vmax, $F[$j][$i]);
    }
  }  
  my ($jmax, $imax) = (0, 0);
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      if ($vmax == $F[$j][$i]) {
        &markReachable2(\@B, $i, $j);
        $jmax = $j;
        $imax = $i;
      }
    }
  }

  #DTableOps::printTable($B[1]);
  #DTableOps::printTable($B[2]);
  #DTableOps::printTable($B[3]);

  #print $B[1][0][0];

  return (\@F, $vmax, &traceback2($x, $y, \@B, $imax, $jmax));
}


# ----------------------------------------------------------------------
# Common subroutines for linear gap cost routines

# Reconstruct the alignment from the traceback, backwards, from ($i, $j)

sub traceback2 {
  my ($x, $y, $B, $i, $j) = @_;         # B by reference
  my ($xAlign, $yAlign) = ("", "");
  while ($$B[$FROMN][$j][$i] || $$B[$FROMW][$j][$i] || $$B[$FROMNW][$j][$i]) {
    if ($$B[$FROMN][$j][$i]) {
      $$B[$FROMN][$j][$i] = $GREEN;
      $xAlign .= "-"; 
      $yAlign .= substr($y, $j-1, 1);
      $j--;
    } elsif ($$B[$FROMW][$j][$i]) {
      $$B[$FROMW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-"; 
      $i--;
    } elsif ($$B[$FROMNW][$j][$i]) {
      $$B[$FROMNW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $i--; $j--;
    }
  }
  # Warning: these expressions cannot be inlined in the list
  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign, $yAlign);
}

# Mark all traceback arrows reachable from a ($i, $j)

sub markReachable2 {
  my ($B, $i, $j) = @_;         # B by reference
  if ($$B[$FROMN][$j][$i] == $RED) {
    $$B[$FROMN][$j][$i] = $BLUE;
    &markReachable2($B, $i, $j-1);
  } 
  if ($$B[$FROMW][$j][$i] == $RED) {
    $$B[$FROMW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j);
  } 
  if ($$B[$FROMNW][$j][$i] == $RED) {
    $$B[$FROMNW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j-1);
  }
}




# ----------------------------------------------------------------------
# The Smith-Waterman local alignment algorithm, affine gap costs
# Input: The sequences $x and $y to align
# Output: references to the matrices M, Ix, Iy, B, and the aligned sequences

sub localAlignAffine {
  my ($x,$posOffsetx,$y,$posOffsety,$d,$e) = @_;
  my ($n, $m) = (length($x), length($y));
  # Initialize upper and left-hand borders
  # M represent an aa/aa match; 
  # Ix represents insertions in x (gaps in y); 
  # Iy represents insertions in y (gaps in x); 
  # The traceback now points to the matrix (M, Ix, Iy) from which the
  # maximum was obtained: $FROMM=1, $FROMIX=2, $FROMIY=3
  # B[$dir][1] is the traceback for M; 
  # B[$dir][2] is the traceback for Ix; 
  # B[$dir][3] is the traceback for Iy


  #check if the inputs have the correct dimensionality
  if ( (length($x)!= $#{$posOffsetx}+1) or (length($y)!= $#{$posOffsety}+1) ){
      die "Position Score Alignment error, wrong dimensionality";
  }
  # convert DNA to RNA and convert to uppercase
  $x=formatSeq($x);
  $y=formatSeq($y);
  
  # reverse the SRNA so it can bind to the target
  $x= reverse $x;


  my (@M, @Ix, @Iy, @B);
  for (my $j=0; $j<=$m; $j++) {
    $M[$j] = [(0) x ($n+1)];
    $Ix[$j] = [($minusInf) x ($n+1)];
    $Iy[$j] = [($minusInf) x ($n+1)];
  }
  # The traceback matrix; also correctly initializes borders
  foreach my $dir (@DIR) {
    for (my $j=0; $j<=$m; $j++) {
      for (my $k=1; $k<=3; $k++) {
        $B[$dir][$k][$j] = [($STOP) x ($n+1)];
      }
    }
  }
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      my $s =  DBioConst::nucScore($matrix, substr($x, $i-1, 1), substr($y, $j-1, 1));
      $s+= ${$posOffsetx}[$i-1]+ ${$posOffsety}[$j-1];
      my $val =  DBioConst::max(0, 
                     $M[$j-1][$i-1]+$s, 
                     $Ix[$j-1][$i-1]+$s, 
                     $Iy[$j-1][$i-1]+$s);
      $M[$j][$i] = $val;
      if ($val == $M[$j-1][$i-1]+$s) {
        $B[$FROMM][1][$j][$i] = $RED; 
      } 
      if ($val == $Ix[$j-1][$i-1]+$s) {
        $B[$FROMIX][1][$j][$i] = $RED; 
      } 
      if ($val == $Iy[$j-1][$i-1]+$s) {
        $B[$FROMIY][1][$j][$i] = $RED; 
      } 
      $val = DBioConst::max($M[$j][$i-1]-$d, $Ix[$j][$i-1]-$e, $Iy[$j][$i-1]-$d);  
      $Ix[$j][$i] = $val;
      if ($val == $M[$j][$i-1]-$d) {
        $B[$FROMM][2][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j][$i-1]-$e) {
        $B[$FROMIX][2][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j][$i-1]-$d) {
        $B[$FROMIY][2][$j][$i] = $RED;
      }      
      $val = DBioConst::max($M[$j-1][$i]-$d, $Iy[$j-1][$i]-$e, $Ix[$j-1][$i]-$d);  
      $Iy[$j][$i] = $val;
      if ($val == $M[$j-1][$i]-$d) {
        $B[$FROMM][3][$j][$i] = $RED;
      } 
      if ($val == $Iy[$j-1][$i]-$e) {
        $B[$FROMIY][3][$j][$i] = $RED;
      } 
      if ($val == $Ix[$j-1][$i]-$d) {
        $B[$FROMIX][3][$j][$i] = $RED;
      }      
    }
  }
  # Find maximal score in matrices
  my $vmax = 0;
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      $vmax = DBioConst::max($vmax, $M[$j][$i], $Ix[$j][$i], $Iy[$j][$i]);
    }
  }  
  my ($kmax, $jmax, $imax) = (0, 0);
  for (my $i=1; $i<=$n; $i++) {
    for (my $j=1; $j<=$m; $j++) {
      if ($vmax == $M[$j][$i]) {
        &markReachable3(\@B, 1, $i, $j);
        $kmax = 1;
        $jmax = $j;
        $imax = $i;
      }
      if ($vmax == $Ix[$j][$i]) {
        &markReachable3(\@B, 2, $i, $j);
        $kmax = 2;
        $jmax = $j;
        $imax = $i;
      }
      if ($vmax == $Iy[$j][$i]) {
        &markReachable3(\@B, 3, $i, $j);
        $kmax = 3;
        $jmax = $j;
        $imax = $i;
      }
    }
  }  
  #return (\@M, \@Ix, \@Iy, \@B, &traceback3($x, $y, \@B, $kmax, $imax, $jmax));
  return (\@M,$vmax,&traceback3($x, $y, \@B, $kmax, $imax, $jmax));
}


# ------------------------------------------------------------
# Common subroutines for affine gap cost alignment
# Reconstruct the alignment from the traceback, backwards, 
# and mark green the path actually taken

sub traceback3 {
  my ($x, $y, $B, $k, $i, $j) = @_;   # B by reference
  my ($xAlign, $yAlign) = ("", "");
  while ($$B[$FROMM][$k][$j][$i] != 0 
         || $$B[$FROMIX][$k][$j][$i] != 0 
         || $$B[$FROMIY][$k][$j][$i] != 0) {
    my $nextk;
    # Colour green the path that was actually taken
    if ($$B[$FROMIY][$k][$j][$i]) {
      $$B[$FROMIY][$k][$j][$i] = $GREEN;
      $nextk = 3;       # From Iy
    } elsif ($$B[$FROMIX][$k][$j][$i]) {
      $$B[$FROMIX][$k][$j][$i] = $GREEN;
      $nextk = 2;       # From Ix
    } elsif ($$B[$FROMM][$k][$j][$i]) {
      $$B[$FROMM][$k][$j][$i] = $GREEN;
      $nextk = 1;       # From M
    } 
    if ($k == 1) {              # We're in the M matrix
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $i--; $j--;
    } elsif ($k == 2) {         # We're in the Ix matrix
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-"; 
      $i--;
    } elsif ($k == 3) {         # We're in the Iy matrix
      $xAlign .= "-"; 
      $yAlign .= substr($y, $j-1, 1);
      $j--;
    }       
    $k = $nextk;
  }
  # Warning: these expressions cannot be inlined in the list
  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign, $yAlign);
}


# Mark blue all (affine) traceback arrows reachable from ($k, $i, $j)

sub markReachable3 {
  my ($B, $k, $i, $j) = @_;             # B by reference
  foreach my $dir (@DIR) {
    if ($$B[$dir][$k][$j][$i] == $RED) {
      $$B[$dir][$k][$j][$i] = $BLUE;
      if ($k == 1) {                    # We're in the M matrix
        &markReachable3($B, $dir, $i-1, $j-1);
      } elsif ($k == 2) {               # We're in the Ix matrix
        &markReachable3($B, $dir, $i-1, $j);
      } elsif ($k == 3) {               # We're in the Iy matrix
        &markReachable3($B, $dir, $i,   $j-1);
      }
    }
  }
}








# convert RNA to DNA and convert to uppercase
sub formatSeq{
    my($seq)=@_;

    $seq=~s/\s//g;
    $seq=~s/a/A/g;
    $seq=~s/g/G/g;
    $seq=~s/t/T/g;
    $seq=~s/c/C/g;
    $seq=~s/u/U/g;
    $seq=~s/T/U/g;

    my $L=length($seq);

    if(!($seq=~/(A|G|U|C){$L}/)){
	die "Error in Position Score Alignment: input sequences are not nucleotide sequences";
    }
    return $seq;

}


#----------------------------------------------------------------------
# For debugging only
# Print a given matrix (array of array of number)
# The matrix must be passed by reference as in printmatrix(\@blosum50)

sub printmatrix {
  my ($title, $matrix) = @_;
  print "-" x 70, "\n$title:\n";
  foreach my $row (@$matrix) {
    foreach my $score (@$row) {
      if ($score == $minusInf) {
        $score = "-Inf";
        printf("%4s", $score);
      } else {
        printf("%4d",$score); 
      }
    }
    print "\n";
  }
  print "-" x 70, "\n";
}

1;


