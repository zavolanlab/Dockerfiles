use strict;
use DComputation::DMath;
package DSplitFileFormat;


# these functions split a given file into many small ones
# they try to split according to the given number
# $com  :command
# $path :path of the input files
# $mainSub : place where the small files should be copied in




sub splitRNAhybFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;


    # parse the filename(s)
    my $target="";
    my $query="";
    if($com=~/\-t (\S+)/){
	$target=$1;
    }
    if($com=~/\-q (\S+)/){
	$query=$1;
    }
    if (($target eq "") or( $query eq "")){die "RNAhybrid params cannot be parsed\n";}

    # split the target
    my $outFiles=generalSplit("$path/$target",$pieces,'>');
    generalSavePieces($outFiles,$mainSub,$target);


    # copy the query
    generalCopyFileToAllTaskDirs($outFiles,$path,$mainSub,$query);

    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}





sub splitRNAhybFilterFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;

    my $sRegExp='^target: ';

    # parse the filename(s)
    my $filename;
    if($com=~/filterRNAhybReport\.pl (\S+)/){
	$filename=$1;
    }

    my $outFiles=generalSplit("$path/$filename",$pieces,$sRegExp);
    generalSavePieces($outFiles,$mainSub,$filename);

    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}

sub splitRNAlocalFoldProbFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;

    
    # parse the filename(s)
    my $mrnas="";
    my $hsps="";

    if ($com=~/rnaLocalFoldProb\.pl\s+(\S+)\s+(\S+)/){
	$mrnas=$1;
	$hsps=$2;
    }else{die "cannot split the input of rnaLocalFoldProb\n"}

    # split the rnaHybrid report
    my $outFiles=generalSplit("$path/$hsps",$pieces,'^target: ');
    generalSavePieces($outFiles,$mainSub,$hsps);
    
    # copy symbolic links for the mrnas to the task directories
    # this avoids copying these large file for every process

    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,$mrnas);

    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}

sub splitBlastpFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;

    
    # parse the filename(s)
    my $target="";
    my $query="";

    if ($com=~/blastp\s+(\S+)\s+(\S+)/){
	$target=$1;
	$query=$2;

    }else{die "cannot split the input of blastp\n"}

    # split the rnaHybrid report
    my $outFiles=generalSplit("$path/$query",$pieces,'>');
    generalSavePieces($outFiles,$mainSub,$query);
    
    # copy symbolic links for the mrnas to the task directories
    # this avoids copying these large file for every process

    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,$target);
    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,"$target.xpd");
    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,"$target.xps");
    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,"$target.xpt");

    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}


sub splitSpaFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;


    # parse the filename(s)
    my $inputFa="";
    if($com=~/(\S+.fa)/){
	$inputFa=$1;
    }

    if (($inputFa eq "")){die "Spa params cannot be parsed\n";}

    # split the target
    my $outFiles=generalSplit("$path/$inputFa",$pieces,'>');
    generalSavePieces($outFiles,$mainSub,$inputFa);


    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}


sub splitTailHybrConsFileFormat{
    my ($com,$pieces,$path,$mainSub)=@_;

    
    # parse the filename(s)
    my $tr="";
    my $alns="";
    my $nMerOcc="";
    my $miRNAs="";

    if ($com=~/tailHybrCons.pl\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	$tr=$1;
	$alns=$2;
	$nMerOcc=$3;
	$miRNAs=$4;

    }else{die "cannot split the input of tailHybrCons\n"}

    # split the miRNA file
    my $outFiles=generalSplit("$path/$miRNAs",$pieces,'>');
    generalSavePieces($outFiles,$mainSub,$miRNAs);
    
    # copy symbolic links for the mrnas to the task directories
    # this avoids copying these large file for every process

    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,$tr);
    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,$alns);
    generalCreateSymlinksToAllTaskDirs($outFiles,$path,$mainSub,$nMerOcc);

    my $realNumberOfTasks=$#{$outFiles}+1;
    return $realNumberOfTasks;
}




sub generalSplit{
    my($filename,$pieces,$regExp)=@_;


    my @outFiles;
    for my $i(1..$pieces){
	$outFiles[$i-1]="";
    }

    open(FHANDLE,$filename) or die "Cannot open $filename \n";
  
    my $roundInd=-1;
    my $maxInd=0;
    while(<FHANDLE>){
	if ($_=~/$regExp/){
	    $roundInd+=1;
	    $roundInd %= $pieces;
	    $maxInd=DMath::max($maxInd,$roundInd);
	}
	if ($roundInd != -1){
	    $outFiles[$roundInd].=$_;
	}
    }
    close(FHANDLE);

    # check if any of the pieces are empty -> delete
    for (my $j=$#outFiles; $j>$maxInd; $j--){
	pop(@outFiles);
    }

    return \@outFiles;
}


# save the splitted files into a different directories (taskX)
sub generalSavePieces{
    my($outFiles,$mainSub,$filename)=@_;
   
    for my $i(1..$#{$outFiles}+1){
	system("mkdir $mainSub/task.$i");
	open(FHANDLE, ">$mainSub/task.$i/$filename");
	print FHANDLE "${$outFiles}[$i-1]";
	close(FHANDLE);
    } 
}


sub generalCopyFileToAllTaskDirs{
    my($outFiles,$path,$mainSub,$filename)=@_;
    
    for my $i(1..$#{$outFiles}+1){
	system("cp $path/$filename $mainSub/task.$i/$filename");
    }
}


sub generalCreateSymlinksToAllTaskDirs{
    my($outFiles,$path,$mainSub,$filename)=@_;
    
    for my $i(1..$#{$outFiles}+1){
	system("ln -s $path/$filename $mainSub/task.$i/$filename");
    }
}




1;
