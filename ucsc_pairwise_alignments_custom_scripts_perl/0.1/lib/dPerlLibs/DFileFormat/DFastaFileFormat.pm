use strict;
package DFastaFileFormat;

# loads a fasta file with one or more sequences in it
# 1 parameter: filename
#
# returns table with 2 columns:
# header|sequence
# usage: my @tab=DFastaFileFormat::loadFile($ARGV[0]);

sub loadFile{
    my ($filename)=@_;
    open(FHANDLE,$filename) or die "Cannot open $filename \n";
  
    my $nrOfSeq=0;
    my @headers;
    my @sequences;
    my @tabOut;

    while(<FHANDLE>){
	if((/\>/)){
	    # if there is a > char then a new sequence starts
	    s/\n//g;
	    s/^\>//g;
	    $nrOfSeq++;
	    $tabOut[$nrOfSeq-1][0]=$_;
	}else{
	    if ($nrOfSeq>0){
		s/\n//g;
		#tr/acgtu/ACGTU/;
		$tabOut[$nrOfSeq-1][1].=$_;
	    }
        }

    }

    close(FHANDLE);
    return (@tabOut);


}

# same as loadFile but with References
sub loadFileFast{
    my ($filename)=@_;
    open(FHANDLE,$filename) or die "Cannot open $filename \n";
  
    my $nrOfSeq=0;
    my @headers;
    my @sequences;
    my @tabOut;

    while(<FHANDLE>){
	if((/\>/)){
	    # if there is a > char then a new sequence starts
	    s/\n//g;
	    s/^\>//g;
	    $nrOfSeq++;
	    $tabOut[$nrOfSeq-1][0]=$_;
	}else{
	    if ($nrOfSeq>0){
		s/\n//g;
		#tr/acgtu/ACGTU/;
		$tabOut[$nrOfSeq-1][1].=$_;
	    }
        }

    }

    close(FHANDLE);
    return (\@tabOut);
}

# loads only the first sequence from afasta file
sub loadFirstSequence{
    my ($filename)=@_;
    # read genome fasta file with only one sequence in it and return the sequence string:
    
    my @tab=loadFile($filename);

    return ($tab[0][0],$tab[0][1]);

}

sub saveSingleSeq{
    my ($sequence,$header,$filename)=@_;
     
    open(FHANDLE, ">".$filename);
    
    print FHANDLE ">".$header,"\n",insertNLintoSeq($sequence,70);
    close(FHANDLE);
}

sub saveSeqsFromTBL{
    my ($tabRef,$filename)=@_;
    my @tab=@{$tabRef};

    open(FHANDLE, ">".$filename);
    
    for my $i(0..$#tab){
	my $header=$tab[$i][0];
	my $sequence=$tab[$i][1];

        print FHANDLE ">".$header,"\n",insertNLintoSeq($sequence,70);
    }

    close(FHANDLE);
}


# converts a fasta file table with many sequneces to a hash
# 2 inputs
# 1. table with sequences my @tab=DFastaFileFormat::loadFile($filename)
# 2. regexp which extracts the id from the header of each sequence
# usage:
# DFastaFileFormat::convertSeqsToHash(\@tab,'.+\|\d+\.\.\d+\|');
# ... for an id like "NM_152899|1783..1798|"

sub convertSeqsToHash{
    my ($tab,$idRegExp)=@_;
    my %outHash;
 
    # convert table to a hash
    # ids can be in different forms:

    for my $i(0..$#{$tab}){
	my $header=$tab->[$i][0];
	my $seq=$tab->[$i][1];

	if ($header=~/($idRegExp)/){
	    if (defined $outHash{$1}){die "sequence $1 is multiple times in the file\n";}

	    $outHash{$1}=$seq;
	}else{
	    die "Error in DFastaFileFormat. the fasta input $i($header) contains a".
		"header that cannot be converted into an Id";
	}
    }
  
    #foreach(keys %outHash){print "$_\n";}
  
    return \%outHash;

}







# creates a fasta file from a feature table and a sequences array
# input1: feature table start|end|strand|note
# input2: seqence array 
sub saveFromFeatAndSeqs{
    my($feat,$seqs,$filename)=@_;

    my @tab;
    for my $i(0..$#{$feat}){
	$tab[$i][0]="$feat->[$i][0]..$feat->[$i][1]|$feat->[$i][2]|$feat->[$i][3]";
	$tab[$i][1]=$seqs->[$i];
    }

    saveSeqsFromTBL(\@tab,$filename);

}




# saves a hash containing the headers and the sequences to a file (filename)
sub saveFileFromHash{
    my ($sequencesRef,$filename)=@_;
    my %sequences=%{$sequencesRef};
    
    open(FHANDLE, ">".$filename);

    foreach (keys %sequences) {
	print FHANDLE ">".$_,"\n",insertNLintoSeq($sequences{$_},70);
    }
    close(FHANDLE);
}
1;


sub insertNLintoSeq{
    my ($seq,$width)=@_;
    my $outSeq="";

    for (my $count=0; $count<length($seq); $count+=$width){
	$outSeq.=substr($seq,$count,$width)."\n";
    }
    return $outSeq;
}
