use strict;
package DFeatureFormat;


# saves a feature file from a table
# start|end|strand|note

sub saveFeatFromTbl{
    my ($tabRef,$featureName,$filename)=@_;
    my @tab=@{$tabRef};

    open(FHANDLE, ">".$filename);
    
    for my $i(0..$#tab){
	my $start=$tab[$i][0];
	my $end=$tab[$i][1];
	my $strand=$tab[$i][2]; # (+/-) !!
	my $note=$tab[$i][3];
	my $hLine="     $featureName          ";

        print FHANDLE $hLine;
	if ($strand eq "+"){
	    print FHANDLE "$start..$end\n";
	}else{
	    print FHANDLE "complement($start..$end)\n";
	}
	print FHANDLE spaces(length($hLine)),"/note=\"",$note,"\"\n";


    }

    close(FHANDLE);


}

sub loadFeatureFile{
    my($filename)=@_;
    my @outTab;

    open(INFILE,"$filename");
    my $line= <INFILE>;
    while($line){
	if (($line=~/\s+\w+\s+\d+\.\.\d+/) or ($line=~/\s+\w+\s+complement\(\d+\.\.\d+\)/) ) {
	    my $start;
	    my $end;
	    my $strand;
	    if ($line=~/\s+\w+\s+complement\((\d+)\.\.(\d+)\)/){
		$start=$1;
		$end=$2;
		$strand="-";
		#print "$1 $2 - \n";
            }
	    if ($line=~/\s+\w+\s+(\d+)\.\.(\d+)/){
		$start=$1;
		$end=$2;
		$strand="+";
		#print "$1 $2 + \n";
            }
	    if ($line=<INFILE>){}else{$line="";}
	    my $note="";
	    while($line=~/\/\w+\=\"(.*)\"/){
		$note.=$1;
		if ($line=<INFILE>){}else{$line="";}
	    }
	    push(@outTab,[$start,$end,$strand,$note]);
	}else{  
	    if ($line=<INFILE>){}else{$line="";}
	}
   }
   close(INFILE);

    return @outTab;


}



sub spaces{
    my ($nrOfSpaces)=@_;	
    my $spacesLine=""; 
    for my $i(0..$nrOfSpaces-1){
	$spacesLine.=" ";
    }
    return $spacesLine;
}

1;
