use strict;
use DHash::DNestedHashOps;
use DString::DStringOps;
use DComputation::DMath;
package DRNAhybridFileFormat;


sub loadFile{
    my ($filename)=@_;
    my @outTab;

    open(INFILE,$filename);
    my $line= <INFILE>;

    while ($line){
	if ($line=~/target: (.+)/){
	    my $targetName=$1;

	    $line= <INFILE>;
	    $line=~/length: (\d+)/;
	    my $targetLength=$1;

	    $line= <INFILE>;
	    $line=~/miRNA : (.+)/;
	    my $mirnaName=$1;
	
	    $line= <INFILE>;
	    $line=~/length: (\d+)/;
	    my $mirnaLength=$1;
	
	    $line= <INFILE>;
	    $line= <INFILE>;
	    $line=~/mfe: (\-?\d+\.?\d*|\-?\.\d+)/;
	    my $mfe=$1;

	    $line= <INFILE>;
	    $line=~/p-value: (.+)/;
	    my $p=$1;

	    $line= <INFILE>;
	    $line= <INFILE>;
	    $line=~/position  (\d+)/;
	    my $position=$1;

	    $line= <INFILE>;
	    $line=~/target 5\' (.+) 3\'/;
	    my $aln1=$1;

	    $line= <INFILE>;
	    $line=~/          (.+)/;
	    my $aln2=$1;

	    $line= <INFILE>;
	    $line=~/          (.+)/;
	    my $aln3=$1;

	    $line= <INFILE>;
	    $line=~/miRNA  3\' (.+) 5\'/;
	    my $aln4=$1;

	    my $note="";
	    $line= <INFILE>;
	    if ($line=~/note: (.+)/){
		$note=$1;
	    }

	    $line= <INFILE>;
	    $line= <INFILE>;

	    # add missing spaces to $aln2 and $aln3
	    my $divSp=DMath::max(length($aln1),length($aln4))-length($aln2);
	    $aln2=$aln2.DStringOps::charChain(" ",$divSp);
	    $aln3=$aln3.DStringOps::charChain(" ",$divSp);
  

	    #print " $targetName $targetLength $mirnaName $mirnaLength $mfe $p $position \n";
	    #print "$aln1\n$aln2\n$aln3\n$aln4\n";

	    push(@outTab,[$targetName,$targetLength,$mirnaName,$mirnaLength,$mfe,$p,$position,$aln1,$aln2,$aln3,$aln4,$note]);
	}
    }
    close(INFILE);
    return @outTab;
}

sub loadFileLite{
    my ($filename)=@_;
    my @outTab;

    open(INFILE,$filename);
    my $line= <INFILE>;

    while ($line){
	if ($line=~/target: (.+)/){
	    my $targetName=$1;

	    $line= <INFILE>;
	    $line=~/length: (\d+)/;
	    my $targetLength=$1;

	    $line= <INFILE>;
	    $line=~/miRNA : (.+)/;
	    my $mirnaName=$1;
	
	    $line= <INFILE>;
	    $line=~/length: (\d+)/;
	    my $mirnaLength=$1;
	
	    $line= <INFILE>;
	    $line= <INFILE>;
	    $line=~/mfe: (\-?\d+\.?\d*|\-?\.\d+)/;
	    my $mfe=$1;

	    $line= <INFILE>;
	    $line=~/p-value: (.+)/;
	    my $p=$1;

	    $line= <INFILE>;
	    $line= <INFILE>;
	    $line=~/position  (\d+)/;
	    my $position=$1;

	    $line= <INFILE>;
	    $line=~/target 5\' (.+) 3\'/;
	    my $aln1=$1;

	    $line= <INFILE>;
	    $line=~/          (.+)/;
	    my $aln2=$1;

	    $line= <INFILE>;
	    $line=~/          (.+)/;
	    my $aln3=$1;

	    $line= <INFILE>;
	    $line=~/miRNA  3\' (.+) 5\'/;
	    my $aln4=$1;

	    my $note="";
	    $line= <INFILE>;
	    if ($line=~/note: (.+)/){
		$note=$1;
	    }

	    $line= <INFILE>;
	    $line= <INFILE>;

	    # add missing spaces to $aln2 and $aln3
	    my $divSp=DMath::max(length($aln1),length($aln4))-length($aln2);
	    $aln2=$aln2.DStringOps::charChain(" ",$divSp);
	    $aln3=$aln3.DStringOps::charChain(" ",$divSp);
  

	    #print " $targetName $targetLength $mirnaName $mirnaLength $mfe $p $position \n";
	    #print "$aln1\n$aln2\n$aln3\n$aln4\n";

	    push(@outTab,[$targetName,$targetLength,$mirnaName,$mirnaLength,$mfe,$p,$position,$note]);
	}
    }
    return @outTab;
}


sub loadIntoNestedHash{
    my ($filename,$option)=@_;

    my %nestedHash;

    my @tab=loadFileLite($filename);

    #DTableOps::printTable(\@tab);

    for my $i(0..$#tab){
	my $targetName=$tab[$i][0];
	my $targetLenght=$tab[$i][1]; 
	my $mirnaName=$tab[$i][2];
	my $mirnaLength=$tab[$i][3]; 
	my $mfe=$tab[$i][4];
	my $p=$tab[$i][5]; 
	my $position=$tab[$i][6];
	my $noteIn=$tab[$i][7];
	
	my $hspID="";
	if ($noteIn=~/HspID=([^\;]+)/){$hspID=$1;}	

	# set top element in tree structure
	if ($option==1){
	    $targetName=$tab[$i][2];
	    $mirnaName=$tab[$i][0];
	}

	#print " $targetName $targetLenght $mirnaName $mirnaLength $mfe $p $position \n";
	#print "$aln1\n$aln2\n$aln3\n$aln4\n";
	my $note="pos=$position;$noteIn";

	# insert an entry into the nested hash

	# check if the target already exists
	if (exists  $nestedHash{$targetName }){
	    # check if the query exists
	    my $query=${$nestedHash{$targetName}}[1];

	    if (exists ${$query}{$mirnaName}){
		# insert new hsp key   
		my $hsp= ${${$query}{$mirnaName}}[1] ;
		${$hsp}{$hspID}=[$note,undef];

		#print %{$hsp},"\n";
	    }else{
		# insert new query key
		my %hspHash;
		$hspHash{$hspID}=[$note,undef];
		${$query}{$mirnaName}=["",\%hspHash]; 

	    }
	    #print $targetName,"\n"; 
	    #DNestedHashOps::printNHT(\%query);

	}else{
	    # target does not exist. insert it"
	    my %hspHash;
	    $hspHash{$hspID}=[$note,undef];

	    my %queryHash;
	    $queryHash{$mirnaName}=["",\%hspHash];
	    $nestedHash{$targetName}=["",\%queryHash];

	    #DNestedHashOps::printNHT(\%nestedHash);

	    #die;
	}
    }
    return %nestedHash;
}


sub loadIntoNestedHashWithALN{
    my ($filename,$option)=@_;

    my %nestedHash;

    my @tab=loadFile($filename);

    #DTableOps::printTable(\@tab);

    for my $i(0..$#tab){
	my $targetName=$tab[$i][0];
	my $targetLenght=$tab[$i][1]; 
	my $mirnaName=$tab[$i][2];
	my $mirnaLength=$tab[$i][3]; 
	my $mfe=$tab[$i][4];
	my $p=$tab[$i][5]; 
	my $position=$tab[$i][6];
	my $aln1=$tab[$i][7];
	my $aln2=$tab[$i][8];
	my $aln3=$tab[$i][9];
	my $aln4=$tab[$i][10];
	my $noteIn=$tab[$i][11];
	
	my $hspID="";
	if ($noteIn=~/HspID=([^\;]+)/){$hspID=$1;}	

	my $aln="mRNA  5' $aln1 3'\t         $aln2\t         $aln3\tmiRNA 3' $aln4 5'\t";
	#$aln="GAGA";

	# set top element in tree structure
	if ($option==1){
	    $targetName=$tab[$i][2];
	    $mirnaName=$tab[$i][0];
	}

	#print " $targetName $targetLenght $mirnaName $mirnaLength $mfe $p $position \n";
	#print "$aln1\n$aln2\n$aln3\n$aln4\n";
	my $note="pos=$position;$noteIn"."HYB=$aln;";

	# insert an entry into the nested hash

	# check if the target already exists
	if (exists  $nestedHash{$targetName }){
	    # check if the query exists
	    my $query=${$nestedHash{$targetName}}[1];

	    if (exists ${$query}{$mirnaName}){
		# insert new hsp key   
		my $hsp= ${${$query}{$mirnaName}}[1] ;
		${$hsp}{$hspID}=[$note,undef];

		#print %{$hsp},"\n";
	    }else{
		# insert new query key
		my %hspHash;
		$hspHash{$hspID}=[$note,undef];
		${$query}{$mirnaName}=["",\%hspHash]; 

	    }
	    #print $targetName,"\n"; 
	    #DNestedHashOps::printNHT(\%query);

	}else{
	    # target does not exist. insert it"
	    my %hspHash;
	    $hspHash{$hspID}=[$note,undef];

	    my %queryHash;
	    $queryHash{$mirnaName}=["",\%hspHash];
	    $nestedHash{$targetName}=["",\%queryHash];

	    #DNestedHashOps::printNHT(\%nestedHash);

	    #die;
	}
    }
    return %nestedHash;
}



1;
