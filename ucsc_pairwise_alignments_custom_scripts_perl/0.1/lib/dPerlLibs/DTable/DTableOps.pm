use strict;
package DTableOps;


# loads file and returns an array of strings
# 2 parameters: 
# 1: filename
# 2: regexp for the line to be included '' for no matching
# usage: 
# my @lines=DTableOps::loadFile("input.txt",'\.\.');
# my $targets=DTableOps::convertLinesToTableRef(DTableOps::loadFileRef($ARGV[1],''),'\s+');

sub loadFile{
    my ($filename,$rg)=@_;

    my @fileContent;
    
    open(FHANDLE,$filename) or die "Cannot open $filename \n";
    while(<FHANDLE>){
	if(/$rg/){
	    # delete newline char at the end of the line
	    s/\n$//;
	    push(@fileContent, $_);
	}
    }
    close(FHANDLE);

    return @fileContent;
}

sub loadFileRef{
    my ($filename,$rg)=@_;

    my @fileContent;
    
    open(FHANDLE,$filename) or die "Cannot open $filename \n";
    while(<FHANDLE>){
	if(/$rg/){
	    # delete newline char at the end of the line
	    s/\n$//;
	    push(@fileContent, $_);
	}
    }
    close(FHANDLE);
 
    return \@fileContent;
}


# save a table to a file
# 4 parameters:
# 1. table to save
# 2. delimiter for each columns e.g. "|"
# 3. delimiter for each line e.g. "\n"
# 4. filename
# usage:
# DTableOps::saveTable(\@table,"|","\n","output.txt");
sub saveTable{
    my($refTab,$delimLine,$delimCol,$filename)=@_;
    my @tab=@{$refTab};
    my $i;
    my $j;

    open(FHANDLE, ">".$filename);
    for $i ( 0 .. $#tab ) {
	 for $j ( 0 .. $#{ $tab[$i] } ) {
	     print FHANDLE  $tab[$i][$j],$delimLine;
	 }
	 print FHANDLE $delimCol;
    }
    close(FHANDLE);
}




# creates a table from an array of lines
# 2 parameters:
# 1. array of lines from which a table should be built
# 2. delimiter for line dl
# usage:
# my @tab=DTableOps::convertLinesToTable(\@lines,'\.\.'); !!! single quotes for the regexp
sub convertLinesToTable{
    my ($refLines,$dl)=@_;
    my @lines=@{$refLines};
    my @tab;
    my $counter=0;
    foreach my$line (@lines){ 
	# delete newline char at the end of the line
	$line =~ s/\n$//;
	# remove spaces at the beginning
	$line =~ s/^\s+//;
        my @sp=split(/$dl/,$line);  # my: allocate memory each time in the loop
	if (@sp != 0){
	    $tab[$counter]=\@sp;
	    $counter++;
	}
    }

    return @tab;
}

sub convertLinesToTableRef{
    my ($refLines,$dl)=@_;
    my @tab;
    my $counter=0;
    foreach (@{$refLines}){ 
	# delete newline char at the end of the line
	s/\n$//;
	# remove spaces at the beginning
	s/^\s+//;
        my @sp=split(/$dl/,$_);  # my: allocate memory each time in the loop
	if (@sp != 0){
	    $tab[$counter]=\@sp;
	    $counter++;
	}
    }

    return \@tab;
}



# select columns from a table (deepcopy)
# 2 parameters
# 1. table from which the columns will be selected
# 2. which columns to select: list of indices (starting from zero)
# !! its up to the user to check if all the selected columns are filled with data
# usage: my @codingRegionsIndices=DTableOps::selectColumns(\@codingRegionsTbl,[0,1]);
sub selectColumns{
     my($refTab,$indicesRef)=@_;
     my @tab=@{$refTab};
     my @indices=@{$indicesRef};
     my $i;
     my @outTab;
     
     for $i ( 0 .. $#tab ) {
	 #select the columns
	 
	 my @tempRow=@{$tab[$i]}[@indices];
	 $outTab[$i]=\@tempRow;
     }
    
     return @outTab;
}

# select rows from a tablw
# 2 paramters
# 1. table from which the rows will be selected
# 2. array with the indices of the rows to be selected
#    eg [3,2,17,37]

sub selectRows{
    my($tab,$indices)=@_;
    
    my @outTab;
    foreach (@{$indices}){

	push(@outTab,$tab->[$_]);

    }
    return @outTab;
        
}




# select rows from a table
# 3 parameters
# 1. table from which the rows will be selected
# 2. select the column to be filtered for
# 3. match
# usage
# my @outTab=selectRowsWhere(\@tab,0,"+");

sub selectRowsWhere{
    my($refTab,$col,$match)=@_;
    my @tab=@{$refTab};

    my @outTab;
    for my $i(0..$#tab){
	if ($tab[$i][$col] eq $match){
	    push(@outTab,$tab[$i]);
	}
    }
    return @outTab;   
}

sub selectRowsWhichMatchSet{
    my($refTab,$col,$set)=@_;

    my @outTab;
    for (my $i=0;$i<@{$refTab};$i++){
	if (${$set}{$refTab->[$i][$col]}){
	    push(@outTab,$refTab->[$i]);
	}

    }


    return \@outTab;   
}



sub selectRowsWhereValSmallerThan{
    my($refTab,$col,$thresh)=@_;
    my @tab=@{$refTab};

    my @outTab;
    for my $i(0..$#tab){
	if ($tab[$i][$col] < $thresh){
	    push(@outTab,$tab[$i]);
	}
    }
    return @outTab;   
}

sub selectRowsWhereValGreaterThan{
    my($refTab,$col,$thresh)=@_;
    my @tab=@{$refTab};

    my @outTab;
    for my $i(0..$#tab){
	if ($tab[$i][$col] > $thresh){
	    push(@outTab,$tab[$i]);
	}
    }
    return @outTab;   
}


#concat two tables one underneath the other (deepcopy)
# Table 2 after Table 1
# 2 parameters:
# 1. Table 1
# 2. Table 2
# usage: my @concatRegionIndices=DTableOps::verticalConcatTables(\@codingRegionsIndices,\@structRNAindices);  
sub verticalConcatTables{
     my($tab1,$tab2)=@_;
     my @outTab;
     my $i;
     my $counter=0;

     for $i ( 0 .. $#{$tab1} ) {
	 #select the columns
	 $outTab[$counter]=$tab1->[$i];
	 $counter++;
     }

     for $i ( 0 .. $#{$tab2} ) {
	 #select the columns
	 $outTab[$counter]=$tab2->[$i];
	 $counter++;
     }
     return @outTab;
}

# sort a table based on the selected column
# 2 input parameters
# 1. table to sort
# 2. number of column
# 
# 2 output parameters
# 1. sorted table 
# 2. list of indices
#
sub sortTableForColumn{
     my($refTab,$colToSort)=@_;
     my @tab=@{$refTab};
     my $i;
     my @listToSort;
     my @outTable;

     # extract the selected column

     for $i ( 0 .. $#tab ) {	
	 push(@listToSort,$tab[$i][$colToSort]);
     }

     #foreach (@listToSort){print $_,"\n"}

     my @sortedIndices = sort { $listToSort[$a] <=> $listToSort[$b] } 0..$#listToSort;
     my @sortedValues= @listToSort[@sortedIndices];

     # foreach (@sortedValues){print $_,"\n"}

     # rearrange the complete table based on sortedIndices
     
     my $counter=0;
     foreach (@sortedIndices){
	 my @tempLine=@{$tab[$_]};
	 $outTable[$counter]=\@tempLine;
	 $counter++;
     }

     return (\@outTable,\@sortedIndices);
}


# prints the table to the standard output
# 1 parameter: table to print
# usage:
# DTableOps::printTable(\@tab);
sub printTable{
    my ($refTab)=@_;

    print "\nPrinting table...\n----------------\n";
    for my $i(0..$#{$refTab}){
	print $i,". ";
	for my $j(0..$#{$refTab->[$i]}){
	    if (defined $refTab->[$i][$j]){
		print $refTab->[$i][$j], "||";
	    }else{
		print "undef", "||";
	    }
	}
	print "\n";
    }
    print "----------------\n\n";  
}


# same as  printTable but formated for floats and ints
sub printMatrix{
    my ($refTab)=@_;
    my @tab=@{$refTab};
    my $i;
    my $j;
    print "\nPrinting table...\n----------------\n";
 for $i ( 0 .. $#tab ) {
     print $i,". ";
     for $j ( 0 .. $#{ $tab[$i] } ) {
         printf "%4.2f ", $tab[$i][$j];
     }
     print "\n";
 }
 print "----------------\n\n";  
}

# print all entries with their coordinates:
sub printTableCoord{
    my ($refTab)=@_;
    my @tab=@{$refTab};
    my $i;
    my $j;
    print "Printing table\n";
 for $i ( 0 .. $#tab ) {
     for $j ( 0 .. $#{ $tab[$i] } ) {
         print "elt $i $j is $tab[$i][$j]\n";
     }
 }  
}

sub printTableClean{
    my ($refTab)=@_;
    my @tab=@{$refTab};
    my $i;
    my $j;
 for $i ( 0 .. $#tab ) {
     for $j ( 0 .. $#{ $tab[$i] } ) {
         print $tab[$i][$j], "  ";
     }
     print "\n";
 }
}



#some helpersubs
sub removeNLineAtEndOfLine{
    my ($line)=@_;
    $line =~ s/\n$//;
    return $line;
}

# saves a matrix as a bitmap file for visualization purpose
# values have to be beween 0 and 1. otherwise use the $factor
# to scale the image
sub saveMatrixAsPGM{
    my($tab,$filename,$factor,$pixPerCell)=@_;

    my$height=$pixPerCell*($#{$tab}+1);
    my$width=$pixPerCell*($#{$tab->[0]}+1);

    open(FHANDLE, ">".$filename);
    print FHANDLE"P2\n$width $height\n255\n"; # this is needed for a PGM file
    for my $i(0..$#{$tab}){
	my $line;
	for my $j(0..$#{$tab->[$i]}){
	    my$val=$tab->[$i][$j];
	    my$valScl=DMath::max(DMath::min($val*$factor*256,255),0);
	    for (my $pix=0;$pix<$pixPerCell;$pix++){
		$line.=int($valScl)." ";
	    }
	    #print FHANDLE int($valScl)," ";
	}
	for (my $pix=0;$pix<$pixPerCell;$pix++){
	    print FHANDLE $line." ";
	}

	#print FHANDLE " ";
    }


    close(FHANDLE); 
}

# test stuff

# create rectangular table with 1s in the diagonal and return it
# 1 input parameter: size of the matrix d (d=m=n)
# the returned table has the dimensionality d*d
sub diagonal{
    my ($d)=@_;
    my @tab;
   
    # create table filled with zeros:
    for my $j (0..$d-1){
	for my $i (0..$d-1){
	    $tab[$j][$i]=0.0;
        }
    }
    # fill in ones in the diagonal
    for my $i (0..$d-1){
	$tab[$i][$i]=1.0;
    }
	
    return @tab
}

# create a table and return it to the caller
# no input parameter
sub favTable{
    
    # my @a=("asdfsdl","sdlfkdj");
    # return @a;

    my @tabInt = ( [11, 13], [14, 15, 17], [0] );
    return @tabInt;

}

# look how to create arrays
sub arrayCheck{
    my @arr;
    $arr[0]="tony";
    $arr[1]="bob";
    $arr[2]="heinz";

    print @arr,"\n";
}

sub tableInfo{
    my ($refTab)=@_;
    my @tab=@{$refTab};

    print $#tab,"\n";
    print  $#{ $tab[2]} ,"\n";
   
    $tab[0][0]=18;
    

}

sub passTable{
    my ($refTab)=@_;

     print ${$refTab}[1][2], "\n";
     #print @{$refTab},"\n";

}


sub pass2Arrays{
    my ($refArr1,$refArr2)=@_;
    my @Arr1=@{$refArr1};
    my @Arr2=@{$refArr2};

    print @{$refArr1},"\n";
    print @{$refArr2},"\n";

    print @Arr1,"\n";
    print @Arr2, "\n";
}



1;
