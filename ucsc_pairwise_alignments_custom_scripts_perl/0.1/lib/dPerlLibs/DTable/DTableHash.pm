use strict;
package DTableHash;


# laods a file into a hash of hash
#my $tt=DTableHash::loadFile($filename,'^\#','\s+');

sub loadFile{
    my ($filename,$headerRG,$delim)=@_;
   
    my %tableHash;
    open(F,$filename) or die "Cannot open $filename \n";
    my $header=<F>;
    if($header =~/^$headerRG\s*(\S+.+)/){
	my @colIDs=split($delim,$1);
	shift(@colIDs); # remove the first id

	while(<F>){
	    if (/\S/){
		my@line=split($delim,$_);
		if (scalar(@line) == scalar(@colIDs)+1){
		    my $rowId=shift(@line);
		    if (defined $tableHash{$rowId}){
			die "the table in $filename has a not unique identifier $rowId\n";
		    }
		    for (my$i=0;$i<@colIDs;$i++){
			$tableHash{$rowId}{$colIDs[$i]}=$line[$i];
		    }
		}else{die "the table in $filename contains different number of elements in different rows\n";}
	    }

	}

	close(F);
	return \%tableHash;


    }else{die "$filename does not contain a header which starts with $headerRG\n";}
   
}

# input: 
# matrix and the annotations for cols and rows
sub convertMatrixToTableHash{
    my($matrix,$rowIDs,$colAttr)=@_;

    my %tableHash;
    # some checks, the dimensions have to fit
    for (my $i=0;$i<@{$matrix};$i++){
	for (my $j=0;$j<@{$matrix->[$i]};$j++){
	    $tableHash{$rowIDs->[$i]}{$colAttr->[$j]}=$matrix->[$i][$j];
	}
    }

    return \%tableHash;
}


sub print{
    (my $tableHash)=@_;

    # print the header
    print "#id\t";
    my @hList=sort keys %{$tableHash->{(sort keys %{$tableHash})[0]}};
    print join("\t",@hList),"\n";

    for my $key (sort keys %{$tableHash}){
	print "$key";
	for my $att(sort keys %{$tableHash->{$key}}){ 
	    print "\t$tableHash->{$key}{$att}";

	}
	print "\n";

    }


}



1;
