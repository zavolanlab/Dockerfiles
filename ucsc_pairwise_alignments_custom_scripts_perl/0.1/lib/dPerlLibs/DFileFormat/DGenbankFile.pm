use strict;
use Bio::Perl;
use DTable::DTableOps;
use DTable::DOverlappingSequences;
package DGenbankFile;


sub new {
    my $class = shift;
    my $self = {};
    bless $self;
    if (defined $_[0]) {
	$self->{filename} = shift;
    }

    
    my $in  = Bio::SeqIO->new(-file => $self->{filename},-format => 'GenBank');
    $self->{seq} = $in->next_seq();



    return $self;
}

sub getSeq{
    my $self = shift;
    return $self->{seq}->subseq(1,$self->{seq}->length());
}

sub getSubseq{
    my $self = shift;
    my $from = shift; # number to add to our number	
    my $to = shift;

    if ( ($to<$from) or ($from<1) or ($to > $self->{seq}->length() )){
	return "";
    }
    return $self->{seq}->subseq($from,$to);
}

sub getSeqLength{
    my $self = shift;
    return $self->{seq}->length();
}

sub getHeader{
    my $self = shift;

    my $header=$self->{seq}->id()."|".$self->{seq}->desc();
    return $header;
}

sub getGenes{
    my $self = shift;
    my @tab;

    my @feats = $self->{seq}->get_SeqFeatures();
    
    foreach my $feat (@feats) {
	# check if we have a coding sequence "CDS" tag
	if ($feat->primary_tag =~  /gene/){ 
	    # extract the information needed
	    my $start=$feat->start;
	    my $end=$feat->end;
	    my $strand=$feat->strand;
	    $strand=~s/\-1/\-/g;
	    $strand=~s/1/\+/g;

	    my $geneName="";

	    foreach my $tag ($feat->get_all_tags) {
		if ($tag =~ /gene/){
		    my @gl=$feat->get_tag_values($tag);
		    $geneName=$gl[0];
		}
	    }
	    push(@tab,[$start,$end,$strand,$geneName]);
        }
    }
    
    return @tab;
}

sub getCDS{
    my $self = shift;
    
   
    # print "Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n"; 
    my @feats = $self->{seq}->get_SeqFeatures();
    my @tab;
    # extract the coding regions with their translation
    foreach my $feat (@feats) {
	# check if we have a coding sequence "CDS" tag
	if ($feat->primary_tag =~  /CDS/){ 
	    #print "primary tag: ", $feat->primary_tag, "\n";
	    # extract the information needed
	    my $primKey=$feat->primary_tag;
	    my $start=$feat->start;
	    my $end=$feat->end;
	    my $strand=$feat->strand;
	    $strand=~s/\-1/\-/g;
	    $strand=~s/1/\+/g;
	    my $geneName=" ";
	    my $protein="";
	    foreach my $tag ($feat->get_all_tags) {
		if ($tag =~ /translation/){
		    my @proteins=$feat->get_tag_values($tag);
		    $protein=$proteins[0];
		}
		if ($tag =~ /gene/){
		    my @gl=$feat->get_tag_values($tag);
		    $geneName=$gl[0];
	        }
	    }

	    # combine line;
	    if (length($protein) != 0){
		my @line=($start,$end,$strand,$geneName,$protein);
		push(@tab,\@line);
	    }
       }
    }

    return @tab;
}

sub getMisc_RNA{
    my $self = shift;
    my @feats = $self->{seq}->get_SeqFeatures();
    my @tab;
    # extract the misc_RNAs
    foreach my $feat (@feats) {
	if ($feat->primary_tag =~  /misc_RNA/){ 
	    # extract the information needed
	    my $start=$feat->start;
	    my $end=$feat->end;
	    my $strand=$feat->strand;
	    $strand=~s/\-1/\-/g;
	    $strand=~s/1/\+/g;
	    my $geneName=" ";
	    foreach my $tag ($feat->get_all_tags) {
		if ($tag =~ /gene/){
		    my @gl=$feat->get_tag_values($tag);
		    $geneName=$gl[0];
	        }
	    }

	    # combine line
	    my @line=($start,$end,$strand,$geneName);
	    push(@tab,\@line);
       }
    }

    return @tab;
}




sub getRegions{
     my $self = shift;
     my $tabRef = shift;
     my @outTab;

     my @tab=@{$tabRef};
     foreach (@tab){
	 
	 my $cds=$self->{seq}->subseq(${$_}[0],${$_}[1]);
	 push(@outTab,$cds);
     }

     return @outTab;
}

# from 1 to N
sub getOrientedRegion{
    my $self = shift;
    my $start = shift;
    my $end = shift;
    my $strand = shift;
    
     my $cds=$self->{seq}->subseq($start,$end);
     if ($strand eq "-"){
	 $cds=~tr/ACGT/TGCA/;
	 $cds=reverse $cds;
     }
    return $cds;
}


# input table has three columns: start end strand(+/-)
# from 1 to N
sub getOrientedRegions{
     my $self = shift;
     my $tabRef = shift;
     my @outTab;

     my @tab=@{$tabRef};
     foreach (@tab){
	 
	 my $cds=$self->{seq}->subseq(${$_}[0],${$_}[1]);
	 if (${$_}[2] eq "-"){
	    $cds=~tr/ACGT/TGCA/;
	    $cds=reverse $cds;
	 }

	 push(@outTab,$cds);
     }

     return @outTab;
}




# returns all regions that are surounded by genes
# it also returns the transcription direction of the 
# neighbouring genes

sub getIgRegions{
    my $self = shift;
    
    my @genes=$self->getGenes();
  
     # find clusters for overlapping sequences 
    my($totalCodingRegionsRef,$codingRegionsInfoRef)=DOverlappingSequences::mapSeqsIntoClusters(\@genes,
						    \&DOverlappingSequences::doRegionsOverlap);


    my @totalCodingRegions=@{$totalCodingRegionsRef};
  

    # compute the complement of the transcribed regions = intergenic regions
    my @toCoReg=DTableOps::selectColumns(\@totalCodingRegions,[0,1]); #select only the first two columns
    my @intergenicRegionsIndices=DOverlappingSequences::complementRegions(\@toCoReg);

    # extract the transcription direction of the neighbouring genes
    my @intGenInfo;
    for my $i(0.. $#intergenicRegionsIndices){
	push(@intGenInfo,[ $intergenicRegionsIndices[$i][0], 
			   $intergenicRegionsIndices[$i][1],
			   $genes[$totalCodingRegions[$i][-1]][2],
			   $genes[$totalCodingRegions[$i+1][2]][2]
			   ]) ;

    }

    #DTableOps::printTable(\@totalCodingRegions);
    #DTableOps::printTable(\@intergenicRegionsIndices); 
    #DTableOps::printTable(\@intGenInfo); 
    #DTableOps::printTable(\@genes);

    return @intGenInfo;
}


sub getVerifiedSigma70Promoters{

    my $self = shift;

    my @feats = $self->{seq}->get_SeqFeatures();
    my @tab;
    # extract the experimentally verified promoters
    foreach my $feat (@feats) {
	# check if we have a coding sequence "CDS" tag
	if ($feat->primary_tag =~  /promoter/){ 
	    # extract the information needed
	    my $start=$feat->start;
	    my $end=$feat->end;
	    my $strand=$feat->strand;
	    $strand=~s/\-1/\-/g;
	    $strand=~s/1/\+/g;
	    my $note=" ";
	    my $seq="";
	    foreach my $tag ($feat->get_all_tags) {
		if ($tag =~ /note/){
		    my @notes=$feat->get_tag_values($tag);
		    $note=$notes[0];
		}
	    }

	    if ($note=~/Sigma70.+documented/){
		$seq=$self->{seq}->subseq($start,$end);
		#reverse complement if on - strand
		if ($strand eq "-"){
		    $seq=~tr/ACGT/TGCA/;
		    $seq=reverse $seq;
		}

		#print $note,"\n";
		# combine line;
		my @line=($start,$end,$strand,$note,$seq);
		push(@tab,\@line);
	    }
       }
    }

    return @tab;


}


sub loadFile{
    my ($filename)=@_;
    
    my $in  = Bio::SeqIO->new(-file => $filename,-format => 'GenBank');

    my $seq = $in->next_seq(); 
    # print "Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n"; 
    my @feats = $seq->get_SeqFeatures();
    my @tab;
    # extract the coding regions with their translation
    foreach my $feat (@feats) {
	# check if we have a coding sequence "CDS" tag
	if ($feat->primary_tag =~  /CDS/){ 
	    #print "primary tag: ", $feat->primary_tag, "\n";
	    # extract the information needed
	    my $primKey=$feat->primary_tag;
	    my $start=$feat->start;
	    my $end=$feat->end;
	    my $statusLine="";
	    my $protein="";
	    foreach my $tag ($feat->get_all_tags) {
		if ($tag =~ /translation/){
		    my @proteins=$feat->get_tag_values($tag);
		    $protein=$proteins[0];
		}else{
		    #print "  tag: ", $tag, "\n";
		    $statusLine.=$tag."=";
		    foreach my $value ($feat->get_tag_values($tag)) {
			$statusLine.=$value;
		        #print "    value: ", $value, "\n";
		    }
		    $statusLine.= "|";
	        }
	    }

	    # combine line;
	    if (length($protein) != 0){
		my @line=($seq->id,$start,$end,$statusLine,$protein);
		push(@tab,\@line);
	    }
       }
    }

    return @tab;

}



# test stuff

sub test{

    my $in  = Bio::SeqIO->new(-file => "ML.fa",-format => 'GenBank');

    while ( my $seq = $in->next_seq() ) { 
	print "Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n"; 
    } 

}


sub modifyAndSave{

    my $self = shift;
    my ($tabRef)=shift;
    my ($filename)=shift;
    my @tab=@{$tabRef};

    my @feats = $self->{seq}->get_SeqFeatures();

    DTableOps::printTable(\@tab);
    my $counter=0;
    foreach my $feat (@feats) {
       
	if ($feat->primary_tag =~  /gene/){ 
	    # extract the information needed
	    my($start,$end,$strand)=($tab[$counter][0],$tab[$counter][1],$tab[$counter][2]);
	    $feat->start($start);
	    $feat->end($end);
	    $feat->strand($strand);
	    $counter++;
        }

	if ($feat->primary_tag =~  /CDS/){ 
	    # extract the information needed
	    my($start,$end,$strand)=($tab[$counter-1][0],$tab[$counter-1][1],$tab[$counter-1][2]);
	    $feat->start($start);
	    $feat->end($end);
	    $feat->strand($strand);
        }


    }

    my $out = Bio::SeqIO->new(-file => ">$filename" , '-format' => 'GenBank');
    
    $out->write_seq($self->{seq});
}



1;
