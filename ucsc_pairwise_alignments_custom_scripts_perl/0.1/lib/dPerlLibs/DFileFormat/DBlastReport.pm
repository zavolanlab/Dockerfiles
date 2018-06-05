use strict;
use Bio::Perl;
use Bio::SearchIO;

package DBlastReport;



sub simpleBlastSummary{
    my ($blast_report)=@_;
    my @outTab;
 
    my $searchio = new Bio::SearchIO (-format => 'blast',-file   => $blast_report);
    my $in = new Bio::SearchIO(-format => 'blast', -file   => "$blast_report");
    while( my $result = $in->next_result ) {
	#print  $result->query_name,"\n";
	my @line=($result->query_name);

	while( my $hit = $result->next_hit ) {
	    while( my $hsp = $hit->next_hsp ) {
		# iterate trough all hits
		#my @hitnamel=split(/\|/,$hit->name);
		push(@line,$hit->name);

		#print $hit->name,"\n";
		#print $hit->description, "\n";
		#print $hsp->score(),"\n";
		#print $hsp->length('total'), "\n";
		#print $hsp->percent_identity, "\n";			   
	    }  
	}
	push(@outTab,\@line);
    }
    return @outTab;
}

1;
