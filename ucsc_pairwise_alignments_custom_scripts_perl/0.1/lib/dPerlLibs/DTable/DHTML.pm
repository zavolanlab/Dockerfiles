use strict;
package DHTML;

sub printRegTable{
    my($tab,$caption,$colNames)=@_;

    # print the first line of the table with the field descriptions
    printHeader($caption);
    print "<TR>";
    for my $c(0..$#{$colNames}){
	print "<TD>${$colNames}[$c]";
    }
    print "\n";

    #print all the table
    for my $i(0..$#{$tab}){
	print "<TR>";
	for my $j(0..$#{$tab->[$i]}){
	    print "<TD>$tab->[$i][$j]";
	}
	print "\n";
    } 


    printFooter();
}


sub printHeader{
    my ($caption)=@_;
    #print "<FONT face=\"Courier\">\n";
    print "<TABLE border=\"1\"\n";
    print "       summary=\"TABLE\">\n";
    print "<CAPTION><EM>$caption</EM></CAPTION>\n";
}
sub printFooter{
    print "</TABLE>\n";
    print "</FONT>\n";
}
1;
