# implementation of a table without any restrictions for its columns
# it should be used to read data from files
# works for all input files that are delimited by two diffentent characters

use strict;
package DDelimTable;


#public functions

# constructor
#loads a text file into an unformated table-form.
#1 parameter:
#filename: file to load

sub new {
        my $class = shift;
        my $self = {};
        bless $self;
        if (defined $_[0]) {
                $self->{filename} = shift;
        }
	#try to open the file
	




        return $self;
}




#formats the data based on the chosen delimiter
#2 parameters: 
#delimCol: the delimiter of the column
#delimLine: the delimiter of the line

sub formatDelimTable{

}
1;
