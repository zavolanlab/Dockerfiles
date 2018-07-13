INSTALLATIONFOLDER=$1
MIRNAEXPRESSIONFILE=$2
MIRNASEQUENCEFILE=$3
MIRNAEXPOUTPUTFILE=$4
MIRNASEQOUTPUTFILE=$5

# cut everything that is longer than 21 in the miRNA sequence file
cut -c -21 $MIRNASEQUENCEFILE >> $INSTALLATIONFOLDER/tmpmirnaseq.txt

# parse the miRNA expression file and replace "\t" or "   " with " "
awk -v OFS=" " '$1=$1' $MIRNAEXPRESSIONFILE >> $INSTALLATIONFOLDER/tmpmirnaexp.txt

# call perl script for parsing miRNA files
$INSTALLATIONFOLDER/bin/parseInput.pl $INSTALLATIONFOLDER/tmpmirnaseq.txt  $INSTALLATIONFOLDER/tmpmirnaexp.txt $MIRNASEQOUTPUTFILE $MIRNAEXPOUTPUTFILE 

#remove tmp files
rm $INSTALLATIONFOLDER/tmpmirnaseq.txt
rm $INSTALLATIONFOLDER/tmpmirnaexp.txt

