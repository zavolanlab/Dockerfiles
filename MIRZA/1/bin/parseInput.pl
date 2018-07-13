#!/usr/bin/perl -w
#use strict; 
#use warnings;

my %mirnaseq = ();
my %mirnaexp = ();
my $buffer;


open MIRNASEQ, "<", $ARGV[0] or die $!;
open MIRNAEXP, "<", $ARGV[1] or die $!;

while(<MIRNAEXP>) {
  #read lines in the expression file, delimiter is a single space
  chomp($_);
  my $delimiter = " ";
  my @mirnaexp_line = split($delimiter, $_);
  # and put them in an associative array
  $mirnaexp{$mirnaexp_line[0]} = $mirnaexp_line[1];
}

while (my $mirnaseq_line = <MIRNASEQ>) {
  # read lines in the sequence file
  if ($. % 2) {
    chomp $mirnaseq_line;
    #remove the first character of the name, if it is ">" 
    if ($mirnaseq_line =~ m/^>/) {
      $mirnaseq_line =~ s/^.//;
    }
    $buffer = $mirnaseq_line;
  }
  else {
    chomp $mirnaseq_line;
    #put the lines in an associative array
    if (length($mirnaseq_line)==21) {
      $mirnaseq{$buffer} = $mirnaseq_line;
    }
  }
}

#check which entries the two associatives array have in common

open OUTPUTSEQ, ">", $ARGV[2] or die $!;
open OUTPUTEXP, ">", $ARGV[3] or die $!;

foreach $nameexp(keys %mirnaexp) {
  if(exists($mirnaexp{$nameexp}) && exists($mirnaseq{$nameexp})){
    #and write two new files
    print OUTPUTSEQ ">$nameexp\n$mirnaseq{$nameexp}\n";
    print OUTPUTEXP "$nameexp\t$mirnaexp{$nameexp}\n";
  }
}

close (OUTPUTSEQ);
close (OUTPUTEXP);
close (MIRNASEQ);
close (MIRNAEXP);


