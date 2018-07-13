#!/usr/bin/perl -w
use strict;
my %h=();


open(IN,$ARGV[1]) || die;
my %HYBRIDS;
while(<IN>){
    chomp;
    if($_=~/^>/){
	#>24204|chr15|+|39379971|39380009|5865|5691	<<< VERSUS >>>	>hsa-miR-143	----> MIRZA Target Quality= 3.73981e-07
	#miRNA             5' 	GAG-----A-UG---A-A----GCACUG---UAG---------------------- 3'	GAGAUGAAGCACUGUAG
	#A L I G N M E N T	|O|vvvvv|v||vvv|v|vvvv|^^^^|vvv|||######################
	#mRNA   (reversed) 3' 	TATAAAGATAACGTGTATAGGGC----CAAAATTGAAAAGTT-------------- 5'	ATATAAAGATAACGTGTATAGGGCCAAAATTGAAAAGTT

	$_=~s/\(//g;
	$_=~s/\>//g;
	$_=~s/\)//g;
	my @R=split(/\s+/);
	my @e=split(/[^0-9a-zA-Z\_\.]/,$R[0]);
	my $mirna = <IN>;
	my @d=split(/\s+/,$mirna);
	$HYBRIDS{$R[0]}{$R[3]}{"miRNA_ALIGNMENT"}= join(" ",$d[1],$d[2],$d[3])."<br>";
	my $alignment =<IN>;
	@d=split(/\s+/,$alignment);
#	$HYBRIDS{$e[0]}{$R[3]}{"ALIGNMENT"} ="-->".$d[9]."<br>";
	$HYBRIDS{$R[0]}{$R[3]}{"ALIGNMENT"} = $d[9];
	my $mrna = <IN>;
	@d=split(/\s+/,$mrna);
	$HYBRIDS{$R[0]}{$R[3]}{"mRNA_ALIGNMENT"}= join(" ",$d[2],$d[3],$d[4])."<br>";
    }
}

open(H2, ">./results/best_results.html");
print H2 "<HTML><head><script src=\"../bin/sorttable.js\"></script></head><BODY><h1>miRNA Target Prediction Table (sortable)<br></h1>\n";
print H2 "<table class=\"sortable\">\n";

print H2 "<thead><tr><th><u>TargetID</u></th><th><u>miRNA</u></th><th><u>miRNA Sequence</u></th><th><u>miRNA Activity</u></th><th><u>MIRZA Target Quality</u></th><th>Hybridization</th></tr></thead>\n";

print H2 "<tbody>\n";

open(IN,$ARGV[0]) || die;
while(<IN>){
    #(>NM_170695|1658|1697|3|Ago2_MNase_CLIP43.fa,Ago2_365nm_CLIP11.fa,Ago2_254nm_CLIP31.fa|98.0472731522276)	TTGCATTATTTTATATATTTTTTATTAATATTTGCACAT	(>hsa-miR-19a)	UGUGCAAAUCUAUGCAAAACU	11407.3	695.34	0.0609558
    #(>16745|chr21|+|33574041|33574079|60|53)        ATAACTCATGGACTTATAATGTGCAATACTGGAAAAACG (>hsa-miR-32)   AUUGCACAUUACUAAGU       120546  7416.3  0.0615228
    chomp;
    $_=~s/\(//g;
    $_=~s/\>//g;
    $_=~s/\)//g;
    my @d=split(/\s+/);
    my @e=split(/[^0-9a-zA-Z\_\.]/,$d[0]);
    $h{$d[0]}++;

    if($h{$d[0]} <=1){
	print STDERR "ALL\t$d[0]\t".$d[4]."\t".$d[5]."\n";	
	print H2 "<tr><td>".$d[0]."</td><td><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=".$d[2]."\"><center>".$d[2]."</center></a></td><td><center>".$d[3]."</center></td><td><center>".$d[6]."</center></td><td><center>".$d[4]."</center></td><td><font face=\"Courier\">".$HYBRIDS{$d[0]}{$d[2]}{"miRNA_ALIGNMENT"}."</font><font face=\"Courier\" color=\"red\">-->".$HYBRIDS{$d[0]}{$d[2]}{"ALIGNMENT"}."<br></font><font face=\"Courier\">".$HYBRIDS{$d[0]}{$d[2]}{"mRNA_ALIGNMENT"}."</font><br></td></tr>\n";

	my $aln = $HYBRIDS{$d[0]}{$d[2]}{"ALIGNMENT"};
	my $mRNA_aln = $HYBRIDS{$d[0]}{$d[2]}{"mRNA_ALIGNMENT"};
	my @mrna_data =split(/\s+/,$mRNA_aln);
	$mRNA_aln =$mrna_data[1];
	@mrna_data =split(//,$mRNA_aln);
	
        my $miRNA_aln = $HYBRIDS{$d[0]}{$d[2]}{"miRNA_ALIGNMENT"};
        my @mirna_data =split(/\s+/,$miRNA_aln);
        $miRNA_aln =$mirna_data[1];
	@mirna_data =split(//,$miRNA_aln);

	my @data =split(//,$aln);

	my $flag="N";
	my $begin_mirna=0;
	for(my $i=0;$i<@mirna_data;$i++){
		if($mirna_data[$i] eq "A" ||$mirna_data[$i] eq "C" ||$mirna_data[$i] eq "G"||$mirna_data[$i] eq "U"){
			if($flag eq "N"){
				$flag = "Y";
			}	
		}
		elsif($flag eq "N"){	
			$begin_mirna++;
		}
	}
	$flag="N";
        my $end_mirna=@mirna_data;
	for(my $i=@mirna_data-1;$i>=0;$i--){
                if($mirna_data[$i] eq "A" ||$mirna_data[$i] eq "C" ||$mirna_data[$i] eq "G"||$mirna_data[$i] eq "U"){
                        if($flag eq "N"){
                                $flag = "Y";
                        }
                }
                elsif($flag eq "N"){
                        $end_mirna--;
                }
        }
	$flag="N";
	my $begin_mrna=0;
        for(my $i=0;$i<@mrna_data;$i++){
                if($mrna_data[$i] eq "A" ||$mrna_data[$i] eq "C" ||$mrna_data[$i] eq "G"||$mrna_data[$i] eq "U"){
                        if($flag eq "N"){
                                $flag = "Y";
                        }
                }
                elsif($flag eq "N"){
                        $begin_mrna++;
                }
        }
        $flag="N";
        my $end_mrna=@mrna_data;
        for(my $i=@mrna_data-1;$i>=0;$i--){
                if($mrna_data[$i] eq "A" ||$mrna_data[$i] eq "C" ||$mrna_data[$i] eq "G"||$mrna_data[$i] eq "U"){
                        if($flag eq "N"){
                                $flag = "Y";
                        }
                }
                elsif($flag eq "N"){
                        $end_mrna--;
                }
        }
	my $check=substr($aln,$begin_mirna,8);
        my $chk1=substr($check,1,7);
        my $chk2=substr($check,1,6);

	my $chk3=substr($aln,$begin_mirna,$end_mirna-$begin_mirna);
	my $chk4=substr($miRNA_aln,$begin_mirna,$end_mirna-$begin_mirna);

	my $chk5=substr($aln,$begin_mrna,$end_mrna-$begin_mrna);
        my $chk6=substr($mRNA_aln,$begin_mirna,$end_mirna-$begin_mirna);

	my @mrna=split(//,$chk6);

	if($chk1 eq "|||||||" || ($chk2 eq "||||||" && $mrna[0] eq 'A')){
	    print STDERR join("\t","CANONICAL",$d[4],$d[5],$chk3,$chk4,$chk5,$chk6,$d[0],$d[2],$e[0],$aln)."\n";
	}
	else{
	   print STDERR join("\t","NON_CANONICAL",$d[4],$d[5],$chk3,$chk4,$chk5,$chk6,$d[0],$d[2], $e[0], $aln)."\n";
	}
    }
}



